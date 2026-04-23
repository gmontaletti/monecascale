# moneca_sbm(): hierarchical degree-corrected SBM backend
# ========================================================
#
# Alternative clustering backend for MONECA based on hierarchical
# degree-corrected stochastic block models (DC-SBM). The theoretical bridge
# (MONECA's RR = O/E is the DC-SBM residual under identity block interaction)
# makes SBM inference a principled drop-in for the clique-enumeration step.
#
# Two backends are supported:
#   * "greed" (R-native, via greed::DcSbm); scales to ~10^4-10^5 nodes.
#   * "graphtool" (Python via reticulate); scales to 10^6+ nodes.

#' Hierarchical DC-SBM Backend for MONECA Segmentation
#'
#' Alternative to \code{moneca::moneca_fast()} that replaces clique
#' enumeration with hierarchical degree-corrected stochastic block model
#' (DC-SBM) inference. MONECA's weight matrix \code{RR = O / E} is the
#' expected-value residual under a DC-SBM with identity block interaction,
#' so DC-SBM clustering is a principled scalable replacement for the clique
#' step.
#'
#' The number of hierarchical levels is selected automatically by the
#' backend (ICL for \code{greed}; MDL for \code{graphtool}), which
#' eliminates MONECA's manual \code{segment.levels} requirement. Passing
#' \code{segment.levels} truncates the backend's full hierarchy to that
#' many approximately log-spaced levels, preserving the finest and coarsest.
#'
#' @param mx A square mobility matrix, optionally with margins as the last
#'   row and column (auto-detected via \code{has_margins = "auto"}). Sparse
#'   \code{dgCMatrix} input is supported end-to-end through the
#'   \code{greed} backend.
#' @param backend One of \code{"auto"}, \code{"greed"}, \code{"graphtool"}.
#'   \code{"auto"} prefers \code{graphtool} when available and the matrix
#'   exceeds 5000 nodes; otherwise uses \code{greed}.
#' @param deg_corr Logical. If \code{TRUE} (default), fits a degree-corrected
#'   SBM. Set to \code{FALSE} for plain SBM.
#' @param edge_model One of \code{"poisson"} (default, appropriate for
#'   count-valued mobility) or \code{"bernoulli"} (binary edges).
#' @param small.cell.reduction Numeric. Cells in the core matrix below this
#'   value are zeroed before SBM inference. Mirrors the
#'   \code{moneca::moneca_fast()} parameter.
#' @param symmetric_method One of \code{"sum"} (default; \code{mx + t(mx)}),
#'   \code{"min"} (\code{2 * pmin(mx, t(mx))}, down-weights one-way flows),
#'   or \code{"none"} (keep directed). Applied to the core before SBM.
#' @param segment.levels Integer, \code{NULL}, or the string \code{"auto"}.
#'   If \code{NULL} (default), all hierarchy levels produced by the backend
#'   (from K_final down to 2) are returned. If an integer, the hierarchy is
#'   truncated to that many approximately log-spaced levels. If
#'   \code{"auto"}, the full hierarchy is fit and then
#'   \code{\link{auto_segment_levels}} picks a single preferred level via
#'   \code{auto_method}; the object is trimmed to that level and carries
#'   the picker result under \code{$auto_level}.
#' @param auto_method One of \code{"mdl"} (default) or \code{"mi_plateau"}.
#'   Ignored unless \code{segment.levels = "auto"}.
#' @param max_K Integer or \code{NULL}. Upper bound on number of blocks at
#'   the finest level (greed backend only). Defaults to
#'   \code{max(20, n_core / 5)}.
#' @param isolates Logical. If \code{TRUE}, attaches an
#'   \code{$isolates_summary} slot mirroring \code{moneca::moneca_fast()}.
#' @param has_margins One of \code{"auto"} (default), \code{TRUE}, or
#'   \code{FALSE}. Controls whether to detect and/or append row/column
#'   margins to the input matrix.
#' @param seed Integer seed for reproducible inference.
#' @param n_init Integer. Number of independent backend runs; the
#'   best-scoring fit is retained. Defaults to 1 for \code{greed} (which is
#'   quasi-deterministic) and 5 is typical for \code{graphtool} (MCMC).
#' @param verbose Logical. If \code{TRUE}, prints backend progress.
#' @param ... Additional arguments forwarded to the backend.
#'
#' @return An object of class \code{"moneca"} with the same structure as
#'   \code{moneca::moneca_fast()}, plus an extra slot
#'   \code{$sbm_diagnostics} describing the backend used, the number of
#'   blocks per level, and the MDL/ICL trace.
#'
#' @seealso \code{moneca::moneca_fast()} for the clique-based reference
#'   implementation; \code{\link{moneca_sbm_install_graphtool}} to set up
#'   the \code{graphtool} backend.
#'
#' @examples
#' \dontrun{
#' mob <- moneca::generate_mobility_data(
#'   n_classes = 20, n_total = 2000, seed = 42
#' )
#' seg <- moneca_sbm(mob, backend = "greed", seed = 1L)
#' moneca::segment.membership(seg)
#' moneca::plot_moneca_hierarchical(seg)
#' }
#'
#' @export
moneca_sbm <- function(
  mx,
  backend = c("auto", "greed", "graphtool"),
  deg_corr = TRUE,
  edge_model = c("poisson", "bernoulli"),
  small.cell.reduction = 0,
  symmetric_method = c("sum", "min", "none"),
  segment.levels = NULL,
  auto_method = c("mdl", "mi_plateau"),
  max_K = NULL,
  isolates = FALSE,
  has_margins = "auto",
  seed = NULL,
  n_init = NULL,
  verbose = FALSE,
  ...
) {
  backend <- match.arg(backend)
  edge_model <- match.arg(edge_model)
  symmetric_method <- match.arg(symmetric_method)
  auto_method <- match.arg(auto_method)

  # 0. Auto-level gate ------------------------------------------------------
  auto_mode <- identical(segment.levels, "auto")
  segment.levels_for_fit <- if (auto_mode) NULL else segment.levels

  # 1. Margin handling ------------------------------------------------------
  mx_info <- .ensure_margins_sbm(mx, has_margins = has_margins)
  mx_full <- mx_info$mx
  margins_added <- mx_info$margins_added
  n_core <- nrow(mx_full) - 1L

  if (n_core < 3L) {
    stop("Need at least 3 non-margin rows/cols.", call. = FALSE)
  }

  # 2. Backend selection ----------------------------------------------------
  resolved_backend <- .resolve_sbm_backend(backend, n_core)
  if (is.null(n_init)) {
    n_init <- if (resolved_backend == "graphtool") 5L else 1L
  }

  # 3. Preprocess core for SBM ---------------------------------------------
  core_processed <- .preprocess_for_sbm(
    mx_full,
    small.cell.reduction = small.cell.reduction,
    symmetric_method = symmetric_method,
    zero_diagonal = TRUE
  )

  core_names <- rownames(mx_full)[1:n_core]

  # 4. Backend fit ----------------------------------------------------------
  sbm_fit <- switch(
    resolved_backend,
    "greed" = .sbm_fit_greed(
      core_mx = core_processed,
      deg_corr = deg_corr,
      edge_model = edge_model,
      seed = seed,
      n_init = n_init,
      max_K = max_K,
      node_names = core_names,
      verbose = verbose,
      ...
    ),
    "graphtool" = .sbm_fit_graphtool(
      core_mx = core_processed,
      deg_corr = deg_corr,
      edge_model = edge_model,
      seed = seed,
      n_init = n_init,
      node_names = core_names,
      verbose = verbose,
      ...
    )
  )

  # 5. Adapter to moneca S3 object -----------------------------------------
  out <- .sbm_fit_to_moneca(
    sbm_fit = sbm_fit,
    mx_full = mx_full,
    small.cell.reduction = small.cell.reduction,
    margins_added = margins_added,
    symmetric_method = symmetric_method,
    segment.levels = segment.levels_for_fit,
    isolates = isolates
  )

  # 6. Post-hoc auto-level selection ---------------------------------------
  if (auto_mode) {
    picked <- auto_segment_levels(out, method = auto_method)
    out <- .trim_moneca_to_level(out, picked$level)
    out$auto_level <- picked
  }

  out
}

# 6. Backend auto-resolution --------------------------------------------------

.resolve_sbm_backend <- function(backend, n_core) {
  if (backend == "greed") {
    if (!requireNamespace("greed", quietly = TRUE)) {
      stop(
        "backend = 'greed' but the 'greed' package is not installed. ",
        "Install with install.packages('greed').",
        call. = FALSE
      )
    }
    return("greed")
  }
  if (backend == "graphtool") {
    .check_graphtool_available(hard_stop = TRUE)
    return("graphtool")
  }
  # auto
  large <- n_core > 5000L
  gt_ok <- .check_graphtool_available(hard_stop = FALSE)
  greed_ok <- requireNamespace("greed", quietly = TRUE)

  if (large && gt_ok) {
    return("graphtool")
  }
  if (greed_ok) {
    if (large && !gt_ok) {
      message(
        "moneca_sbm(): n_core = ",
        n_core,
        " exceeds the greed comfort zone (~5000). ",
        "Consider installing graph-tool via moneca_sbm_install_graphtool() ",
        "for better scaling."
      )
    }
    return("greed")
  }
  if (gt_ok) {
    return("graphtool")
  }
  stop(
    "Neither 'greed' nor 'graphtool' backend is available. ",
    "Install greed with install.packages('greed'), or set up graph-tool ",
    "via moneca_sbm_install_graphtool().",
    call. = FALSE
  )
}

.check_graphtool_available <- function(hard_stop = FALSE) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    if (hard_stop) {
      stop(
        "backend = 'graphtool' requires reticulate. ",
        "Install with install.packages('reticulate').",
        call. = FALSE
      )
    }
    return(FALSE)
  }
  ok <- tryCatch(
    reticulate::py_module_available("graph_tool"),
    error = function(e) FALSE
  )
  if (hard_stop && !ok) {
    stop(
      "backend = 'graphtool' requested but the 'graph_tool' Python module ",
      "is not available in the active reticulate environment. ",
      "Run moneca_sbm_install_graphtool() to set it up.",
      call. = FALSE
    )
  }
  ok
}
