# moneca_bipartite(): hierarchical degree-corrected LBM backend
# ==============================================================
#
# Bipartite (rectangular) sibling of moneca_sbm() for person x employer,
# worker x occupation, and similar two-mode mobility data. Fits a
# hierarchical degree-corrected Latent Block Model (DcLbm) on the
# rectangular relative-risk matrix and packages each view (row, col) as a
# full moneca-class object consumable by moneca's analysis / plotting
# stack. The joint co-clustering diagnostics are retained in a dedicated
# $bipartite_diagnostics slot.
#
# Output class is "moneca_bipartite" (disjoint from "moneca"); dispatch on
# out$rows / out$cols for moneca consumers.

#' Bipartite (Rectangular) DC-LBM Backend for MONECA Segmentation
#'
#' Rectangular counterpart of [moneca_sbm()]. Accepts a two-mode mobility
#' matrix (rows and columns are distinct entity types, e.g. persons x
#' employers, workers x occupations) and fits a hierarchical
#' degree-corrected Latent Block Model (`greed::DcLbm`). The rectangular
#' relative-risk matrix `RR = O / E` acts as the DcLbm residual under
#' identity block interaction, mirroring the square RR = DC-SBM bridge
#' that underpins `moneca_sbm()`.
#'
#' Because `moneca::weight.matrix()` and the downstream `segment.*()` stack
#' assume square inputs, the bipartite fit is exposed to moneca through
#' two linked moneca-class objects built from one-mode projections of the
#' rectangular RR. The row view uses `RR %*% D_c^{-1} %*% t(RR)` and the
#' column view the transpose, with `D_c = diag(col_margin / grand_total)`
#' (analogously for rows). Each projection is zero-diagonal, symmetric,
#' sparse, and carries moneca-style margins appended as the last row and
#' column.
#'
#' @param mx A rectangular mobility matrix of shape `n_row x n_col`.
#'   Optional moneca-style margins in the last row and column are
#'   auto-detected when `has_margins = "auto"`. Sparse `dgCMatrix` input
#'   is supported end to end.
#' @param backend Currently only `"greed"` (DcLbm via the `greed`
#'   package). Argument kept for forward compatibility with future
#'   bipartite backends.
#' @param deg_corr Logical. Reserved. `greed::DcLbm` is always degree
#'   corrected; a non-default value is silently accepted but has no
#'   effect.
#' @param edge_model One of `"poisson"` (default, appropriate for
#'   count-valued mobility) or `"bernoulli"`. Reserved. DcLbm uses the
#'   Poisson edge model; non-default values emit a warning and are
#'   otherwise ignored.
#' @param small.cell.reduction Numeric. Cells in the core rectangle with
#'   values below this threshold are zeroed before RR computation.
#'   Mirrors the `moneca::moneca_fast()` parameter.
#' @param margin_policy One of `"row_col"` (default; append margins to
#'   each one-mode projection, matching moneca's convention) or `"none"`
#'   (skip appending margins). `"row_col"` is the load-bearing choice
#'   for downstream moneca compatibility.
#' @param segment.levels Integer, `NULL`, or the string `"auto"`. If `NULL`
#'   (default), all hierarchy levels produced by the backend (from the
#'   finest joint `K` down to `2`) are returned. If an integer, the
#'   hierarchy is truncated to that many approximately log-spaced levels.
#'   If `"auto"`, the full hierarchy is fit and then
#'   [auto_segment_levels()] picks a preferred level **independently per
#'   side**; `$rows` and `$cols` are trimmed accordingly and each carries
#'   the picker result under `$auto_level`. The two sides may end up at
#'   different levels.
#' @param auto_method One of `"mdl"` (default) or `"mi_plateau"`. Ignored
#'   unless `segment.levels = "auto"`.
#' @param max_K Integer or `NULL`. Reserved for backend-specific caps; the
#'   `greed` backend currently uses its internal default.
#' @param seed Integer seed for reproducible inference.
#' @param n_init Integer. Number of independent backend runs; the
#'   best-scoring fit (highest ICL) is retained. Defaults to `1L`.
#' @param has_margins One of `"auto"` (default), `TRUE`, or `FALSE`.
#'   Controls whether to detect and/or append row/column margins to the
#'   input matrix.
#' @param verbose Logical. If `TRUE`, prints backend progress.
#' @param ... Additional arguments forwarded to the backend.
#'
#' @return An object of class `"moneca_bipartite"` with slots:
#'   * `rows` - full moneca-class object for the row view (one-mode
#'     projection on rows).
#'   * `cols` - full moneca-class object for the column view.
#'   * `bipartite_diagnostics` - list with backend metadata, per-level
#'     row / column block counts, joint ICL / MDL traces, and
#'     per-level `Krow x Kcol` block-interaction matrices computed
#'     directly from the rectangular core.
#'
#' @details
#' The row view and column view each expose the full moneca output
#' contract (`segment.list`, `mat.list`, `segment_metadata`, and
#' companion scalars). `moneca::segment.membership()`,
#' `moneca::segment.quality()`, `moneca::plot_moneca_ggraph()`, and
#' `moneca::plot_moneca_hierarchical()` run unchanged on `out$rows` or
#' `out$cols`.
#'
#' The joint `Krow x Kcol` co-clustering is not representable as a
#' single square moneca object. It is preserved in
#' `bipartite_diagnostics$block_interaction_matrix_per_level`, which
#' cross-tabulates the rectangular core by row and column memberships at
#' each level.
#'
#' @seealso [moneca_sbm()] for the square (one-mode) counterpart;
#'   `moneca::moneca_fast()` for the reference clique-based algorithm on
#'   square matrices.
#'
#' @examples
#' if (requireNamespace("greed", quietly = TRUE)) {
#'   # Toy rectangular mobility fixture
#'   set.seed(1)
#'   rect <- matrix(stats::rpois(40 * 25, lambda = 2), 40, 25)
#'   rownames(rect) <- paste0("R", seq_len(40))
#'   colnames(rect) <- paste0("C", seq_len(25))
#'   out <- moneca_bipartite(rect, seed = 1L, verbose = FALSE)
#'   moneca::segment.membership(out$rows)
#'   moneca::segment.membership(out$cols)
#' }
#'
#' @export
moneca_bipartite <- function(
  mx,
  backend = "greed",
  deg_corr = TRUE,
  edge_model = "poisson",
  small.cell.reduction = 0,
  margin_policy = c("row_col", "none"),
  segment.levels = NULL,
  auto_method = c("mdl", "mi_plateau"),
  max_K = NULL,
  seed = NULL,
  n_init = 1L,
  has_margins = "auto",
  verbose = FALSE,
  ...
) {
  margin_policy <- match.arg(margin_policy)
  auto_method <- match.arg(auto_method)

  # 0. Auto-level gate ------------------------------------------------------
  auto_mode <- identical(segment.levels, "auto")
  segment.levels_for_fit <- if (auto_mode) NULL else segment.levels

  # 1. Reserved-arg sanity --------------------------------------------------
  if (!isTRUE(deg_corr)) {
    warning(
      "deg_corr is reserved: the greed DcLbm backend is always degree ",
      "corrected; the supplied value has no effect.",
      call. = FALSE
    )
  }
  if (!identical(edge_model, "poisson")) {
    warning(
      "edge_model is reserved: the greed DcLbm backend uses the Poisson ",
      "edge model; the supplied value has no effect.",
      call. = FALSE
    )
  }

  # 2. Margin handling ------------------------------------------------------
  rect_info <- .ensure_margins_bipartite(mx, has_margins = has_margins)
  rect_core <- rect_info$rect_core
  n_row <- nrow(rect_core)
  n_col <- ncol(rect_core)

  if (n_row < 3L || n_col < 3L) {
    stop(
      "moneca_bipartite() needs at least 3 non-margin rows and 3 ",
      "non-margin cols; got ",
      n_row,
      " x ",
      n_col,
      ".",
      call. = FALSE
    )
  }

  # 3. Preprocess core ------------------------------------------------------
  rect_core <- .preprocess_for_bipartite(
    rect_core,
    small.cell.reduction = small.cell.reduction
  )
  rect_info$rect_core <- rect_core

  # 4. Backend dispatch -----------------------------------------------------
  resolved_backend <- .resolve_bipartite_backend(backend)

  # 5. Rectangular RR -------------------------------------------------------
  rr <- .bipartite_rr(
    rect_core = rect_core,
    row_margin = rect_info$row_margin,
    col_margin = rect_info$col_margin,
    grand_total = rect_info$grand_total
  )

  # 6. Backend fit ----------------------------------------------------------
  fit <- switch(
    resolved_backend,
    "greed" = .bipartite_fit_greed(
      rr = rr,
      rect_core_sum = sum(rect_core),
      seed = seed,
      n_init = n_init,
      max_K = max_K,
      segment.levels = segment.levels_for_fit,
      verbose = verbose,
      ...
    )
  )

  # 7. Adapter to moneca_bipartite S3 object --------------------------------
  out <- .fit_to_moneca_bipartite(
    fit = fit,
    rect_info = rect_info,
    rr = rr,
    small.cell.reduction = small.cell.reduction,
    margin_policy = margin_policy,
    verbose = verbose
  )

  # 8. Post-hoc per-side auto-level selection -------------------------------
  if (auto_mode) {
    picked_rows <- auto_segment_levels(out$rows, method = auto_method)
    picked_cols <- auto_segment_levels(out$cols, method = auto_method)
    out$rows <- .trim_moneca_to_level(out$rows, picked_rows$level)
    out$cols <- .trim_moneca_to_level(out$cols, picked_cols$level)
    out$rows$auto_level <- picked_rows
    out$cols$auto_level <- picked_cols
  }

  out
}

# 9. Backend resolution -------------------------------------------------------

.resolve_bipartite_backend <- function(backend) {
  if (!is.character(backend) || length(backend) != 1L) {
    stop("backend must be a single string.", call. = FALSE)
  }
  if (!identical(backend, "greed")) {
    stop(
      "moneca_bipartite() currently supports only backend = 'greed'; got '",
      backend,
      "'.",
      call. = FALSE
    )
  }
  if (!requireNamespace("greed", quietly = TRUE)) {
    stop(
      "backend = 'greed' but the 'greed' package is not installed. ",
      "Install with install.packages('greed').",
      call. = FALSE
    )
  }
  "greed"
}
