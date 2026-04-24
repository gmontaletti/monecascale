# moneca_flow(): flow-based (Infomap / Map Equation) backend
# ==========================================================
#
# Alternative clustering backend for MONECA built on the Map Equation
# (Rosvall & Bergstrom 2008). Mobility is a natural flow process, so
# clustering under the map equation aligns directly with the substantive
# semantics: workers flow through jobs, modules trap flow.
#
# Only one backend is shipped: `"igraph"`, via `igraph::cluster_infomap()`.
# `cluster_infomap()` returns a *flat* partition; a recursive wrapper
# synthesises an approximate 2- or 3-level hierarchy by re-running the
# flat algorithm on each module's induced subgraph. This is not joint-MDL
# optimal (true hierarchical Infomap requires the C++ standalone binary)
# but is zero-cost in dependency terms and produces a usable hierarchy
# for downstream moneca consumers.

#' Flow-based (Infomap / Map Equation) Backend for MONECA Segmentation
#'
#' Alternative to [moneca::moneca_fast()] and [moneca_sbm()] that replaces
#' the clique / SBM step with Infomap (Rosvall & Bergstrom 2008). The Map
#' Equation treats the graph as a stationary random walk and extracts
#' modules that trap flow, which is a natural fit for mobility data
#' where edges are directed transition counts.
#'
#' @details
#' `igraph::cluster_infomap()` returns a flat partition. To expose a
#' hierarchy compatible with moneca's multi-level segment contract, a
#' recursive wrapper re-runs the flat algorithm on each module's induced
#' subgraph for level 2, and optionally once more for level 3 when
#' `depth = 3L`. This is an *approximate* hierarchy: it does not jointly
#' optimise a single hierarchical map equation. For joint-MDL
#' optimisation use the standalone C++ Infomap binary (planned as an
#' optional backend in a future release).
#'
#' Input is raw mobility counts. No relative-risk transformation is
#' applied; users who want RR-weighted flow clustering should preprocess
#' the matrix themselves, or use [moneca_sbm()] which consumes RR
#' natively via the DC-SBM bridge.
#'
#' The codelength per level is computed directly from the map equation
#' on the full graph (PageRank-based stationary distribution plus module
#' exit probabilities), not inferred from the sub-module Infomap runs.
#' This ensures `codelength_per_level` is comparable across levels of
#' the recursive hierarchy.
#'
#' @param mx A square mobility matrix, optionally with margins as the
#'   last row and column (auto-detected via `has_margins = "auto"`).
#'   Sparse `dgCMatrix` input is accepted and converted to dense for the
#'   igraph call.
#' @param backend Currently only `"igraph"`. Reserved for a future
#'   `"infomap_cpp"` backend that invokes the standalone Infomap C++
#'   binary for true hierarchical optimisation.
#' @param depth Integer in `1:3`. Number of recursion levels to build.
#'   `1` is flat; `2` (default) adds one level of within-module
#'   re-clustering; `3` adds a further level.
#' @param nb_trials Integer. Number of Infomap trials per call, passed
#'   through to `igraph::cluster_infomap()`. Higher values reduce
#'   stochastic variation at linear cost.
#' @param directed Logical. If `TRUE` (default), treat the graph as
#'   directed (mobility has intrinsic direction). If `FALSE`, the graph
#'   is undirected; combine with `symmetric_method = "sum"` or `"min"`
#'   to symmetrise the weights explicitly.
#' @param symmetric_method One of `"none"` (default; keep directed),
#'   `"sum"` (`mx + t(mx)`), or `"min"` (`2 * pmin(mx, t(mx))`). Applied
#'   to the core before the igraph fit.
#' @param small.cell.reduction Numeric. Cells in the core matrix below
#'   this value are zeroed before the fit. Mirrors
#'   [moneca::moneca_fast()].
#' @param zero_diagonal Logical. If `TRUE` (default), the diagonal of
#'   the core is zeroed before the fit (self-loops removed). Set to
#'   `FALSE` to preserve self-flows.
#' @param segment.levels Integer, `NULL`, or the string `"auto"`. If
#'   `NULL` (default), all levels produced by the recursion are kept.
#'   If an integer, the hierarchy is truncated to that many levels.
#'   If `"auto"`, the full hierarchy is fit and then
#'   [auto_segment_levels()] picks a single preferred level via
#'   `auto_method`; the object is trimmed to that level and carries the
#'   picker result under `$auto_level`.
#' @param auto_method One of `"mdl"` (default) or `"mi_plateau"`.
#'   Ignored unless `segment.levels = "auto"`. With `"mdl"` the criterion
#'   is applied to the codelength trace under the hood.
#' @param seed Integer seed for reproducible inference.
#' @param has_margins One of `"auto"` (default), `TRUE`, or `FALSE`.
#'   Controls margin detection and appending.
#' @param verbose Logical. If `TRUE`, prints per-level progress messages.
#' @param ... Reserved for forward compatibility.
#'
#' @return An object of class `"moneca"` with the same structure as
#'   [moneca::moneca_fast()] output, plus an extra slot
#'   `$flow_diagnostics` describing the backend used, the number of
#'   blocks per level, and the codelength trace. The diagnostics slot
#'   also carries `mdl_per_level` as an alias of
#'   `codelength_per_level` so downstream tooling that keys on
#'   `$mdl_per_level` can consume flow output unchanged.
#'
#' @seealso [moneca_sbm()], [moneca_bipartite()], [auto_segment_levels()],
#'   [moneca::moneca_fast()].
#'
#' @references
#' Rosvall, M. and Bergstrom, C. T. (2008). Maps of random walks on
#' complex networks reveal community structure. *PNAS*, 105(4),
#' 1118-1123.
#'
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   # Small planted-community fixture
#'   set.seed(42)
#'   n <- 20
#'   block <- rep(1:4, each = 5)
#'   mx <- matrix(rpois(n * n, 0.2), n, n)
#'   for (k in 1:4) {
#'     idx <- which(block == k)
#'     mx[idx, idx] <- mx[idx, idx] + rpois(length(idx)^2, 3)
#'   }
#'   diag(mx) <- 0
#'   fit <- moneca_flow(mx, depth = 2L, seed = 1L)
#'   fit$flow_diagnostics$codelength_per_level
#' }
#'
#' @export
moneca_flow <- function(
  mx,
  backend = "igraph",
  depth = 2L,
  nb_trials = 10L,
  directed = TRUE,
  symmetric_method = c("none", "sum", "min"),
  small.cell.reduction = 0,
  zero_diagonal = TRUE,
  segment.levels = NULL,
  auto_method = c("mdl", "mi_plateau"),
  seed = NULL,
  has_margins = "auto",
  verbose = FALSE,
  ...
) {
  symmetric_method <- match.arg(symmetric_method)
  auto_method <- match.arg(auto_method)

  # 1. Argument validation -------------------------------------------------
  depth <- as.integer(depth)
  if (length(depth) != 1L || is.na(depth) || depth < 1L || depth > 3L) {
    stop(
      "`depth` must be an integer in 1:3; got ",
      paste(depth, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  nb_trials <- as.integer(nb_trials)
  if (length(nb_trials) != 1L || is.na(nb_trials) || nb_trials < 1L) {
    stop("`nb_trials` must be a positive integer.", call. = FALSE)
  }
  stopifnot(
    is.logical(directed),
    length(directed) == 1L,
    is.logical(zero_diagonal),
    length(zero_diagonal) == 1L,
    is.logical(verbose),
    length(verbose) == 1L,
    is.numeric(small.cell.reduction),
    length(small.cell.reduction) == 1L,
    small.cell.reduction >= 0
  )

  # 2. Auto-level gate -----------------------------------------------------
  auto_mode <- identical(segment.levels, "auto")
  levels_for_fit <- if (auto_mode) NULL else segment.levels

  # 3. Margin handling -----------------------------------------------------
  mx_info <- .ensure_margins_sbm(mx, has_margins = has_margins)
  mx_full <- mx_info$mx
  margins_added <- mx_info$margins_added
  n_core <- nrow(mx_full) - 1L

  if (n_core < 3L) {
    stop("Need at least 3 non-margin rows/cols.", call. = FALSE)
  }

  # 4. Backend resolution --------------------------------------------------
  backend <- .resolve_flow_backend(backend)

  # 5. Preprocess core -----------------------------------------------------
  core_processed <- .preprocess_for_sbm(
    mx_full,
    small.cell.reduction = small.cell.reduction,
    symmetric_method = symmetric_method,
    zero_diagonal = zero_diagonal
  )
  # igraph's adjacency constructor doesn't digest dgCMatrix reliably
  # across versions; coerce to dense for the fit.
  if (inherits(core_processed, "sparseMatrix")) {
    core_processed <- as.matrix(core_processed)
  }
  rownames(core_processed) <- rownames(mx_full)[1:n_core]
  colnames(core_processed) <- colnames(mx_full)[1:n_core]

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # 6. Backend fit ---------------------------------------------------------
  fit <- switch(
    backend,
    "igraph" = .flow_fit_igraph(
      core_mx = core_processed,
      depth = depth,
      nb_trials = nb_trials,
      directed = directed,
      zero_diagonal = zero_diagonal,
      seed = seed,
      verbose = verbose
    ),
    stop(
      "Unsupported flow backend: '",
      backend,
      "'. Only 'igraph' is available in this version.",
      call. = FALSE
    )
  )

  # 7. Optional level truncation ------------------------------------------
  if (is.numeric(levels_for_fit)) {
    keep <- .select_sbm_levels(fit$n_blocks_per_level, levels_for_fit)
    fit$memberships <- fit$memberships[keep]
    fit$codelength_per_level <- fit$codelength_per_level[keep]
    fit$n_blocks_per_level <- fit$n_blocks_per_level[keep]
  }

  # 8. Adapter to moneca S3 object -----------------------------------------
  out <- .flow_fit_to_moneca(
    fit = fit,
    mx_full = mx_full,
    margins_added = margins_added,
    small.cell.reduction = small.cell.reduction,
    symmetric_method = symmetric_method,
    directed = directed
  )

  # 9. Post-hoc auto-level selection --------------------------------------
  if (auto_mode) {
    picked <- auto_segment_levels(out, method = auto_method)
    out <- .trim_moneca_to_level(out, picked$level)
    out$auto_level <- picked
  }

  out
}

# 10. Backend resolution ------------------------------------------------------

.resolve_flow_backend <- function(backend) {
  backend <- match.arg(backend, choices = c("igraph"))
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "moneca_flow(backend = 'igraph') requires the 'igraph' package. ",
      "Install with install.packages('igraph').",
      call. = FALSE
    )
  }
  backend
}
