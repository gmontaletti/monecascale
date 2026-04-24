# auto_segment_levels(): auto-select a preferred hierarchy level
# ==============================================================
#
# Two criteria in v0.3.0: MDL elbow and level-to-level MI plateau.
# Operates on moneca_sbm() output and on the $rows / $cols children of
# moneca_bipartite() output. Plain moneca::moneca_fast() output is
# rejected with an informative error and deferred to 0.3.1.

#' Auto-select a Hierarchy Level for a `moneca` Fit
#'
#' @description
#' Picks a preferred level in a hierarchical moneca clustering based on
#' one of two criteria: elbow on the Minimum Description Length (MDL)
#' trace, or plateau detection on level-to-level mutual information.
#' Dispatches on [moneca_sbm()] output directly, and recurses into the
#' `$rows` and `$cols` children of a [moneca_bipartite()] object.
#'
#' @details
#' The `"mdl"` criterion fetches the per-level MDL series from
#' `obj$sbm_diagnostics$mdl_per_level` (for SBM backends) or from the
#' slim per-side slice `obj$bipartite_diagnostics_side$joint_mdl_per_level`
#' (for a `$rows` or `$cols` child of a `moneca_bipartite` object). The
#' series is reordered to coarse-to-fine (ascending block count) and an
#' elbow is selected by either the Kneedle heuristic or the discrete
#' second-difference peak.
#'
#' The `"mi_plateau"` criterion computes level-to-level sample mutual
#' information from `obj$segment.list` directly (no reliance on backend
#' diagnostics). It picks the smallest level at which the gain
#' \eqn{I_l - I_{l-1}} falls below `plateau_tol * max(I)`; if no plateau is
#' detected, it returns the finest level.
#'
#' For `moneca::moneca_fast()` output the function errors out: neither
#' an MDL series nor a principled MI plateau is available without a
#' per-iteration trace, which upstream `moneca_fast()` does not emit.
#' Support is planned for monecascale 0.3.1.
#'
#' @param obj A `moneca` object produced by [moneca_sbm()], a
#'   `moneca_bipartite` object produced by [moneca_bipartite()], or a
#'   `$rows` / `$cols` child of such an object.
#' @param method One of `"mdl"` (default) or `"mi_plateau"`.
#' @param return One of `"full"` (default, returns an
#'   `auto_segment_levels` object) or `"level"` (returns a bare integer
#'   scalar).
#' @param plateau_tol Positive numeric in `(0, 1)`. Relative gain
#'   threshold under which the MI curve is considered flat. Only used
#'   when `method = "mi_plateau"`.
#' @param mdl_elbow One of `"kneedle"` (default) or `"max_second_diff"`.
#'   Only used when `method = "mdl"`.
#' @param verbose Logical. Emit per-level diagnostic messages.
#' @param ... Reserved for forward compatibility.
#'
#' @return When `obj` inherits from `"moneca_bipartite"`, a list with
#'   named components `rows` and `cols`, each an `auto_segment_levels`
#'   object (or integer scalar if `return = "level"`). Otherwise, when
#'   `return = "full"`, an object of class `auto_segment_levels` with
#'   components:
#'   * `level` - integer scalar, the selected `segment.list` index.
#'   * `method` - the criterion used.
#'   * `diagnostics` - a data.frame with columns `level`, `n_blocks`,
#'     `mdl`, `mi_to_next`, `score`.
#'   * `backend` - one of `"sbm"`, `"bipartite_rows"`, `"bipartite_cols"`.
#'   * `call` - the matched call.
#'
#'   When `return = "level"`, a single integer scalar.
#'
#' @seealso [moneca_sbm()], [moneca_bipartite()].
#'
#' @examples
#' if (requireNamespace("greed", quietly = TRUE)) {
#'   mob <- moneca::generate_mobility_data(n_classes = 20, seed = 1L)
#'   fit <- moneca_sbm(mob, backend = "greed", seed = 1L)
#'   auto_segment_levels(fit, method = "mdl")
#' }
#'
#' @export
auto_segment_levels <- function(
  obj,
  method = c("mdl", "mi_plateau"),
  return = c("full", "level"),
  plateau_tol = 0.05,
  mdl_elbow = c("kneedle", "max_second_diff"),
  verbose = FALSE,
  ...
) {
  method <- match.arg(method)
  return_mode <- match.arg(return)
  mdl_elbow <- match.arg(mdl_elbow)
  stopifnot(
    is.numeric(plateau_tol),
    length(plateau_tol) == 1L,
    plateau_tol > 0,
    plateau_tol < 1,
    is.logical(verbose),
    length(verbose) == 1L
  )

  # 1. Bipartite top-level dispatch -----------------------------------------
  if (inherits(obj, "moneca_bipartite")) {
    return(list(
      rows = auto_segment_levels(
        obj$rows,
        method = method,
        return = return_mode,
        plateau_tol = plateau_tol,
        mdl_elbow = mdl_elbow,
        verbose = verbose,
        ...
      ),
      cols = auto_segment_levels(
        obj$cols,
        method = method,
        return = return_mode,
        plateau_tol = plateau_tol,
        mdl_elbow = mdl_elbow,
        verbose = verbose,
        ...
      )
    ))
  }

  # 2. moneca guard ---------------------------------------------------------
  if (!inherits(obj, "moneca")) {
    stop(
      "auto_segment_levels() requires a `moneca` or `moneca_bipartite` ",
      "object; got class '",
      paste(class(obj), collapse = "/"),
      "'.",
      call. = FALSE
    )
  }

  # 3. Backend detection ---------------------------------------------------
  backend <- .detect_backend(obj)

  # 4. Criterion dispatch ---------------------------------------------------
  result <- switch(
    method,
    mdl = .criterion_mdl(
      obj,
      backend = backend,
      elbow = mdl_elbow,
      verbose = verbose
    ),
    mi_plateau = .criterion_mi_plateau(
      obj,
      plateau_tol = plateau_tol,
      verbose = verbose
    )
  )

  # 5. Decorate -------------------------------------------------------------
  result$backend <- backend
  result$call <- match.call()
  class(result) <- "auto_segment_levels"

  if (return_mode == "level") {
    return(as.integer(result$level))
  }
  result
}

# 6. Backend detection -------------------------------------------------------

#' @keywords internal
#' @noRd
.detect_backend <- function(obj) {
  if (!is.null(obj$flow_diagnostics)) {
    return("flow")
  }
  if (!is.null(obj$sbm_diagnostics)) {
    return("sbm")
  }
  if (!is.null(obj$bipartite_origin)) {
    origin <- obj$bipartite_origin
    if (identical(origin, "rows")) {
      return("bipartite_rows")
    }
    if (identical(origin, "cols")) {
      return("bipartite_cols")
    }
  }
  "fast"
}

# 7. Trim a moneca object to a target level ---------------------------------

#' @keywords internal
#' @noRd
.trim_moneca_to_level <- function(obj, level) {
  L <- length(obj$segment.list)
  stopifnot(
    is.numeric(level),
    length(level) == 1L,
    level >= 1L,
    level <= L
  )
  level <- as.integer(level)
  if (level == L) {
    return(obj)
  }
  obj$segment.list <- obj$segment.list[seq_len(level)]
  obj$mat.list <- obj$mat.list[seq_len(level)]
  # $sbm_diagnostics / $bipartite_diagnostics_side are kept whole as a
  # reference trace; they describe the full hierarchy that was fit.
  obj$segment_metadata <- moneca::moneca_segments(obj)
  obj
}
