# S3 methods for moneca_bipartite
# ================================
#
# Minimal console-surface methods for the moneca_bipartite class. They do
# not replace moneca's own print / summary methods on the embedded $rows
# and $cols moneca-class objects; those still dispatch normally.

# 1. format -------------------------------------------------------------------

#' Format a `moneca_bipartite` object
#'
#' One-line summary of a bipartite fit suitable for `cat()` contexts.
#'
#' @param x An object of class `"moneca_bipartite"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return A single character string describing the bipartite fit.
#' @export
format.moneca_bipartite <- function(x, ...) {
  diag <- x$bipartite_diagnostics
  backend <- diag$backend %||% "greed"
  model <- diag$model %||% "DcLbm"
  sprintf(
    "<moneca_bipartite> rows: %d x cols: %d, backend: %s/%s, levels: %d",
    as.integer(diag$n_row %||% NA_integer_),
    as.integer(diag$n_col %||% NA_integer_),
    backend,
    model,
    length(diag$n_blocks_row_per_level %||% integer(0))
  )
}

# 2. print --------------------------------------------------------------------

#' Print a `moneca_bipartite` object
#'
#' Prints a header line followed by a per-level table of row block count,
#' column block count, and joint ICL.
#'
#' @param x An object of class `"moneca_bipartite"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.moneca_bipartite <- function(x, ...) {
  cat(format.moneca_bipartite(x), "\n", sep = "")
  diag <- x$bipartite_diagnostics
  krow <- diag$n_blocks_row_per_level
  kcol <- diag$n_blocks_col_per_level
  icl <- diag$joint_icl_per_level
  n_lvl <- length(krow)
  if (n_lvl == 0L) {
    cat("No non-trivial hierarchy levels recovered.\n")
    return(invisible(x))
  }
  tbl <- data.frame(
    level = seq_len(n_lvl),
    Krow = as.integer(krow),
    Kcol = as.integer(kcol),
    joint_ICL = icl,
    stringsAsFactors = FALSE
  )
  print(tbl, row.names = FALSE)
  invisible(x)
}

# 3. summary ------------------------------------------------------------------

#' Summary of a `moneca_bipartite` object
#'
#' Returns a structured summary of the bipartite fit, including the
#' per-level table, backend metadata, and side-level block counts.
#'
#' @param object An object of class `"moneca_bipartite"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return An object of class `"summary.moneca_bipartite"` with elements
#'   `header`, `levels`, and `meta`.
#' @export
summary.moneca_bipartite <- function(object, ...) {
  diag <- object$bipartite_diagnostics
  krow <- diag$n_blocks_row_per_level
  kcol <- diag$n_blocks_col_per_level
  icl <- diag$joint_icl_per_level
  n_lvl <- length(krow)

  levels_tbl <- if (n_lvl > 0L) {
    data.frame(
      level = seq_len(n_lvl),
      Krow = as.integer(krow),
      Kcol = as.integer(kcol),
      joint_ICL = icl,
      joint_MDL = -icl,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      level = integer(0),
      Krow = integer(0),
      Kcol = integer(0),
      joint_ICL = numeric(0),
      joint_MDL = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  out <- list(
    header = format.moneca_bipartite(object),
    levels = levels_tbl,
    meta = list(
      backend = diag$backend %||% "greed",
      model = diag$model %||% "DcLbm",
      n_row = diag$n_row,
      n_col = diag$n_col,
      n_init = diag$n_init %||% 1L,
      rr_scale = diag$rr_scale %||% 1L,
      margins_added = diag$margins_added %||% NA
    )
  )
  class(out) <- "summary.moneca_bipartite"
  out
}

# 4. print.summary.moneca_bipartite ------------------------------------------

#' Print a `summary.moneca_bipartite` object
#'
#' @param x An object of class `"summary.moneca_bipartite"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.summary.moneca_bipartite <- function(x, ...) {
  cat(x$header, "\n", sep = "")
  cat(
    sprintf(
      "  backend: %s   model: %s   n_init: %d   rr_scale: %d\n",
      x$meta$backend,
      x$meta$model,
      as.integer(x$meta$n_init),
      as.integer(x$meta$rr_scale)
    )
  )
  if (nrow(x$levels) == 0L) {
    cat("No non-trivial hierarchy levels recovered.\n")
  } else {
    cat("Hierarchy levels:\n")
    print(x$levels, row.names = FALSE)
  }
  invisible(x)
}
