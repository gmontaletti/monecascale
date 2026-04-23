# S3 methods for auto_segment_levels
# ==================================
#
# Minimal console-surface methods. `format()` yields the one-line
# header; `print()` emits the header plus the per-level diagnostics
# table.

# 1. format ------------------------------------------------------------------

#' Format an `auto_segment_levels` Result
#'
#' One-line summary suitable for `cat()` contexts.
#'
#' @param x An object of class `"auto_segment_levels"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return A single character string.
#' @export
format.auto_segment_levels <- function(x, ...) {
  n_rows <- if (is.null(x$diagnostics)) 0L else nrow(x$diagnostics)
  sprintf(
    "<auto_segment_levels> method=%s backend=%s picked level=%d of %d",
    x$method %||% "NA",
    x$backend %||% "NA",
    as.integer(x$level),
    as.integer(n_rows)
  )
}

# 2. print -------------------------------------------------------------------

#' Print an `auto_segment_levels` Result
#'
#' Prints a header line followed by the per-level diagnostics table.
#'
#' @param x An object of class `"auto_segment_levels"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return `x`, invisibly.
#' @export
print.auto_segment_levels <- function(x, ...) {
  cat(format.auto_segment_levels(x), "\n", sep = "")
  if (!is.null(x$diagnostics) && nrow(x$diagnostics) > 0L) {
    print(x$diagnostics, row.names = FALSE)
  }
  invisible(x)
}
