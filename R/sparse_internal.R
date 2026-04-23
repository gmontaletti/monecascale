# Local copies of two internal sparse helpers from moneca.
# ========================================================
#
# Both helpers are `@keywords internal @noRd` in moneca (not exported),
# so we ship local copies rather than reach into the moneca namespace
# via `:::` (which is frowned upon by CRAN).

# 1. Null-coalescing helper ---------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# 2. Symmetric sparse pmin ----------------------------------------------------
# Pairwise minimum of a sparse matrix and its transpose, retaining only
# structurally bilateral positions.

.sparse_pmin_symmetric <- function(m) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Sparse symmetrisation requires the 'Matrix' package.", call. = FALSE)
  }
  if (!inherits(m, "sparseMatrix")) {
    stop("Input must be a sparseMatrix.", call. = FALSE)
  }

  ma <- as(m, "TsparseMatrix")
  mb <- as(Matrix::t(m), "TsparseMatrix")

  dims <- dim(m)
  dn <- dimnames(m)

  if (length(ma@x) == 0L || length(mb@x) == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = dims,
      dimnames = dn
    ))
  }

  df_a <- data.frame(i = ma@i, j = ma@j, a = ma@x)
  df_b <- data.frame(i = mb@i, j = mb@j, b = mb@x)
  merged <- merge(df_a, df_b, by = c("i", "j"))

  if (nrow(merged) == 0L) {
    return(Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = dims,
      dimnames = dn
    ))
  }

  Matrix::sparseMatrix(
    i = merged$i + 1L,
    j = merged$j + 1L,
    x = pmin(merged$a, merged$b),
    dims = dims,
    dimnames = dn
  )
}

# 3. Sparse segment-matrix aggregation ----------------------------------------
# Aggregates a sparse mobility matrix by membership vector, mirroring the
# dense rowsum path in moneca_fast's segment.matrix.fast (margin row labelled
# length(membership) + 1).

.segment_matrix_sparse <- function(mx, segments) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Sparse aggregation requires the 'Matrix' package.", call. = FALSE)
  }
  if (!inherits(mx, "sparseMatrix")) {
    mx <- Matrix::Matrix(mx, sparse = TRUE)
  }

  memb_int <- as.integer(segments$membership)
  groups.1 <- c(memb_int, length(memb_int) + 1L)
  f <- factor(groups.1, levels = sort(unique(groups.1)))

  G <- Matrix::fac2sparse(f)
  res <- G %*% mx %*% Matrix::t(G)

  grp_names <- levels(f)
  dimnames(res) <- list(grp_names, grp_names)

  if (!inherits(res, "dgCMatrix")) {
    res <- as(res, "generalMatrix")
    res <- as(res, "CsparseMatrix")
  }
  res
}
