# Shared utilities for moneca_bipartite() backends
# ================================================
#
# Adapter chain that turns a rectangular mobility matrix into a pair of
# moneca-class objects (row view / column view) wrapped in a
# moneca_bipartite S3 object. All helpers in this file are internal.

# 1. Margin handling ----------------------------------------------------------

.ensure_margins_bipartite <- function(mx, has_margins = "auto") {
  sparse_input <- inherits(mx, "sparseMatrix")
  if (!is.matrix(mx) && !sparse_input) {
    mx <- as.matrix(mx)
  }

  if (nrow(mx) < 2L || ncol(mx) < 2L) {
    stop(
      "Input matrix is too small: need at least 2 rows and 2 cols.",
      call. = FALSE
    )
  }

  detect <- function(mx) {
    n <- nrow(mx)
    m <- ncol(mx)
    if (n < 3L || m < 3L) {
      return(FALSE)
    }
    use_matrix <- inherits(mx, "sparseMatrix") &&
      requireNamespace("Matrix", quietly = TRUE)
    .rs <- if (use_matrix) Matrix::rowSums else rowSums
    .cs <- if (use_matrix) Matrix::colSums else colSums
    core <- mx[1:(n - 1L), 1:(m - 1L)]
    tol <- max(sum(core) * 1e-6, 1)
    last_row <- as.numeric(mx[n, 1:(m - 1L)])
    last_col <- as.numeric(mx[1:(n - 1L), m])
    all(abs(last_row - .cs(core)) < tol) &&
      all(abs(last_col - .rs(core)) < tol) &&
      abs(as.numeric(mx[n, m]) - sum(core)) < tol
  }

  if (identical(has_margins, "auto")) {
    already_has <- detect(mx)
  } else if (identical(has_margins, TRUE)) {
    already_has <- TRUE
  } else if (identical(has_margins, FALSE)) {
    already_has <- FALSE
  } else {
    stop(
      "has_margins must be 'auto', TRUE, or FALSE.",
      call. = FALSE
    )
  }

  if (already_has) {
    n <- nrow(mx)
    m <- ncol(mx)
    rect_core <- mx[1:(n - 1L), 1:(m - 1L)]
    if (sparse_input) {
      row_margin <- as.numeric(mx[1:(n - 1L), m])
      col_margin <- as.numeric(mx[n, 1:(m - 1L)])
      grand_total <- as.numeric(mx[n, m])
    } else {
      row_margin <- as.numeric(mx[1:(n - 1L), m])
      col_margin <- as.numeric(mx[n, 1:(m - 1L)])
      grand_total <- as.numeric(mx[n, m])
    }
    row_names <- rownames(mx)
    if (!is.null(row_names)) {
      row_names <- row_names[1:(n - 1L)]
    }
    col_names <- colnames(mx)
    if (!is.null(col_names)) {
      col_names <- col_names[1:(m - 1L)]
    }
    margins_added <- FALSE
  } else {
    rect_core <- mx
    if (sparse_input) {
      row_margin <- as.numeric(Matrix::rowSums(rect_core))
      col_margin <- as.numeric(Matrix::colSums(rect_core))
    } else {
      row_margin <- rowSums(rect_core)
      col_margin <- colSums(rect_core)
    }
    grand_total <- sum(rect_core)
    row_names <- rownames(rect_core)
    col_names <- colnames(rect_core)
    margins_added <- TRUE
  }

  # Ensure rect_core is a dgCMatrix for sparse paths.
  if (!inherits(rect_core, "sparseMatrix")) {
    rect_core <- Matrix::Matrix(rect_core, sparse = TRUE)
  }
  if (!inherits(rect_core, "dgCMatrix")) {
    rect_core <- as(rect_core, "generalMatrix")
    rect_core <- as(rect_core, "CsparseMatrix")
  }

  if (is.null(row_names)) {
    row_names <- paste0("Row", seq_len(nrow(rect_core)))
  }
  if (is.null(col_names)) {
    col_names <- paste0("Col", seq_len(ncol(rect_core)))
  }
  dimnames(rect_core) <- list(row_names, col_names)

  list(
    rect_core = rect_core,
    row_margin = row_margin,
    col_margin = col_margin,
    grand_total = grand_total,
    margins_added = margins_added,
    row_names = row_names,
    col_names = col_names
  )
}

# 2. Preprocess rectangular core ----------------------------------------------

.preprocess_for_bipartite <- function(rect_core, small.cell.reduction = 0) {
  if (!inherits(rect_core, "sparseMatrix")) {
    rect_core <- Matrix::Matrix(rect_core, sparse = TRUE)
  }
  if (!inherits(rect_core, "dgCMatrix")) {
    rect_core <- as(rect_core, "generalMatrix")
    rect_core <- as(rect_core, "CsparseMatrix")
  }

  if (small.cell.reduction > 0) {
    below <- rect_core@x < small.cell.reduction
    if (any(below)) {
      rect_core@x[below] <- 0
      rect_core <- Matrix::drop0(rect_core)
    }
  }

  rect_core
}

# 3. Rectangular relative-risk ------------------------------------------------

.bipartite_rr <- function(rect_core, row_margin, col_margin, grand_total) {
  if (!inherits(rect_core, "sparseMatrix")) {
    rect_core <- Matrix::Matrix(rect_core, sparse = TRUE)
  }
  if (!is.finite(grand_total) || grand_total <= 0) {
    stop(
      "grand_total must be a positive finite number.",
      call. = FALSE
    )
  }

  rect_t <- as(rect_core, "TsparseMatrix")
  # Zero-safe division: if either margin is zero at a nonzero cell position,
  # the expected value is zero and RR is ill-defined; drop those entries.
  rm <- row_margin[rect_t@i + 1L]
  cm <- col_margin[rect_t@j + 1L]
  expected <- (rm * cm) / grand_total
  safe <- is.finite(expected) & expected > 0
  new_x <- numeric(length(rect_t@x))
  new_x[safe] <- rect_t@x[safe] / expected[safe]

  rr <- Matrix::sparseMatrix(
    i = rect_t@i[safe] + 1L,
    j = rect_t@j[safe] + 1L,
    x = new_x[safe],
    dims = dim(rect_core),
    dimnames = dimnames(rect_core)
  )
  if (!inherits(rr, "dgCMatrix")) {
    rr <- as(rr, "generalMatrix")
    rr <- as(rr, "CsparseMatrix")
  }
  rr
}

# 4. One-mode RR projection ---------------------------------------------------

.one_mode_project_rr <- function(
  rr,
  side = c("row", "col"),
  row_margin,
  col_margin,
  grand_total,
  normalize = "margin_inv"
) {
  side <- match.arg(side)
  if (!identical(normalize, "margin_inv")) {
    stop(
      "Only normalize = 'margin_inv' is currently supported.",
      call. = FALSE
    )
  }

  if (side == "row") {
    # D_c^{-1} with grand_total scaling; pseudo-inverse on zero margins.
    inv <- ifelse(col_margin > 0, grand_total / col_margin, 0)
    D_inv <- Matrix::Diagonal(x = inv)
    P <- rr %*% D_inv %*% Matrix::t(rr)
  } else {
    inv <- ifelse(row_margin > 0, grand_total / row_margin, 0)
    D_inv <- Matrix::Diagonal(x = inv)
    P <- Matrix::t(rr) %*% D_inv %*% rr
  }

  if (!inherits(P, "dgCMatrix")) {
    P <- as(P, "generalMatrix")
    P <- as(P, "CsparseMatrix")
  }
  Matrix::diag(P) <- 0
  P <- Matrix::drop0(P)
  P
}

# 5. Append moneca-style margins to a sparse square projection ---------------

.append_margins_sparse <- function(P, side_names) {
  if (!inherits(P, "sparseMatrix")) {
    P <- Matrix::Matrix(P, sparse = TRUE)
  }
  rs <- as.numeric(Matrix::rowSums(P))
  cs <- as.numeric(Matrix::colSums(P))
  gt <- sum(rs)

  rs_col <- Matrix::Matrix(rs, ncol = 1L, sparse = TRUE)
  cs_row <- Matrix::Matrix(c(cs, gt), nrow = 1L, sparse = TRUE)
  P2 <- Matrix::cbind2(P, rs_col)
  P2 <- Matrix::rbind2(P2, cs_row)
  if (!inherits(P2, "dgCMatrix")) {
    P2 <- as(P2, "generalMatrix")
    P2 <- as(P2, "CsparseMatrix")
  }
  full_names <- c(side_names, "Total")
  dimnames(P2) <- list(full_names, full_names)
  P2
}

# 6. Per-level block-interaction matrix on the rectangular core --------------

.bipartite_block_interaction <- function(rect_core, memb_row, memb_col) {
  if (!inherits(rect_core, "sparseMatrix")) {
    rect_core <- Matrix::Matrix(rect_core, sparse = TRUE)
  }
  fr <- factor(
    as.integer(memb_row),
    levels = sort(unique(as.integer(memb_row)))
  )
  fc <- factor(
    as.integer(memb_col),
    levels = sort(unique(as.integer(memb_col)))
  )
  Gr <- Matrix::fac2sparse(fr)
  Gc <- Matrix::fac2sparse(fc)
  out <- as.matrix(Gr %*% rect_core %*% Matrix::t(Gc))
  dimnames(out) <- list(
    paste0("R", levels(fr)),
    paste0("C", levels(fc))
  )
  out
}

# 7. Build one side's moneca-class object ------------------------------------

.build_side_moneca <- function(
  memberships_side,
  P_full,
  side_names,
  small.cell.reduction,
  margin_policy,
  side_name = c("rows", "cols"),
  diag_side = NULL
) {
  side_name <- match.arg(side_name)
  n_side <- length(side_names)
  segment_list <- vector("list", length(memberships_side) + 1L)
  segment_list[[1]] <- as.list(seq_len(n_side))
  for (lvl in seq_along(memberships_side)) {
    segment_list[[lvl + 1L]] <- .memberships_to_cliques(
      memberships_side[[lvl]]
    )
  }

  mat_list <- vector("list", length(memberships_side) + 1L)
  mat_list[[1]] <- P_full
  for (lvl in seq_along(memberships_side)) {
    memb <- as.integer(memberships_side[[lvl]])
    mat_list[[lvl + 1L]] <- .aggregate_matrix_by_membership(P_full, memb)
  }

  obj <- list(
    segment.list = segment_list,
    mat.list = mat_list,
    small.cell.reduction = small.cell.reduction,
    margins_added = identical(margin_policy, "row_col"),
    density_reduction = NULL,
    symmetric_method = "sum"
  )
  class(obj) <- "moneca"
  obj$bipartite_origin <- side_name
  if (!is.null(diag_side)) {
    obj$bipartite_diagnostics_side <- diag_side
  }
  obj$segment_metadata <- moneca::moneca_segments(obj)
  obj
}

# 8. Main adapter: fit + rect_info => moneca_bipartite S3 object --------------

.fit_to_moneca_bipartite <- function(
  fit,
  rect_info,
  rr,
  small.cell.reduction,
  margin_policy,
  verbose = FALSE
) {
  rect_core <- rect_info$rect_core
  row_margin <- rect_info$row_margin
  col_margin <- rect_info$col_margin
  grand_total <- rect_info$grand_total
  row_names <- rect_info$row_names
  col_names <- rect_info$col_names

  # 8a. One-mode RR projections per side ----------------------------------
  P_row <- .one_mode_project_rr(
    rr = rr,
    side = "row",
    row_margin = row_margin,
    col_margin = col_margin,
    grand_total = grand_total,
    normalize = "margin_inv"
  )
  P_col <- .one_mode_project_rr(
    rr = rr,
    side = "col",
    row_margin = row_margin,
    col_margin = col_margin,
    grand_total = grand_total,
    normalize = "margin_inv"
  )

  # 8b. Append moneca-style margins ---------------------------------------
  if (identical(margin_policy, "row_col")) {
    P_row_full <- .append_margins_sparse(P_row, row_names)
    P_col_full <- .append_margins_sparse(P_col, col_names)
  } else {
    # No margins appended; still ensure dimnames preserved on the square.
    if (!inherits(P_row, "dgCMatrix")) {
      P_row <- as(P_row, "CsparseMatrix")
    }
    if (!inherits(P_col, "dgCMatrix")) {
      P_col <- as(P_col, "CsparseMatrix")
    }
    dimnames(P_row) <- list(row_names, row_names)
    dimnames(P_col) <- list(col_names, col_names)
    P_row_full <- P_row
    P_col_full <- P_col
  }

  # 8c. Per-side slim diagnostics slice ------------------------------------
  joint_icl <- fit$joint_icl_per_level
  joint_mdl <- -joint_icl
  diag_row <- list(
    joint_icl_per_level = joint_icl,
    joint_mdl_per_level = joint_mdl,
    n_blocks_per_level = fit$n_blocks_row
  )
  diag_col <- list(
    joint_icl_per_level = joint_icl,
    joint_mdl_per_level = joint_mdl,
    n_blocks_per_level = fit$n_blocks_col
  )

  # 8d. Build each side's moneca-class object ------------------------------
  rows_obj <- .build_side_moneca(
    memberships_side = fit$memberships_row,
    P_full = P_row_full,
    side_names = row_names,
    small.cell.reduction = small.cell.reduction,
    margin_policy = margin_policy,
    side_name = "rows",
    diag_side = diag_row
  )
  cols_obj <- .build_side_moneca(
    memberships_side = fit$memberships_col,
    P_full = P_col_full,
    side_names = col_names,
    small.cell.reduction = small.cell.reduction,
    margin_policy = margin_policy,
    side_name = "cols",
    diag_side = diag_col
  )

  # 8e. Per-level block-interaction matrices -------------------------------
  n_levels <- length(fit$memberships_row)
  block_interaction <- vector("list", n_levels)
  for (lvl in seq_len(n_levels)) {
    block_interaction[[lvl]] <- .bipartite_block_interaction(
      rect_core,
      fit$memberships_row[[lvl]],
      fit$memberships_col[[lvl]]
    )
  }

  # 8f. Diagnostics slot ----------------------------------------------------
  bipartite_diagnostics <- list(
    backend = fit$backend %||% "greed",
    model = fit$model %||% "DcLbm",
    n_row = length(row_names),
    n_col = length(col_names),
    n_blocks_row_per_level = fit$n_blocks_row,
    n_blocks_col_per_level = fit$n_blocks_col,
    joint_icl_per_level = fit$joint_icl_per_level,
    joint_mdl_per_level = -fit$joint_icl_per_level,
    block_interaction_matrix_per_level = block_interaction,
    rr_rect = rr,
    row_margin = row_margin,
    col_margin = col_margin,
    grand_total = grand_total,
    rr_scale = fit$rr_scale %||% 1L,
    n_init = fit$n_init %||% 1L,
    margins_added = rect_info$margins_added,
    row_names = row_names,
    col_names = col_names
  )

  out <- list(
    rows = rows_obj,
    cols = cols_obj,
    bipartite_diagnostics = bipartite_diagnostics
  )
  class(out) <- "moneca_bipartite"
  out
}
