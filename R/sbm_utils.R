# Shared utilities for moneca_sbm() backends
# ==========================================
#
# Adapter that turns backend-agnostic SBM output into a moneca-class S3
# object, plus margin / symmetry / diagonal preprocessing shared across
# backends. All helpers in this file are internal.

# 1. Margin handling ----------------------------------------------------------

.ensure_margins_sbm <- function(mx, has_margins = "auto") {
  sparse_input <- inherits(mx, "sparseMatrix")
  if (!is.matrix(mx) && !sparse_input) {
    mx <- as.matrix(mx)
  }
  if (nrow(mx) != ncol(mx)) {
    stop(
      "moneca_sbm() currently requires a square matrix (rectangular ",
      "person-employer input is deferred to a future moneca_bipartite(); ",
      "see the MONECA scaling roadmap direction D1).",
      call. = FALSE
    )
  }

  detect <- function(mx) {
    n <- nrow(mx)
    if (n < 3) {
      return(FALSE)
    }
    use_matrix <- inherits(mx, "sparseMatrix") &&
      requireNamespace("Matrix", quietly = TRUE)
    .rs <- if (use_matrix) Matrix::rowSums else rowSums
    .cs <- if (use_matrix) Matrix::colSums else colSums
    core <- mx[1:(n - 1), 1:(n - 1)]
    tol <- max(sum(core) * 1e-6, 1)
    last_row <- as.numeric(mx[n, 1:(n - 1)])
    last_col <- as.numeric(mx[1:(n - 1), n])
    all(abs(last_row - .cs(core)) < tol) &&
      all(abs(last_col - .rs(core)) < tol) &&
      abs(as.numeric(mx[n, n]) - sum(core)) < tol
  }

  needs_margins <- if (identical(has_margins, FALSE)) {
    TRUE
  } else if (identical(has_margins, "auto")) {
    !detect(mx)
  } else {
    FALSE
  }

  if (!needs_margins) {
    return(list(mx = mx, margins_added = FALSE))
  }

  core <- mx
  if (sparse_input) {
    rs <- as.numeric(Matrix::rowSums(core))
    cs <- as.numeric(Matrix::colSums(core))
  } else {
    rs <- rowSums(core)
    cs <- colSums(core)
  }
  gt <- sum(core)

  cat_names <- rownames(core)
  if (is.null(cat_names)) {
    cat_names <- paste0("Cat", seq_len(nrow(core)))
  }
  full_names <- c(cat_names, "Total")

  if (sparse_input) {
    rs_col <- Matrix::Matrix(rs, ncol = 1, sparse = TRUE)
    cs_row <- Matrix::Matrix(c(cs, gt), nrow = 1, sparse = TRUE)
    mx2 <- Matrix::cbind2(core, rs_col)
    mx2 <- Matrix::rbind2(mx2, cs_row)
    if (!inherits(mx2, "dgCMatrix")) {
      mx2 <- as(mx2, "generalMatrix")
      mx2 <- as(mx2, "CsparseMatrix")
    }
    dimnames(mx2) <- list(full_names, full_names)
  } else {
    mx2 <- rbind(cbind(core, rs), c(cs, gt))
    rownames(mx2) <- full_names
    colnames(mx2) <- full_names
  }

  list(mx = mx2, margins_added = TRUE)
}

# 2. Core preprocessing for SBM backends --------------------------------------

.preprocess_for_sbm <- function(
  mx_full,
  small.cell.reduction = 0,
  symmetric_method = "sum",
  zero_diagonal = TRUE
) {
  n <- nrow(mx_full)
  core <- mx_full[1:(n - 1), 1:(n - 1)]

  if (small.cell.reduction > 0) {
    core[core < small.cell.reduction] <- 0
  }

  if (zero_diagonal) {
    if (inherits(core, "sparseMatrix")) {
      Matrix::diag(core) <- 0
      core <- Matrix::drop0(core)
    } else {
      diag(core) <- 0
    }
  }

  core <- switch(
    symmetric_method,
    "sum" = core + t(core),
    "min" = {
      if (inherits(core, "sparseMatrix")) {
        sym <- .sparse_pmin_symmetric(core)
        sym * 2
      } else {
        pmin(core, t(core)) * 2
      }
    },
    "none" = core,
    stop(
      "symmetric_method must be one of 'sum', 'min', 'none'; got '",
      symmetric_method,
      "'",
      call. = FALSE
    )
  )

  core
}

# 3. Convert a per-level membership vector to a moneca clique list ------------

.memberships_to_cliques <- function(memb) {
  memb <- as.integer(memb)
  blocks <- split(seq_along(memb), memb)
  blocks <- lapply(blocks, as.integer)
  blocks <- blocks[lengths(blocks) > 1L]
  unname(blocks)
}

# 4. Aggregate a matrix by a membership vector -------------------------------

.aggregate_matrix_by_membership <- function(mx, membership) {
  if (inherits(mx, "sparseMatrix")) {
    return(.segment_matrix_sparse(mx, list(membership = membership)))
  }
  membership <- as.integer(membership)
  groups.1 <- c(membership, length(membership) + 1L)
  mx.2_r <- rowsum(mx, groups.1)
  mx.2_rc <- rowsum(t(mx.2_r), groups.1)
  t(mx.2_rc)
}

# 5. Select a subset of SBM levels when user caps segment.levels -------------

.select_sbm_levels <- function(available_K, segment.levels = NULL) {
  if (is.null(segment.levels)) {
    return(seq_along(available_K))
  }
  n_avail <- length(available_K)
  if (n_avail <= segment.levels) {
    return(seq_len(n_avail))
  }
  idx <- unique(round(seq(1, n_avail, length.out = segment.levels)))
  sort(idx)
}

# 6. Build isolates_summary matching moneca_fast() ----------------------------

.build_isolates_summary <- function(mx_full, segment_list) {
  n <- nrow(mx_full) - 1L
  category_names <- rownames(mx_full)[1:n]
  if (is.null(category_names)) {
    category_names <- as.character(seq_len(n))
  }

  final_level <- length(segment_list)
  final_segments <- segment_list[[final_level]]

  membership <- rep("altri", n)
  for (seg_idx in seq_along(final_segments)) {
    membership[final_segments[[seg_idx]]] <- paste0("Segment_", seg_idx)
  }

  membership_df <- data.frame(
    name = category_names,
    group = membership,
    stringsAsFactors = FALSE
  )

  groups <- unique(membership)
  group_factor <- factor(membership, levels = groups)
  core_mx <- mx_full[1:n, 1:n]

  if (
    inherits(core_mx, "sparseMatrix") &&
      requireNamespace("Matrix", quietly = TRUE)
  ) {
    G <- Matrix::fac2sparse(group_factor)
    mobility_matrix <- as.matrix(G %*% core_mx %*% Matrix::t(G))
    dimnames(mobility_matrix) <- list(groups, groups)
  } else {
    if (inherits(core_mx, "sparseMatrix")) {
      core_mx <- as.matrix(core_mx)
    }
    mobility_matrix <- matrix(
      0,
      nrow = length(groups),
      ncol = length(groups),
      dimnames = list(groups, groups)
    )
    for (i in seq_along(groups)) {
      for (j in seq_along(groups)) {
        rows_i <- which(group_factor == groups[i])
        cols_j <- which(group_factor == groups[j])
        mobility_matrix[i, j] <- sum(core_mx[rows_i, cols_j])
      }
    }
  }

  list(membership = membership_df, mobility_matrix = mobility_matrix)
}

# 7. Main adapter: sbm_fit + mx_full => moneca S3 object ----------------------

.sbm_fit_to_moneca <- function(
  sbm_fit,
  mx_full,
  small.cell.reduction,
  margins_added,
  symmetric_method,
  segment.levels = NULL,
  isolates = FALSE
) {
  n_core <- nrow(mx_full) - 1L

  picked <- .select_sbm_levels(sbm_fit$n_blocks, segment.levels)
  memberships_picked <- sbm_fit$memberships[picked]
  n_blocks_picked <- sbm_fit$n_blocks[picked]
  mdl_picked <- sbm_fit$mdl[picked]

  segment_list <- vector("list", length(memberships_picked) + 1L)
  segment_list[[1]] <- as.list(seq_len(n_core))
  for (lvl in seq_along(memberships_picked)) {
    segment_list[[lvl + 1L]] <- .memberships_to_cliques(
      memberships_picked[[lvl]]
    )
  }

  mat_list <- vector("list", length(memberships_picked) + 1L)
  mat_list[[1]] <- mx_full
  for (lvl in seq_along(memberships_picked)) {
    memb <- as.integer(memberships_picked[[lvl]])
    mat_list[[lvl + 1L]] <- .aggregate_matrix_by_membership(
      mat_list[[1]],
      memb
    )
  }

  out <- list(
    segment.list = segment_list,
    mat.list = mat_list,
    small.cell.reduction = small.cell.reduction,
    margins_added = margins_added,
    density_reduction = NULL,
    symmetric_method = symmetric_method,
    sbm_diagnostics = list(
      backend = sbm_fit$backend,
      n_blocks_per_level = n_blocks_picked,
      mdl_per_level = mdl_picked,
      n_init = sbm_fit$n_init %||% 1L,
      node_names = sbm_fit$node_names
    )
  )

  if (isolates) {
    out$isolates_summary <- .build_isolates_summary(mx_full, segment_list)
  }

  class(out) <- "moneca"
  # Reuse moneca's canonical metadata builder for downstream plotting
  out$segment_metadata <- moneca::moneca_segments(out)
  out
}
