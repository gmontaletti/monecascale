# greed backend for moneca_bipartite()
# ====================================
#
# Wraps greed::greed() with DcLbm() on the rectangular relative-risk matrix
# and exposes its hierarchical co-clustering as per-level row/column
# membership vectors consumable by the bipartite adapter.
#
# greed fits a flat DcLbm at the ICL-optimal joint K (Krow + Kcol), then
# builds a full nested partition tree. We expose levels from sol@K down to
# 2 (the K = 1 level is trivial).

# 1. Main fit function --------------------------------------------------------

.bipartite_fit_greed <- function(
  rr,
  rect_core_sum,
  seed = NULL,
  n_init = 1L,
  max_K = NULL,
  segment.levels = NULL,
  verbose = FALSE,
  ...
) {
  if (!requireNamespace("greed", quietly = TRUE)) {
    stop(
      "moneca_bipartite(backend = 'greed') requires the 'greed' package. ",
      "Install it with install.packages('greed').",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # 2. Scale rectangular RR to integer counts for Poisson DcLbm --------------
  rr_sum <- if (inherits(rr, "sparseMatrix")) {
    sum(rr@x)
  } else {
    sum(rr)
  }
  scale <- if (is.finite(rr_sum) && rr_sum > 0) {
    as.integer(round(rect_core_sum / rr_sum))
  } else {
    1L
  }
  if (!is.finite(scale) || scale < 1L) {
    scale <- 1L
  }

  rr_scaled <- rr
  if (inherits(rr_scaled, "sparseMatrix")) {
    rr_scaled@x <- as.numeric(round(scale * rr_scaled@x))
    rr_scaled <- Matrix::drop0(rr_scaled)
  } else {
    rr_scaled <- round(scale * rr_scaled)
  }

  # 3. Best-of-n_init fit ---------------------------------------------------
  n_init <- as.integer(n_init)
  if (is.na(n_init) || n_init < 1L) {
    n_init <- 1L
  }

  fits <- vector("list", n_init)
  icls <- numeric(n_init)
  for (i in seq_len(n_init)) {
    sol_i <- suppressWarnings(suppressMessages(greed::greed(
      rr_scaled,
      model = greed::DcLbm(),
      alg = greed::Hybrid(),
      verbose = verbose
    )))
    fits[[i]] <- sol_i
    icls[i] <- tryCatch(
      suppressWarnings(greed::ICL(sol_i)),
      error = function(e) NA_real_
    )
  }
  best_idx <- which.max(icls)
  if (length(best_idx) == 0L) {
    best_idx <- 1L
  }
  sol <- fits[[best_idx]]

  K_final <- sol@K
  if (is.null(K_final) || length(K_final) == 0L || K_final < 2L) {
    return(list(
      memberships_row = list(),
      memberships_col = list(),
      n_blocks_row = integer(0),
      n_blocks_col = integer(0),
      joint_icl_per_level = numeric(0),
      backend = "greed",
      model = "DcLbm",
      n_init = n_init,
      rr_scale = scale
    ))
  }

  # 4. Enumerate per-level partitions ---------------------------------------
  k_seq <- seq.int(from = as.integer(K_final), to = 2L, by = -1L)

  memb_row <- vector("list", length(k_seq))
  memb_col <- vector("list", length(k_seq))
  krow_vec <- integer(length(k_seq))
  kcol_vec <- integer(length(k_seq))
  icl_vec <- numeric(length(k_seq))

  for (j in seq_along(k_seq)) {
    k <- k_seq[j]
    cut_sol <- tryCatch(
      suppressWarnings(greed::cut(sol, k)),
      error = function(e) NULL
    )
    if (is.null(cut_sol)) {
      memb_row[[j]] <- rep(1L, nrow(rr_scaled))
      memb_col[[j]] <- rep(1L, ncol(rr_scaled))
      krow_vec[j] <- 1L
      kcol_vec[j] <- 1L
      icl_vec[j] <- NA_real_
      next
    }
    memb_row[[j]] <- as.integer(cut_sol@clrow)
    memb_col[[j]] <- as.integer(cut_sol@clcol)
    krow_vec[j] <- as.integer(cut_sol@Krow)
    kcol_vec[j] <- as.integer(cut_sol@Kcol)
    icl_vec[j] <- tryCatch(
      suppressWarnings(greed::ICL(cut_sol)),
      error = function(e) NA_real_
    )
  }

  # 5. Drop degenerate levels ----------------------------------------------
  # A level is kept if at least one side has a non-singleton block and at
  # least one side has more than one block (otherwise both views collapse).
  keep <- vapply(
    seq_along(k_seq),
    function(j) {
      has_blk <- krow_vec[j] > 1L || kcol_vec[j] > 1L
      has_clique <- any(table(memb_row[[j]]) >= 2L) ||
        any(table(memb_col[[j]]) >= 2L)
      has_blk && has_clique
    },
    logical(1)
  )

  memb_row <- memb_row[keep]
  memb_col <- memb_col[keep]
  krow_vec <- krow_vec[keep]
  kcol_vec <- kcol_vec[keep]
  icl_vec <- icl_vec[keep]

  # 6. Optional truncation via user-requested segment.levels ---------------
  if (length(memb_row) > 0L && !is.null(segment.levels)) {
    joint_k <- as.integer(krow_vec) + as.integer(kcol_vec)
    picked <- .select_sbm_levels(joint_k, segment.levels)
    memb_row <- memb_row[picked]
    memb_col <- memb_col[picked]
    krow_vec <- krow_vec[picked]
    kcol_vec <- kcol_vec[picked]
    icl_vec <- icl_vec[picked]
  }

  list(
    memberships_row = memb_row,
    memberships_col = memb_col,
    n_blocks_row = as.integer(krow_vec),
    n_blocks_col = as.integer(kcol_vec),
    joint_icl_per_level = icl_vec,
    backend = "greed",
    model = "DcLbm",
    n_init = n_init,
    rr_scale = scale
  )
}
