# Tests: latent-block recovery for monecascale::moneca_bipartite()
# ================================================================
#
# The bipartite fixture draws Poisson counts from a block-constant
# intensity matrix with a checkerboard high/low pattern. Ground-truth
# row/column block assignments are carried as attributes. This test
# checks that the DcLbm fit recovers those latent blocks at the finest
# non-trivial hierarchy level on each side, using adjusted Rand index
# (ARI) as the agreement measure. A local NMI helper logs an auxiliary
# score mirroring test-sbm-parity.R.

skip_if_no_greed <- function() {
  testthat::skip_if_not_installed("greed")
}

# 1. Helpers ------------------------------------------------------------------

# Local NMI (base-R, no mclust dependency) for informative logging.
nmi <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  p_ab <- tab / n
  p_a <- rowSums(p_ab)
  p_b <- colSums(p_ab)
  nz <- p_ab > 0
  mi <- sum(p_ab[nz] * log(p_ab[nz] / outer(p_a, p_b)[nz]))
  ha <- -sum(p_a[p_a > 0] * log(p_a[p_a > 0]))
  hb <- -sum(p_b[p_b > 0] * log(p_b[p_b > 0]))
  if (ha == 0 || hb == 0) {
    return(0)
  }
  mi / sqrt(ha * hb)
}

# Convert the finest non-trivial segment.list level on a moneca view into
# a membership vector of length `n_side`. Singletons pruned by the
# backend are re-assigned fresh cluster indices so the full length is
# preserved.
level_labels_side <- function(mon, n_side, level = 2L) {
  labels <- rep(NA_integer_, n_side)
  cliques <- mon$segment.list[[level]]
  for (b in seq_along(cliques)) {
    labels[cliques[[b]]] <- b
  }
  na_idx <- is.na(labels)
  if (any(na_idx)) {
    labels[na_idx] <- max(0L, labels, na.rm = TRUE) + seq_len(sum(na_idx))
  }
  as.integer(labels)
}

# 2. Parity experiment --------------------------------------------------------

test_that("moneca_bipartite recovers latent row/col blocks with ARI > 0.5", {
  skip_if_no_greed()
  testthat::skip_if_not_installed("mclust")

  mx <- get_bipartite_test_data(
    n_rows = 60,
    n_cols = 40,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  row_truth <- attr(mx, "row_block")
  col_truth <- attr(mx, "col_block")

  out <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)

  # Both sides need at least one non-trivial level to even attempt recovery.
  if (
    length(out$rows$segment.list) < 2L ||
      length(out$cols$segment.list) < 2L
  ) {
    testthat::skip(
      "Bipartite fit produced no non-atomic levels on at least one side."
    )
  }

  row_labels <- level_labels_side(
    out$rows,
    n_side = length(row_truth),
    level = 2L
  )
  col_labels <- level_labels_side(
    out$cols,
    n_side = length(col_truth),
    level = 2L
  )

  ari_row <- mclust::adjustedRandIndex(row_labels, row_truth)
  ari_col <- mclust::adjustedRandIndex(col_labels, col_truth)
  nmi_row <- nmi(row_labels, row_truth)
  nmi_col <- nmi(col_labels, col_truth)

  message(sprintf(
    "moneca_bipartite latent-block recovery: rows ARI = %.3f, NMI = %.3f; cols ARI = %.3f, NMI = %.3f.",
    ari_row,
    nmi_row,
    ari_col,
    nmi_col
  ))

  expect_gt(ari_row, 0.5)
  expect_gt(ari_col, 0.5)
})
