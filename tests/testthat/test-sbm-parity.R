# Tests: empirical parity between monecascale::moneca_sbm() and moneca::moneca_fast()
# ===================================================================================
#
# The theoretical bridge (RR = O/E is the DC-SBM residual under identity
# block interaction) predicts that moneca_sbm(backend="greed") should
# recover segments similar to moneca::moneca_fast() on a well-structured
# mobility matrix. This test records the empirical agreement but does NOT
# gate acceptance on a strict threshold — the parity is an empirical
# question whose answer is itself a scientific finding reported in the
# vignette.

skip_if_no_greed <- function() {
  testthat::skip_if_not_installed("greed")
}

# 1. Helpers ------------------------------------------------------------------

ari <- function(a, b) {
  stopifnot(length(a) == length(b))
  tab <- table(a, b)
  n <- sum(tab)
  sum_comb <- function(x) sum(choose(x, 2))
  a_sum <- sum_comb(rowSums(tab))
  b_sum <- sum_comb(colSums(tab))
  tab_sum <- sum_comb(as.vector(tab))
  expected <- a_sum * b_sum / choose(n, 2)
  max_idx <- (a_sum + b_sum) / 2
  if (max_idx == expected) {
    return(1)
  }
  (tab_sum - expected) / (max_idx - expected)
}

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

level_labels <- function(mon, level) {
  n_core <- nrow(mon$mat.list[[1]]) - 1L
  labels <- rep(NA_integer_, n_core)
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

test_that("moneca_sbm vs moneca::moneca_fast: best-match ARI / NMI", {
  skip_if_no_greed()

  mob <- get_custom_test_data(
    n_classes = 20,
    n_total = 5000,
    immobility_strength = 0.5,
    class_clustering = 0.75,
    noise_level = 0.05,
    seed = 2026
  )

  fast <- moneca::moneca_fast(mob, segment.levels = 3, progress = FALSE)
  sbm <- moneca_sbm(mob, backend = "greed", seed = 2026, verbose = FALSE)

  fast_levels <- seq.int(2L, length(fast$segment.list))
  sbm_levels <- seq.int(2L, length(sbm$segment.list))
  if (length(fast_levels) == 0L || length(sbm_levels) == 0L) {
    testthat::skip(
      "Either moneca_fast or moneca_sbm produced no non-atomic levels."
    )
  }

  scores <- expand.grid(fast_lvl = fast_levels, sbm_lvl = sbm_levels)
  scores$ari <- NA_real_
  scores$nmi <- NA_real_
  for (i in seq_len(nrow(scores))) {
    lf <- level_labels(fast, scores$fast_lvl[i])
    ls <- level_labels(sbm, scores$sbm_lvl[i])
    scores$ari[i] <- ari(lf, ls)
    scores$nmi[i] <- nmi(lf, ls)
  }

  best <- scores[which.max(scores$ari), , drop = FALSE]

  message(sprintf(
    "moneca_sbm vs moneca_fast parity: ARI = %.3f, NMI = %.3f (fast lvl %d vs sbm lvl %d).",
    best$ari,
    best$nmi,
    best$fast_lvl,
    best$sbm_lvl
  ))

  expect_gt(best$ari, 0)
  expect_gt(best$nmi, 0)
})
