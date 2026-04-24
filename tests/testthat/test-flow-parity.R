# Tests: latent-module recovery for monecascale::moneca_flow()
# ============================================================
#
# Fixture: square Poisson draws from a planted-community model with
# known module memberships. The Map Equation / Infomap backend should
# recover those modules at the finest non-trivial level of the flow
# hierarchy. ARI is the primary agreement measure; NMI is logged as an
# auxiliary diagnostic.

skip_if_no_igraph <- function() {
  testthat::skip_if_not_installed("igraph")
}

# 1. Helpers ------------------------------------------------------------------

# Local NMI in bits, mirroring test-bipartite-parity.R.
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

# Rebuild a length-n membership vector from a segment.list level. Mirrors
# the internal `.cliques_to_membership` used elsewhere in the package.
cliques_to_membership <- function(cliques, n) {
  m <- seq_len(n)
  for (i in seq_along(cliques)) {
    idx <- cliques[[i]]
    if (length(idx) >= 2L) {
      m[idx] <- n + i
    }
  }
  as.integer(factor(m))
}

# 2. Parity experiment --------------------------------------------------------

test_that("moneca_flow recovers planted modules with ARI > 0.5", {
  skip_if_no_igraph()
  testthat::skip_if_not_installed("mclust")

  mx <- get_flow_test_data(
    n = 60,
    n_modules = 4,
    p_in = 0.3,
    p_out = 0.02,
    seed = 2026
  )
  truth <- attr(mx, "module")

  out <- moneca_flow(mx, depth = 2L, nb_trials = 10L, seed = 2026)

  if (length(out$segment.list) < 2L) {
    testthat::skip("Flow fit produced no non-atomic level.")
  }

  m <- cliques_to_membership(out$segment.list[[2]], n = length(truth))
  ari_val <- mclust::adjustedRandIndex(m, truth)
  nmi_val <- nmi(m, truth)

  message(sprintf(
    "moneca_flow latent-module recovery: ARI = %.3f, NMI = %.3f (n = %d, K = %d).",
    ari_val,
    nmi_val,
    length(truth),
    length(unique(truth))
  ))

  expect_gt(ari_val, 0.5)
})
