# Tests: moneca_flow() — core behavior and contract
# ==================================================
#
# Covers the output contract for the Infomap / Map Equation backend:
# class, diagnostics, square mat.list, segment.list length for variable
# depth, downstream moneca consumer and plot compatibility, margin
# auto-detection, undirected mode, and the flat depth = 1L edge case.

# 1. Skip guard: igraph is a moneca dependency but guard anyway --------------

skip_if_no_igraph <- function() {
  testthat::skip_if_not_installed("igraph")
}

# 2. Class ---------------------------------------------------------------------

test_that("moneca_flow() returns a plain moneca-class object", {
  skip_if_no_igraph()
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  out <- moneca_flow(mx, depth = 2L, nb_trials = 5L, seed = 1L)

  expect_s3_class(out, "moneca")
  expect_false(inherits(out, "moneca_bipartite"))
})

# 3. Diagnostics slot ---------------------------------------------------------

test_that("$flow_diagnostics carries the documented fields", {
  skip_if_no_igraph()
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  out <- moneca_flow(mx, depth = 2L, nb_trials = 5L, seed = 1L)

  diag <- out$flow_diagnostics
  expect_true(is.list(diag))
  expect_equal(diag$backend, "igraph")
  expect_equal(diag$model, "infomap")
  expect_true(all(
    c(
      "depth",
      "n_blocks_per_level",
      "codelength_per_level",
      "mdl_per_level",
      "n_init",
      "directed",
      "node_names"
    ) %in%
      names(diag)
  ))
  # mdl_per_level is an alias of codelength_per_level
  expect_equal(diag$mdl_per_level, diag$codelength_per_level)
})

# 4. Square mat.list at every level -------------------------------------------

test_that("every mat.list[[l]] is square", {
  skip_if_no_igraph()
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  out <- moneca_flow(mx, depth = 2L, nb_trials = 5L, seed = 1L)

  for (lvl in seq_along(out$mat.list)) {
    m <- out$mat.list[[lvl]]
    expect_equal(
      nrow(m),
      ncol(m),
      info = sprintf("mat.list level %d", lvl)
    )
  }
})

# 5. segment.list length tracks depth -----------------------------------------

test_that("length(segment.list) == depth + 1 at depth 2L and 3L", {
  skip_if_no_igraph()
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)

  out2 <- moneca_flow(mx, depth = 2L, nb_trials = 5L, seed = 1L)
  out3 <- moneca_flow(mx, depth = 3L, nb_trials = 5L, seed = 1L)

  expect_equal(length(out2$segment.list), 3L)
  expect_equal(length(out3$segment.list), 4L)
  expect_equal(length(out2$flow_diagnostics$codelength_per_level), 2L)
  expect_equal(length(out3$flow_diagnostics$codelength_per_level), 3L)
})

# 6. Downstream moneca consumers run without error ---------------------------

test_that("moneca::segment.quality and segment.membership run on flow output", {
  skip_if_no_igraph()
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  out <- moneca_flow(mx, depth = 2L, nb_trials = 5L, seed = 1L)

  expect_no_error(moneca::segment.quality(out))
  expect_no_error(moneca::segment.membership(out))
})

# 7. Plot compatibility --------------------------------------------------------

test_that("moneca::plot_moneca_ggraph constructs on flow output at level 2", {
  skip_if_no_igraph()
  skip_if_not_installed("ggraph")
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  out <- moneca_flow(mx, depth = 2L, nb_trials = 5L, seed = 1L)

  if (length(out$segment.list) < 2L) {
    testthat::skip("Flow fit produced no non-atomic level for plotting.")
  }

  p <- moneca::plot_moneca_ggraph(out, level = 2)
  expect_s3_class(p, "ggplot")
})

# 8. has_margins = "auto" is idempotent on pre-appended margins ---------------

test_that("has_margins='auto' on pre-appended margins matches raw input", {
  skip_if_no_igraph()
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)

  rs <- rowSums(mx)
  cs <- colSums(mx)
  gt <- sum(mx)
  mx_aug <- rbind(cbind(mx, rs), c(cs, gt))
  rownames(mx_aug) <- c(rownames(mx), "Total")
  colnames(mx_aug) <- c(colnames(mx), "Total")

  out_raw <- moneca_flow(mx, depth = 2L, nb_trials = 5L, seed = 1L)
  out_aug <- moneca_flow(
    mx_aug,
    depth = 2L,
    nb_trials = 5L,
    seed = 1L,
    has_margins = "auto"
  )

  # Level-1 (post-atomic) membership must match under a fixed seed.
  expect_equal(out_raw$segment.list[[2]], out_aug$segment.list[[2]])
})

# 9. Undirected mode runs and produces a valid output -------------------------

test_that("moneca_flow(directed = FALSE) runs and yields a valid moneca", {
  skip_if_no_igraph()
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  out <- moneca_flow(
    mx,
    depth = 2L,
    nb_trials = 5L,
    directed = FALSE,
    seed = 1L
  )

  expect_s3_class(out, "moneca")
  expect_false(out$flow_diagnostics$directed)
  expect_equal(length(out$segment.list), 3L)
  for (lvl in seq_along(out$mat.list)) {
    expect_equal(nrow(out$mat.list[[lvl]]), ncol(out$mat.list[[lvl]]))
  }
})

# 10. Flat depth = 1L ---------------------------------------------------------

test_that("depth = 1L produces a flat hierarchy with one backend level", {
  skip_if_no_igraph()
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  out <- moneca_flow(mx, depth = 1L, nb_trials = 5L, seed = 1L)

  expect_equal(length(out$segment.list), 2L)
  expect_equal(length(out$flow_diagnostics$codelength_per_level), 1L)
  expect_equal(out$flow_diagnostics$depth, 1L)
})
