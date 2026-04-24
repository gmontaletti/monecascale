# Tests: segment.levels = "auto" round-trip on wrappers
# =====================================================
#
# Exercises the `segment.levels = "auto"` path on moneca_sbm() and
# moneca_bipartite(). Asserts that each wrapper trims its output to the
# auto-selected level and attaches the picker result under $auto_level.
# Also guards against regressions on the explicit-integer path, which
# must not attach $auto_level.

# 1. Skip guard: greed is optional -------------------------------------------

skip_if_no_greed <- function() {
  testthat::skip_if_not_installed("greed")
}

# 2. SBM wrapper, default auto_method = "mdl" --------------------------------

test_that("moneca_sbm(segment.levels='auto') trims and attaches $auto_level", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(
    mob,
    segment.levels = "auto",
    seed = 2026,
    verbose = FALSE
  )

  expect_s3_class(out, "moneca")
  expect_s3_class(out$auto_level, "auto_segment_levels")
  expect_equal(length(out$segment.list), out$auto_level$level)
})

# 3. SBM wrapper, auto_method = "mi_plateau" ---------------------------------

test_that("moneca_sbm(segment.levels='auto', auto_method='mi_plateau') records the method", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(
    mob,
    segment.levels = "auto",
    auto_method = "mi_plateau",
    seed = 2026,
    verbose = FALSE
  )

  expect_s3_class(out$auto_level, "auto_segment_levels")
  expect_equal(out$auto_level$method, "mi_plateau")
})

# 4. Bipartite wrapper -------------------------------------------------------

test_that("moneca_bipartite(segment.levels='auto') trims and attaches $auto_level per side", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  out <- moneca_bipartite(
    mx,
    segment.levels = "auto",
    seed = 2026,
    verbose = FALSE
  )

  expect_s3_class(out$rows$auto_level, "auto_segment_levels")
  expect_s3_class(out$cols$auto_level, "auto_segment_levels")
  expect_equal(
    length(out$rows$segment.list),
    out$rows$auto_level$level
  )
  expect_equal(
    length(out$cols$segment.list),
    out$cols$auto_level$level
  )
})

# 5. Regression guard: explicit integer does not attach $auto_level ---------

test_that("moneca_sbm(segment.levels=3L) does not attach $auto_level", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(
    mob,
    segment.levels = 3L,
    seed = 2026,
    verbose = FALSE
  )

  expect_s3_class(out, "moneca")
  expect_null(out$auto_level)
})

# 6. Flow wrapper, default auto_method = "mdl" -------------------------------

test_that("moneca_flow(segment.levels = 'auto') attaches $auto_level", {
  skip_if_not_installed("igraph")
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  out <- moneca_flow(mx, depth = 3L, segment.levels = "auto", seed = 2026)
  expect_s3_class(out$auto_level, "auto_segment_levels")
  expect_equal(length(out$segment.list), out$auto_level$level)
})
