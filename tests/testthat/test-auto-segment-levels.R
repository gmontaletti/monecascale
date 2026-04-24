# Tests: auto_segment_levels() — standalone function contract
# ============================================================
#
# Covers dispatcher guards, backend detection, both criteria (MDL and
# MI-plateau), return-mode behavior, bipartite top-level dispatch and
# per-side detached dispatch, both elbow variants, and the print method.

# 1. Skip guard: greed is optional -------------------------------------------

skip_if_no_greed <- function() {
  testthat::skip_if_not_installed("greed")
}

# 2. Dispatcher error handling ------------------------------------------------

test_that("auto_segment_levels() rejects non-moneca input", {
  expect_error(
    auto_segment_levels(list(a = 1)),
    regexp = "moneca"
  )
})

# 3. Backend detection on a fast-like object ---------------------------------

test_that("auto_segment_levels() handles a trivial fast-like object", {
  skip_if_not_installed("igraph")
  fake_fast_like <- list(
    segment.list = list(as.list(1:5)),
    mat.list = list(matrix(0, 6, 6)),
    small.cell.reduction = 0,
    margins_added = TRUE,
    density_reduction = NULL,
    symmetric_method = "sum"
  )
  class(fake_fast_like) <- "moneca"
  res <- auto_segment_levels(fake_fast_like, method = "mdl")
  expect_s3_class(res, "auto_segment_levels")
  expect_equal(res$backend, "fast")
  expect_equal(res$level, 1L)
})

# 4. MDL criterion on SBM fit -------------------------------------------------

test_that("method='mdl' returns a well-formed auto_segment_levels on SBM fit", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  fit <- moneca_sbm(mob, backend = "greed", seed = 1L, verbose = FALSE)
  res <- auto_segment_levels(fit, method = "mdl")

  expect_s3_class(res, "auto_segment_levels")
  expect_true(is.integer(res$level))
  expect_gte(res$level, 1L)
  expect_s3_class(res$diagnostics, "data.frame")
  expect_true("mdl" %in% colnames(res$diagnostics))
  expect_equal(res$backend, "sbm")
})

# 5. MI-plateau criterion on SBM fit -----------------------------------------

test_that("method='mi_plateau' returns a well-formed auto_segment_levels on SBM fit", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  fit <- moneca_sbm(mob, backend = "greed", seed = 1L, verbose = FALSE)
  res <- auto_segment_levels(fit, method = "mi_plateau")

  expect_s3_class(res, "auto_segment_levels")
  expect_true("mi_to_next" %in% colnames(res$diagnostics))
  expect_gte(res$level, 1L)
})

# 6. Return-mode ---------------------------------------------------------------

test_that("return='level' returns a bare integer scalar", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  fit <- moneca_sbm(mob, backend = "greed", seed = 1L, verbose = FALSE)
  lvl <- auto_segment_levels(fit, return = "level")

  expect_true(is.integer(lvl))
  expect_length(lvl, 1L)
  expect_gte(lvl, 1L)
})

# 7. Bipartite top-level dispatch --------------------------------------------

test_that("auto_segment_levels() on moneca_bipartite returns a per-side list", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  bip_fit <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)
  res <- auto_segment_levels(bip_fit)

  expect_true(is.list(res))
  expect_true(all(c("rows", "cols") %in% names(res)))
  expect_s3_class(res$rows, "auto_segment_levels")
  expect_s3_class(res$cols, "auto_segment_levels")
})

# 8. Detached bipartite side dispatch ----------------------------------------

test_that("auto_segment_levels() runs on a detached bipartite side", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  bip_fit <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)

  res_rows <- auto_segment_levels(bip_fit$rows)
  expect_s3_class(res_rows, "auto_segment_levels")
  expect_equal(res_rows$backend, "bipartite_rows")

  res_cols <- auto_segment_levels(bip_fit$cols)
  expect_s3_class(res_cols, "auto_segment_levels")
  expect_equal(res_cols$backend, "bipartite_cols")
})

# 9. Elbow variants -----------------------------------------------------------

test_that("both mdl_elbow modes produce valid auto_segment_levels objects", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  fit <- moneca_sbm(mob, backend = "greed", seed = 1L, verbose = FALSE)

  res_kneedle <- auto_segment_levels(
    fit,
    method = "mdl",
    mdl_elbow = "kneedle"
  )
  res_msd <- auto_segment_levels(
    fit,
    method = "mdl",
    mdl_elbow = "max_second_diff"
  )

  expect_s3_class(res_kneedle, "auto_segment_levels")
  expect_s3_class(res_msd, "auto_segment_levels")
})

# 10. Print method ------------------------------------------------------------

test_that("print.auto_segment_levels emits a recognizable header", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  fit <- moneca_sbm(mob, backend = "greed", seed = 1L, verbose = FALSE)
  res <- auto_segment_levels(fit, method = "mdl")

  expect_output(print(res), regexp = "auto_segment_levels|level")
})

# 11. Flow backend ------------------------------------------------------------

test_that("auto_segment_levels works on flow backend", {
  skip_if_not_installed("igraph")
  mx <- get_flow_test_data(n = 30, n_modules = 3, seed = 2026)
  fit <- moneca_flow(mx, depth = 3L, nb_trials = 10L, seed = 2026)
  res <- auto_segment_levels(fit, method = "mdl")
  expect_s3_class(res, "auto_segment_levels")
  expect_equal(res$backend, "flow")
  expect_true("codelength" %in% colnames(res$diagnostics))
  expect_true(res$level >= 1L && res$level <= length(fit$segment.list))
})

# 12. moneca_fast() backend --------------------------------------------------

test_that("auto_segment_levels supports moneca_fast() with both criteria", {
  skip_if_not_installed("moneca")
  skip_if_not_installed("igraph")

  set.seed(2026)
  mx <- get_custom_test_data(n_classes = 12, seed = 2026)
  fit <- moneca::moneca_fast(mx, segment.levels = 4)

  res_mdl <- auto_segment_levels(fit, method = "mdl")
  expect_s3_class(res_mdl, "auto_segment_levels")
  expect_equal(res_mdl$backend, "fast")
  expect_true("modularity" %in% colnames(res_mdl$diagnostics))
  expect_true(res_mdl$level >= 1L)

  res_mi <- auto_segment_levels(fit, method = "mi_plateau")
  expect_equal(res_mi$backend, "fast")
  expect_true(res_mi$level >= 1L)
})

test_that("fast error gate is gone", {
  skip_if_not_installed("moneca")
  skip_if_not_installed("igraph")

  set.seed(2026)
  mx <- get_custom_test_data(n_classes = 8, seed = 2026)
  fit <- moneca::moneca_fast(mx, segment.levels = 3)
  expect_error(
    auto_segment_levels(fit, method = "mdl"),
    regexp = NA
  )
})
