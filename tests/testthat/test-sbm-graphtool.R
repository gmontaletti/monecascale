# Tests: moneca_sbm(backend = "graphtool") smoke test
# ====================================================
#
# Skipped unless both reticulate and a conda env with graph-tool are
# available.

skip_if_no_graphtool <- function() {
  testthat::skip_if_not_installed("reticulate")
  ok <- tryCatch(
    reticulate::py_module_available("graph_tool"),
    error = function(e) FALSE
  )
  if (!ok) {
    testthat::skip("graph-tool Python module not available")
  }
}

test_that("moneca_sbm(backend='graphtool') returns a valid moneca object", {
  skip_if_no_graphtool()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(
    mob,
    backend = "graphtool",
    seed = 1L,
    n_init = 2L,
    verbose = FALSE
  )
  expect_s3_class(out, "moneca")
  expect_equal(out$sbm_diagnostics$backend, "graphtool")
  expect_gte(length(out$segment.list), 1L)
})

test_that("moneca_sbm_install_graphtool is exported and self-documented", {
  expect_true(exists("moneca_sbm_install_graphtool"))
  expect_true(is.function(moneca_sbm_install_graphtool))
})
