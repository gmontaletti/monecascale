# Tests: moneca_sbm() — core behavior and contract
# =================================================

# 1. Skip guard: greed is optional -------------------------------------------

skip_if_no_greed <- function() {
  testthat::skip_if_not_installed("greed")
}

# 2. Basic greed-backend invocation ------------------------------------------

test_that("moneca_sbm(backend='greed') returns a moneca-class object", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "greed", seed = 1L, verbose = FALSE)
  expect_s3_class(out, "moneca")
  expect_true(all(
    c(
      "segment.list",
      "mat.list",
      "small.cell.reduction",
      "margins_added",
      "density_reduction",
      "symmetric_method",
      "segment_metadata",
      "sbm_diagnostics"
    ) %in%
      names(out)
  ))
})

test_that("segment.list[[1]] is atomic and subsequent levels prune singletons", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "greed", seed = 1L)
  n_core <- nrow(mob) - 1L
  expect_length(out$segment.list[[1]], n_core)
  expect_equal(
    vapply(out$segment.list[[1]], length, integer(1)),
    rep(1L, n_core)
  )
  for (lvl in seq.int(2L, length(out$segment.list))) {
    cliques <- out$segment.list[[lvl]]
    if (length(cliques) > 0L) {
      expect_true(all(vapply(cliques, is.integer, logical(1))))
      expect_true(all(vapply(cliques, length, integer(1)) >= 2L))
    }
  }
})

test_that("mat.list length matches segment.list and mat.list[[1]] has margins", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "greed", seed = 1L)
  expect_equal(length(out$mat.list), length(out$segment.list))
  mx1 <- out$mat.list[[1]]
  n <- nrow(mx1)
  core <- mx1[1:(n - 1), 1:(n - 1)]
  expect_equal(unname(as.numeric(mx1[n, 1:(n - 1)])), unname(colSums(core)))
  expect_equal(unname(as.numeric(mx1[1:(n - 1), n])), unname(rowSums(core)))
})

test_that("sbm_diagnostics slot records backend + per-level K and MDL", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "greed", seed = 1L)
  diag <- out$sbm_diagnostics
  expect_equal(diag$backend, "greed")
  expect_true(is.integer(diag$n_blocks_per_level))
  expect_true(is.numeric(diag$mdl_per_level))
  expect_equal(
    length(diag$n_blocks_per_level),
    length(out$segment.list) - 1L
  )
})

# 3. Downstream consumer compatibility with moneca ---------------------------

test_that("moneca::segment.membership.dataframe works on monecascale output", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "greed", seed = 1L)
  df <- moneca::segment.membership.dataframe(out)
  expect_s3_class(df, "data.frame")
  expect_true(nrow(df) >= (nrow(mob) - 1L))
})

test_that("moneca::segment.quality works on monecascale output", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "greed", seed = 1L)
  q <- moneca::segment.quality(out)
  expect_true(is.list(q) || is.data.frame(q))
})

test_that("print.moneca runs on monecascale output", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "greed", seed = 1L)
  expect_output(print(out))
})

# 4. Parameter handling ------------------------------------------------------

test_that("segment.levels truncates backend hierarchy", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  full <- moneca_sbm(mob, backend = "greed", seed = 1L)
  trimmed <- moneca_sbm(
    mob,
    backend = "greed",
    seed = 1L,
    segment.levels = 2L
  )
  expect_lte(length(trimmed$segment.list), 3L)
  expect_lte(length(trimmed$segment.list), length(full$segment.list))
})

test_that("symmetric_method is honored and recorded", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out_sum <- moneca_sbm(
    mob,
    backend = "greed",
    seed = 1L,
    symmetric_method = "sum"
  )
  out_min <- moneca_sbm(
    mob,
    backend = "greed",
    seed = 1L,
    symmetric_method = "min"
  )
  expect_equal(out_sum$symmetric_method, "sum")
  expect_equal(out_min$symmetric_method, "min")
})

test_that("isolates = TRUE attaches isolates_summary", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "greed", seed = 1L, isolates = TRUE)
  expect_true("isolates_summary" %in% names(out))
  expect_true(is.data.frame(out$isolates_summary$membership))
  expect_true(is.matrix(out$isolates_summary$mobility_matrix))
})

# 5. Input validation --------------------------------------------------------

test_that("rectangular input is rejected with pointer to D1", {
  skip_if_no_greed()
  rect <- matrix(rpois(10 * 6, 2), 10, 6)
  expect_error(
    moneca_sbm(rect, backend = "greed"),
    regexp = "square|bipartite",
    ignore.case = TRUE
  )
})

test_that("invalid backend triggers match.arg error", {
  expect_error(moneca_sbm(get_test_data("small"), backend = "nope"))
})

# 6. Auto-backend dispatch ---------------------------------------------------

test_that("backend='auto' resolves to greed on modest-size input", {
  skip_if_no_greed()
  mob <- get_test_data(size = "medium", seed = 42)
  out <- moneca_sbm(mob, backend = "auto", seed = 1L)
  expect_equal(out$sbm_diagnostics$backend, "greed")
})
