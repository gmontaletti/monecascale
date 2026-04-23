# Tests: moneca_bipartite() — core behavior and contract
# =======================================================
#
# The bipartite entry point must produce:
#   (a) a moneca_bipartite S3 object with $rows, $cols, $bipartite_diagnostics;
#   (b) two moneca-class objects (one per view) each carrying the standard
#       moneca slots and a square mat.list at every level;
#   (c) downstream compatibility with moneca's segment.* and plot stack.
#
# The "every mat.list is square" check is the load-bearing invariant that
# justifies exposing a rectangular fit to the square-only moneca consumer
# stack.

# 1. Skip guard: greed is optional -------------------------------------------

skip_if_no_greed <- function() {
  testthat::skip_if_not_installed("greed")
}

# 2. Top-level class and slots ------------------------------------------------

test_that("moneca_bipartite() returns a moneca_bipartite with $rows, $cols, diagnostics", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  out <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)

  expect_s3_class(out, "moneca_bipartite")
  expect_true(all(
    c("rows", "cols", "bipartite_diagnostics") %in% names(out)
  ))
  expect_false(inherits(out, "moneca"))
})

# 3. Each side exposes the full moneca output contract -----------------------

test_that("$rows and $cols are moneca-class objects with the standard slots", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  out <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)

  required <- c(
    "segment.list",
    "mat.list",
    "segment_metadata",
    "small.cell.reduction",
    "margins_added",
    "symmetric_method"
  )

  expect_s3_class(out$rows, "moneca")
  expect_s3_class(out$cols, "moneca")
  expect_true(all(required %in% names(out$rows)))
  expect_true(all(required %in% names(out$cols)))
})

# 4. Every mat.list on both sides is square ----------------------------------

test_that("every mat.list[[l]] on $rows and $cols is square", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  out <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)

  for (lvl in seq_along(out$rows$mat.list)) {
    m <- out$rows$mat.list[[lvl]]
    expect_equal(
      nrow(m),
      ncol(m),
      info = sprintf("rows mat.list level %d", lvl)
    )
  }
  for (lvl in seq_along(out$cols$mat.list)) {
    m <- out$cols$mat.list[[lvl]]
    expect_equal(
      nrow(m),
      ncol(m),
      info = sprintf("cols mat.list level %d", lvl)
    )
  }
})

# 5. Downstream moneca consumers run without error ---------------------------

test_that("moneca::segment.quality and segment.membership run on $rows and $cols", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  out <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)

  expect_no_error(moneca::segment.quality(out$rows))
  expect_no_error(moneca::segment.membership(out$rows))
  expect_no_error(moneca::segment.quality(out$cols))
  expect_no_error(moneca::segment.membership(out$cols))
})

# 6. moneca plotting builds a ggplot on the row view -------------------------

test_that("moneca::plot_moneca_ggraph constructs on $rows at level 2", {
  skip_if_no_greed()
  skip_if_not_installed("ggraph")
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  out <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)
  # Guard against an output that produced only the trivial level 1.
  if (length(out$rows$segment.list) < 2L) {
    testthat::skip("No non-atomic row level available for plotting.")
  }

  p <- moneca::plot_moneca_ggraph(out$rows, level = 2)
  expect_s3_class(p, "ggplot")
})

# 7. has_margins = "auto" is idempotent on pre-appended input ----------------

test_that("has_margins='auto' on pre-appended margins matches raw rectangular input", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )

  # Build the moneca-style margin-appended counterpart explicitly.
  rs <- rowSums(mx)
  cs <- colSums(mx)
  gt <- sum(mx)
  mx_aug <- rbind(cbind(mx, rs), c(cs, gt))
  rownames(mx_aug) <- c(rownames(mx), "Total")
  colnames(mx_aug) <- c(colnames(mx), "Total")

  out_raw <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)
  out_aug <- moneca_bipartite(
    mx_aug,
    seed = 2026,
    has_margins = "auto",
    verbose = FALSE
  )

  # Dimensionality
  expect_equal(
    out_raw$bipartite_diagnostics$n_row,
    out_aug$bipartite_diagnostics$n_row
  )
  expect_equal(
    out_raw$bipartite_diagnostics$n_col,
    out_aug$bipartite_diagnostics$n_col
  )

  # segment.list structure per side, level by level
  expect_equal(
    length(out_raw$rows$segment.list),
    length(out_aug$rows$segment.list)
  )
  expect_equal(
    length(out_raw$cols$segment.list),
    length(out_aug$cols$segment.list)
  )
  expect_equal(out_raw$rows$segment.list, out_aug$rows$segment.list)
  expect_equal(out_raw$cols$segment.list, out_aug$cols$segment.list)

  # Per-level memberships in diagnostics (deterministic under fixed seed)
  expect_equal(
    out_raw$bipartite_diagnostics$n_blocks_row_per_level,
    out_aug$bipartite_diagnostics$n_blocks_row_per_level
  )
  expect_equal(
    out_raw$bipartite_diagnostics$n_blocks_col_per_level,
    out_aug$bipartite_diagnostics$n_blocks_col_per_level
  )
})

# 8. Diagnostics: per-level block-interaction matrix shape -------------------

test_that("block_interaction_matrix_per_level matches Krow x Kcol per level", {
  skip_if_no_greed()
  mx <- get_bipartite_test_data(
    n_rows = 30,
    n_cols = 20,
    k_row = 3,
    k_col = 2,
    seed = 2026
  )
  out <- moneca_bipartite(mx, seed = 2026, verbose = FALSE)
  diag <- out$bipartite_diagnostics

  expect_equal(
    length(diag$block_interaction_matrix_per_level),
    length(diag$n_blocks_row_per_level)
  )
  for (lvl in seq_along(diag$block_interaction_matrix_per_level)) {
    bim <- diag$block_interaction_matrix_per_level[[lvl]]
    expect_equal(
      nrow(bim),
      as.integer(diag$n_blocks_row_per_level[lvl]),
      info = sprintf("block_interaction nrow at level %d", lvl)
    )
    expect_equal(
      ncol(bim),
      as.integer(diag$n_blocks_col_per_level[lvl]),
      info = sprintf("block_interaction ncol at level %d", lvl)
    )
  }
})
