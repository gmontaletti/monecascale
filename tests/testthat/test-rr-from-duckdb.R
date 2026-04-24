# Tests for rr_from_duckdb()
# ==========================
#
# Contract: public entry point in R/rr_from_duckdb.R. Covers returned
# object shape, RR parity against a manually computed observed / expected
# baseline (the moneca::weight.matrix() arg surface is sidestepped in
# favour of the unambiguous O/E definition), sparse preservation, input
# polymorphism (data.frame, CSV file, DuckDB table name), bipartite mode,
# self-loop and min_count filters, and downstream composition with
# moneca_flow() / moneca_bipartite().

# 1. Structure and basic contract --------------------------------------------

test_that("returns a monecascale_rr object with correct structure", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  edges <- get_edge_list_fixture(n_nodes = 30, n_edges = 200, seed = 2026)
  res <- rr_from_duckdb(edges, weight_col = "weight")

  expect_s3_class(res, "monecascale_rr")
  expect_true(methods::is(res$rr, "dgCMatrix"))
  expect_true(methods::is(res$counts, "dgCMatrix"))
  expect_identical(dim(res$rr), dim(res$counts))
  expect_equal(nrow(res$rr), ncol(res$rr))
  expect_false(isTRUE(res$bipartite))
  expect_identical(length(res$margins$row), nrow(res$rr))
  expect_identical(length(res$margins$col), ncol(res$rr))
  expect_equal(res$margins$total, sum(edges$weight))
  expect_identical(res$node_names$row, res$node_names$col)
  expect_identical(res$n_edges, nrow(edges))
  expect_identical(res$n_nonzero, length(res$rr@x))
})

# 2. RR parity against a manual O/E baseline ---------------------------------

test_that("RR values match manually computed O / E on a small dense fixture", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")
  skip_if_not_installed("moneca")

  mx_full <- get_custom_test_data(n_classes = 8, seed = 2026)
  # generate_mobility_data() returns a matrix with the "Total" margin
  # appended as the last row / column; strip it for the dense baseline.
  n_core <- nrow(mx_full) - 1L
  core <- mx_full[seq_len(n_core), seq_len(n_core), drop = FALSE]

  # Convert the dense core to an edge list (long format).
  idx <- which(core > 0, arr.ind = TRUE)
  edges <- data.frame(
    row = rownames(core)[idx[, 1L]],
    col = colnames(core)[idx[, 2L]],
    weight = core[idx],
    stringsAsFactors = FALSE
  )

  res <- rr_from_duckdb(
    edges,
    weight_col = "weight",
    zero_diagonal = FALSE
  )

  # Manual baseline: E[i,j] = row_sum[i] * col_sum[j] / total,
  # RR = O / E for cells with O > 0.
  row_sum <- rowSums(core)
  col_sum <- colSums(core)
  total <- sum(core)
  E <- outer(row_sum, col_sum) / total
  expected_rr <- core / E

  rr_dense <- as.matrix(res$rr)

  # Align dimnames: the DuckDB path sorts node ids lexicographically.
  rr_dense <- rr_dense[rownames(core), colnames(core), drop = FALSE]

  nonzero <- core > 0
  expect_true(any(nonzero))
  expect_true(
    max(abs(rr_dense[nonzero] - expected_rr[nonzero])) < 1e-8
  )
})

# 3. Sparse preservation ------------------------------------------------------

test_that("sparse output is actually sparse on a 50-node 400-edge fixture", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  edges <- get_edge_list_fixture(n_nodes = 50, n_edges = 400, seed = 2026)
  res <- rr_from_duckdb(edges, weight_col = "weight")

  density <- length(res$rr@x) / prod(dim(res$rr))
  expect_lt(density, 0.5)
  expect_true(methods::is(res$rr, "dgCMatrix"))
})

# 4. File-path input (CSV) ----------------------------------------------------

test_that("CSV file-path input produces the same result as in-memory input", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  edges <- get_edge_list_fixture(n_nodes = 30, n_edges = 200, seed = 2026)
  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  utils::write.csv(edges, file = path, row.names = FALSE)

  res_mem <- rr_from_duckdb(edges, weight_col = "weight")
  res_csv <- rr_from_duckdb(
    path,
    row_col = "row",
    col_col = "col",
    weight_col = "weight"
  )

  expect_identical(dim(res_csv$rr), dim(res_mem$rr))
  expect_identical(res_csv$n_nonzero, res_mem$n_nonzero)
  expect_identical(res_csv$node_names$row, res_mem$node_names$row)
  # Values must match to floating point. Align by dimnames defensively.
  m_csv <- as.matrix(res_csv$rr)
  m_mem <- as.matrix(res_mem$rr)[rownames(m_csv), colnames(m_csv), drop = FALSE]
  expect_true(max(abs(m_csv - m_mem)) < 1e-12)
})

# 5. Pre-existing DuckDB connection ------------------------------------------

test_that("pre-existing DuckDB connection input matches transient connection", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  edges <- get_edge_list_fixture(n_nodes = 30, n_edges = 200, seed = 2026)

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(
    try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE),
    add = TRUE
  )

  DBI::dbWriteTable(con, "my_edges", edges)

  res_con <- rr_from_duckdb(
    "my_edges",
    con = con,
    weight_col = "weight"
  )
  res_mem <- rr_from_duckdb(edges, weight_col = "weight")

  expect_identical(dim(res_con$rr), dim(res_mem$rr))
  expect_identical(res_con$n_nonzero, res_mem$n_nonzero)
  m_con <- as.matrix(res_con$rr)
  m_mem <- as.matrix(res_mem$rr)[rownames(m_con), colnames(m_con), drop = FALSE]
  expect_true(max(abs(m_con - m_mem)) < 1e-12)
})

# 6. Bipartite mode -----------------------------------------------------------

test_that("bipartite mode produces rectangular output with disjoint names", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  edges <- get_edge_list_fixture(
    n_nodes = 30,
    n_edges = 200,
    seed = 2026,
    bipartite = TRUE
  )
  res <- rr_from_duckdb(
    edges,
    weight_col = "weight",
    bipartite = TRUE
  )

  expect_true(isTRUE(res$bipartite))
  expect_identical(nrow(res$rr), length(unique(edges$row)))
  expect_identical(ncol(res$rr), length(unique(edges$col)))
  # Disjoint namespaces: rows begin with "r", cols with "c".
  expect_true(all(grepl("^r", rownames(res$rr))))
  expect_true(all(grepl("^c", colnames(res$rr))))
  expect_length(intersect(rownames(res$rr), colnames(res$rr)), 0L)
})

# 7. zero_diagonal = TRUE drops self-loops ------------------------------------

test_that("zero_diagonal = TRUE drops self-loops", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  edges <- data.frame(
    row = c("a", "a", "b"),
    col = c("a", "b", "b"),
    weight = c(5, 3, 2),
    stringsAsFactors = FALSE
  )
  res <- rr_from_duckdb(edges, weight_col = "weight")

  rr_dense <- as.matrix(res$rr)
  expect_equal(rr_dense["a", "a"], 0)
  expect_equal(rr_dense["b", "b"], 0)
  expect_gt(rr_dense["a", "b"], 0)
})

# 8. zero_diagonal = FALSE keeps self-loops -----------------------------------

test_that("zero_diagonal = FALSE keeps self-loops", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  edges <- data.frame(
    row = c("a", "a", "b"),
    col = c("a", "b", "b"),
    weight = c(5, 3, 2),
    stringsAsFactors = FALSE
  )
  res <- rr_from_duckdb(
    edges,
    weight_col = "weight",
    zero_diagonal = FALSE
  )

  rr_dense <- as.matrix(res$rr)
  expect_gt(rr_dense["a", "a"], 0)
  expect_gt(rr_dense["b", "b"], 0)
})

# 9. min_count filters edges below threshold ----------------------------------

test_that("min_count filters edges below threshold", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")

  edges <- data.frame(
    row = c("a", "b", "c"),
    col = c("b", "c", "a"),
    weight = c(1, 2, 5),
    stringsAsFactors = FALSE
  )
  res <- rr_from_duckdb(
    edges,
    weight_col = "weight",
    min_count = 3
  )

  counts_dense <- as.matrix(res$counts)
  # Only the c -> a edge (weight 5) survives the filter.
  expect_equal(sum(counts_dense > 0), 1L)
  expect_equal(counts_dense["c", "a"], 5)
  expect_equal(counts_dense["a", "b"], 0)
  expect_equal(counts_dense["b", "c"], 0)
})

# 10. Downstream composition with moneca_flow() ------------------------------

test_that("$counts feeds moneca_flow() unchanged", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")
  skip_if_not_installed("igraph")

  edges <- get_edge_list_fixture(n_nodes = 30, n_edges = 200, seed = 2026)
  res <- rr_from_duckdb(edges, weight_col = "weight")
  dense_counts <- as.matrix(res$counts)

  fit <- moneca_flow(
    dense_counts,
    segment.levels = 2,
    seed = 2026,
    verbose = FALSE
  )

  expect_s3_class(fit, "moneca")
  expect_true(!is.null(fit$segment.list))
  expect_true(!is.null(fit$mat.list))
})

# 11. Downstream composition with moneca_bipartite() -------------------------

test_that("$rr feeds moneca_bipartite() unchanged", {
  skip_if_not_installed("duckdb")
  skip_if_not_installed("DBI")
  skip_if_not_installed("greed")

  edges <- get_edge_list_fixture(
    n_nodes = 30,
    n_edges = 200,
    seed = 2026,
    bipartite = TRUE
  )
  res <- rr_from_duckdb(
    edges,
    weight_col = "weight",
    bipartite = TRUE
  )

  fit <- suppressMessages(
    moneca_bipartite(res$rr, seed = 2026L, verbose = FALSE)
  )

  expect_s3_class(fit, "moneca_bipartite")
  expect_true(!is.null(fit$rows))
  expect_true(!is.null(fit$cols))
})
