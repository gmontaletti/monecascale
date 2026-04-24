# Internal helpers for rr_from_duckdb()
# ======================================
#
# SQL-side pipeline that takes a polymorphic edge-list input and produces
# a small R-side data.frame of non-zero (row_id, col_id, count, rr) tuples
# ready to be converted into sparse dgCMatrix output.
#
# Every helper here is internal. The public entry point lives in
# R/rr_from_duckdb.R. No helper touches a dense matrix representation.

# 1. Dependency guard ---------------------------------------------------------

#' @keywords internal
#' @noRd
.require_duckdb <- function() {
  have_duckdb <- requireNamespace("duckdb", quietly = TRUE)
  have_dbi <- requireNamespace("DBI", quietly = TRUE)
  if (!have_duckdb || !have_dbi) {
    stop(
      "rr_from_duckdb() requires the 'duckdb' and 'DBI' packages. ",
      "Install with: install.packages(c('duckdb', 'DBI')).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# 2. Transient connection -----------------------------------------------------

#' @keywords internal
#' @noRd
.open_transient_duckdb <- function() {
  DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
}

# 3. Unique view naming -------------------------------------------------------

#' @keywords internal
#' @noRd
.new_view_name <- function(prefix = "monecascale_v") {
  suffix <- paste(
    sample(c(letters, 0:9), 10L, replace = TRUE),
    collapse = ""
  )
  paste0(prefix, "_", suffix)
}

# 4. Quote identifier (double-quote; escape embedded double quotes) -----------

#' @keywords internal
#' @noRd
.quote_ident <- function(x) {
  paste0("\"", gsub("\"", "\"\"", x, fixed = TRUE), "\"")
}

# 5. Register / identify the source and canonicalise columns -----------------
#
# Accepts:
#   - data.frame / data.table: register via duckdb::duckdb_register() on a
#     temp name, then build a view that casts row/col to VARCHAR and
#     renames to (row_id, col_id, w).
#   - character path to .parquet: CREATE VIEW via read_parquet().
#   - character path to .csv / .csv.gz: CREATE VIEW via read_csv_auto().
#   - character table name: CREATE VIEW that selects / renames from it.
#
# Returns the name of the canonical view.

#' @keywords internal
#' @noRd
.register_edges <- function(con, edges, row_col, col_col, weight_col) {
  if (!is.character(row_col) || length(row_col) != 1L) {
    stop("`row_col` must be a single string.", call. = FALSE)
  }
  if (!is.character(col_col) || length(col_col) != 1L) {
    stop("`col_col` must be a single string.", call. = FALSE)
  }
  if (!is.character(weight_col) || length(weight_col) != 1L) {
    stop("`weight_col` must be a single string.", call. = FALSE)
  }

  canonical <- .new_view_name("edges_canon")

  if (is.data.frame(edges)) {
    missing_cols <- setdiff(c(row_col, col_col, weight_col), names(edges))
    if (length(missing_cols) > 0L) {
      stop(
        "Input data.frame is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    # duckdb_register expects a data.frame; downgrade data.table if needed.
    df <- as.data.frame(edges)
    raw_name <- .new_view_name("edges_raw")
    duckdb::duckdb_register(con, raw_name, df)
    sql_body <- sprintf(
      "SELECT CAST(%s AS VARCHAR) AS row_id, CAST(%s AS VARCHAR) AS col_id, CAST(%s AS DOUBLE) AS w FROM %s",
      .quote_ident(row_col),
      .quote_ident(col_col),
      .quote_ident(weight_col),
      .quote_ident(raw_name)
    )
  } else if (is.character(edges) && length(edges) == 1L) {
    is_parquet <- grepl("\\.parquet$", edges, ignore.case = TRUE)
    is_csv <- grepl("\\.csv(\\.gz)?$", edges, ignore.case = TRUE)
    if (is_parquet) {
      source_sql <- sprintf(
        "read_parquet(%s)",
        DBI::dbQuoteString(con, edges)
      )
    } else if (is_csv) {
      source_sql <- sprintf(
        "read_csv_auto(%s)",
        DBI::dbQuoteString(con, edges)
      )
    } else {
      # Assume a pre-existing DuckDB table / view name in the supplied con.
      source_sql <- .quote_ident(edges)
    }
    sql_body <- sprintf(
      "SELECT CAST(%s AS VARCHAR) AS row_id, CAST(%s AS VARCHAR) AS col_id, CAST(%s AS DOUBLE) AS w FROM %s",
      .quote_ident(row_col),
      .quote_ident(col_col),
      .quote_ident(weight_col),
      source_sql
    )
  } else {
    stop(
      "`edges` must be a data.frame / data.table, a file path, or a ",
      "character table name.",
      call. = FALSE
    )
  }

  create_sql <- sprintf(
    "CREATE OR REPLACE TEMP VIEW %s AS %s",
    .quote_ident(canonical),
    sql_body
  )
  DBI::dbExecute(con, create_sql)

  # Validate schema by probing one row. A zero-row input is allowed; the
  # probe needs only to confirm the column set resolves.
  probe_sql <- sprintf(
    "SELECT row_id, col_id, w FROM %s LIMIT 0",
    .quote_ident(canonical)
  )
  DBI::dbGetQuery(con, probe_sql)

  canonical
}

# 6. Small-cell reduction -----------------------------------------------------
#
# Drops rows and columns whose marginal sums are below `threshold`. The
# threshold is an inclusive lower bound on the sum of w grouped by row_id
# (resp. col_id). The resulting view is again (row_id, col_id, w).

#' @keywords internal
#' @noRd
.apply_scr <- function(con, view_name, threshold) {
  out_name <- .new_view_name("edges_scr")
  sql <- sprintf(
    paste(
      "CREATE OR REPLACE TEMP VIEW %s AS",
      "WITH rsum AS (",
      "  SELECT row_id, SUM(w) AS s FROM %s GROUP BY row_id",
      "),",
      "csum AS (",
      "  SELECT col_id, SUM(w) AS s FROM %s GROUP BY col_id",
      ")",
      "SELECT e.row_id, e.col_id, e.w",
      "FROM %s e",
      "JOIN rsum r ON r.row_id = e.row_id",
      "JOIN csum c ON c.col_id = e.col_id",
      "WHERE r.s >= %f AND c.s >= %f",
      sep = " "
    ),
    .quote_ident(out_name),
    .quote_ident(view_name),
    .quote_ident(view_name),
    .quote_ident(view_name),
    threshold,
    threshold
  )
  DBI::dbExecute(con, sql)
  out_name
}

# 7. Main RR SQL --------------------------------------------------------------
#
# Returns a list with:
#   - `edges`: data.frame of (row_id, col_id, count, rr) restricted to
#     edges that survive the min_count / zero_diagonal filters.
#   - `row_margins`: data.frame of (row_id, s) covering the full input
#     (i.e. before diagonal / min_count filtering), matching the
#     row_m CTE used in the RR denominator.
#   - `col_margins`: data.frame of (col_id, s), analogously.
#   - `total`: numeric scalar, SUM(w) over the full input; aligns with
#     the `tot` CTE.
#
# Margins mirror the RR denominator so R-side consumers can read
# observed marginals without re-aggregating; they intentionally include
# self-loop mass (when zero_diagonal = TRUE) and sub-threshold edges
# (when min_count > 0) so RR stays consistent with its expected values.

#' @keywords internal
#' @noRd
.compute_rr_sql <- function(
  con,
  view_name,
  zero_diagonal,
  min_count,
  bipartite
) {
  min_count <- as.numeric(min_count)
  extra_filters <- character(0)
  if (min_count > 0) {
    extra_filters <- c(
      extra_filters,
      sprintf("e.w >= %f", min_count)
    )
  }
  if (isTRUE(zero_diagonal) && !isTRUE(bipartite)) {
    extra_filters <- c(extra_filters, "e.row_id <> e.col_id")
  }
  where_clause <- if (length(extra_filters) > 0L) {
    paste0("WHERE e.w > 0 AND ", paste(extra_filters, collapse = " AND "))
  } else {
    "WHERE e.w > 0"
  }

  edges_sql <- sprintf(
    paste(
      "WITH row_m AS (SELECT row_id, SUM(w) AS s FROM %s GROUP BY row_id),",
      "col_m AS (SELECT col_id, SUM(w) AS s FROM %s GROUP BY col_id),",
      "tot   AS (SELECT SUM(w) AS total FROM %s)",
      "SELECT e.row_id, e.col_id, e.w AS count,",
      "       e.w * tot.total / (r.s * c.s) AS rr",
      "FROM %s e",
      "JOIN row_m r USING (row_id)",
      "JOIN col_m c USING (col_id)",
      "CROSS JOIN tot",
      "%s",
      sep = " "
    ),
    .quote_ident(view_name),
    .quote_ident(view_name),
    .quote_ident(view_name),
    .quote_ident(view_name),
    where_clause
  )
  edges <- DBI::dbGetQuery(con, edges_sql)

  row_margins <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT row_id, SUM(w) AS s FROM %s GROUP BY row_id",
      .quote_ident(view_name)
    )
  )
  col_margins <- DBI::dbGetQuery(
    con,
    sprintf(
      "SELECT col_id, SUM(w) AS s FROM %s GROUP BY col_id",
      .quote_ident(view_name)
    )
  )
  total <- as.numeric(
    DBI::dbGetQuery(
      con,
      sprintf("SELECT SUM(w) AS total FROM %s", .quote_ident(view_name))
    )$total
  )
  if (length(total) == 0L || is.na(total)) {
    total <- 0
  }

  list(
    edges = edges,
    row_margins = row_margins,
    col_margins = col_margins,
    total = total
  )
}

# 8. Convert the SQL result to sparse dgCMatrix output ------------------------

#' @keywords internal
#' @noRd
.result_to_sparse <- function(sql_result, bipartite) {
  if (
    !is.list(sql_result) ||
      !all(
        c("edges", "row_margins", "col_margins", "total") %in%
          names(sql_result)
      )
  ) {
    stop(
      "Internal error: SQL result does not have the expected shape.",
      call. = FALSE
    )
  }
  df <- sql_result$edges
  row_m_df <- sql_result$row_margins
  col_m_df <- sql_result$col_margins
  total <- as.numeric(sql_result$total)

  # Empty-result branch: still return a well-formed object.
  if (!nrow(df)) {
    empty <- Matrix::sparseMatrix(
      i = integer(0),
      j = integer(0),
      x = numeric(0),
      dims = c(0L, 0L)
    )
    return(
      list(
        counts = empty,
        rr = empty,
        margins = list(row = numeric(0), col = numeric(0), total = total),
        node_names = list(row = character(0), col = character(0))
      )
    )
  }

  row_ids <- as.character(df$row_id)
  col_ids <- as.character(df$col_id)
  counts_x <- as.numeric(df$count)
  rr_x <- as.numeric(df$rr)

  row_m_ids <- as.character(row_m_df$row_id)
  col_m_ids <- as.character(col_m_df$col_id)

  if (isTRUE(bipartite)) {
    # Disjoint namespaces. Use row / col margin tables as the level set
    # so the output matrix covers every node that contributed to the
    # RR denominator, even if all its edges were filtered out.
    row_levels <- sort(unique(row_m_ids))
    col_levels <- sort(unique(col_m_ids))
  } else {
    shared_levels <- sort(unique(c(row_m_ids, col_m_ids)))
    row_levels <- shared_levels
    col_levels <- shared_levels
  }

  row_idx <- match(row_ids, row_levels)
  col_idx <- match(col_ids, col_levels)
  dims <- c(length(row_levels), length(col_levels))
  dimnames_out <- list(row_levels, col_levels)

  counts <- Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = counts_x,
    dims = dims,
    dimnames = dimnames_out
  )
  rr <- Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = rr_x,
    dims = dims,
    dimnames = dimnames_out
  )
  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(as(counts, "generalMatrix"), "CsparseMatrix")
  }
  if (!inherits(rr, "dgCMatrix")) {
    rr <- as(as(rr, "generalMatrix"), "CsparseMatrix")
  }

  # Margins mirror the RR denominator (full input) rather than the
  # filtered `counts` matrix. Missing nodes (no edges) get zero mass.
  row_margin <- numeric(length(row_levels))
  row_margin[match(row_m_ids, row_levels)] <- as.numeric(row_m_df$s)
  col_margin <- numeric(length(col_levels))
  col_margin[match(col_m_ids, col_levels)] <- as.numeric(col_m_df$s)

  list(
    counts = counts,
    rr = rr,
    margins = list(row = row_margin, col = col_margin, total = total),
    node_names = list(row = row_levels, col = col_levels)
  )
}
