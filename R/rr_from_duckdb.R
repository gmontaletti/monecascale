# rr_from_duckdb(): out-of-core sparse relative-risk builder
# ==========================================================
#
# Builds a sparse relative-risk matrix and the companion raw-count matrix
# directly from an edge list, pushing the aggregation into DuckDB SQL. The
# output never materialises as a dense matrix, so memory scales with the
# number of non-zero cells rather than with n_rows * n_cols. Intended for
# person x employer / worker x occupation mobility at 10^5-10^6 node scale
# where `moneca::weight.matrix()` cannot fit in R memory.

#' @title Build a Sparse Relative-Risk Matrix Out-of-Core via DuckDB
#'
#' @description
#' Aggregates an edge list into a sparse relative-risk (`RR = O / E`)
#' matrix plus the companion raw-count matrix, using DuckDB as the
#' aggregation engine. The R side never allocates a dense
#' `n_rows x n_cols` matrix, so the function scales to person x employer
#' or worker x occupation mobility tables that would OOM under
#' `moneca::weight.matrix()`.
#'
#' @details
#' The RR aggregation is pushed into a single SQL query with the template
#'
#' ```sql
#' WITH row_m AS (SELECT row_id, SUM(w) AS s FROM edges GROUP BY row_id),
#'      col_m AS (SELECT col_id, SUM(w) AS s FROM edges GROUP BY col_id),
#'      tot   AS (SELECT SUM(w) AS total FROM edges)
#' SELECT e.row_id, e.col_id, e.w AS count,
#'        e.w * tot.total / (r.s * c.s) AS rr
#' FROM edges e
#' JOIN row_m r USING (row_id)
#' JOIN col_m c USING (col_id)
#' CROSS JOIN tot
#' WHERE e.w > 0;
#' ```
#'
#' The join-and-divide pattern only emits rows whose edge weight is
#' non-zero, so sparsity is preserved end to end: the result set has one
#' row per non-zero input edge, and the R-side `Matrix::sparseMatrix()`
#' call materialises it directly as a `dgCMatrix`. When
#' `bipartite = FALSE`, row and column IDs are unioned into a single
#' factor so that both output matrices are square with identical
#' dimnames. When `bipartite = TRUE`, row and column namespaces stay
#' disjoint and the output is rectangular.
#'
#' Preprocessing is applied in order: optional small-cell reduction
#' (drops rows/cols whose marginal sum falls below a threshold), optional
#' `min_count` filter (drops individual edges), optional zero-diagonal
#' filter (only active when `bipartite = FALSE`). All three filters are
#' applied in SQL.
#'
#' @param edges One of: (a) an in-memory `data.frame` / `data.table` with
#'   columns named by `row_col`, `col_col`, `weight_col`; (b) a single
#'   character file path ending in `.parquet`, `.csv`, or `.csv.gz`; or
#'   (c) a character scalar naming a pre-existing DuckDB table / view on
#'   the supplied `con`.
#' @param row_col,col_col,weight_col Column names inside the source.
#'   Defaults are `"row"`, `"col"`, `"weight"`.
#' @param con Optional `DBIConnection` to DuckDB. Required when `edges`
#'   is a bare table name. When `NULL`, a transient in-memory DuckDB
#'   connection is opened for the call and closed on exit.
#' @param bipartite Logical. If `FALSE` (default), row and column
#'   identifiers live in the same namespace and are unioned into a
#'   shared factor, yielding a square output with identical row and
#'   column names. If `TRUE`, identifiers are disjoint and the output is
#'   rectangular with separate row / column names.
#' @param min_count Numeric. Drop edges with `w < min_count` before the
#'   RR computation. Defaults to `0` (no filter).
#' @param small_cell_reduction Numeric. Drop rows and columns whose
#'   marginal sum is below this threshold before RR computation. Mirrors
#'   the `small.cell.reduction` argument of
#'   `moneca::moneca_fast()`. Defaults to `0` (no filter).
#' @param zero_diagonal Logical. If `TRUE` (default) and
#'   `bipartite = FALSE`, self-loops (`row_id == col_id`) are dropped
#'   before the RR computation. Ignored when `bipartite = TRUE`.
#' @param verbose Logical. If `TRUE`, message the number of non-zero RR
#'   cells returned.
#' @param ... Reserved for forward compatibility.
#'
#' @return A list of class `"monecascale_rr"` with elements
#'   * `rr` - sparse `dgCMatrix` of observed / expected; one entry per
#'     non-zero input cell.
#'   * `counts` - sparse `dgCMatrix` of raw observed counts; same
#'     sparsity pattern as `rr`.
#'   * `margins` - list with numeric vectors `row` and `col`, plus the
#'     scalar `total`.
#'   * `node_names` - list with character vectors `row` and `col`.
#'     Identical when `bipartite = FALSE`.
#'   * `bipartite` - logical echoing the input.
#'   * `n_edges` - integer count of input rows read.
#'   * `n_nonzero` - integer number of non-zero cells in `rr`.
#'
#' @seealso [moneca_bipartite()] for the bipartite clustering path that
#'   consumes `$rr`; [moneca_sbm()] and [moneca_flow()] for the square
#'   clustering paths that consume `$counts`.
#'
#' @examples
#' \dontrun{
#' # In-memory data.frame path
#' edges <- data.frame(
#'   row = c("a", "a", "b", "b", "c"),
#'   col = c("b", "c", "a", "c", "a"),
#'   weight = c(3, 1, 2, 4, 1)
#' )
#' res <- rr_from_duckdb(edges, weight_col = "weight")
#' res$rr
#'
#' # Parquet file path
#' # rr_from_duckdb("mobility.parquet", row_col = "worker",
#' #                col_col = "employer", weight_col = "n")
#'
#' # Pre-existing DuckDB connection
#' # con <- DBI::dbConnect(duckdb::duckdb(), "mobility.duckdb")
#' # rr_from_duckdb("mobility_edges", con = con,
#' #                row_col = "worker", col_col = "employer",
#' #                weight_col = "n")
#' # DBI::dbDisconnect(con, shutdown = TRUE)
#' }
#'
#' @export
rr_from_duckdb <- function(
  edges,
  row_col = "row",
  col_col = "col",
  weight_col = "weight",
  con = NULL,
  bipartite = FALSE,
  min_count = 0,
  small_cell_reduction = 0,
  zero_diagonal = TRUE,
  verbose = FALSE,
  ...
) {
  # 1. Dependency guard ----------------------------------------------------
  .require_duckdb()

  # 2. Argument validation -------------------------------------------------
  stopifnot(
    is.logical(bipartite),
    length(bipartite) == 1L,
    !is.na(bipartite),
    is.logical(zero_diagonal),
    length(zero_diagonal) == 1L,
    !is.na(zero_diagonal),
    is.logical(verbose),
    length(verbose) == 1L,
    !is.na(verbose),
    is.numeric(min_count),
    length(min_count) == 1L,
    min_count >= 0,
    is.numeric(small_cell_reduction),
    length(small_cell_reduction) == 1L,
    small_cell_reduction >= 0
  )

  if (
    is.character(edges) &&
      length(edges) == 1L &&
      !grepl("\\.parquet$|\\.csv(\\.gz)?$", edges, ignore.case = TRUE) &&
      is.null(con)
  ) {
    stop(
      "When `edges` is a bare table name, `con` must be a live DBI ",
      "connection to DuckDB.",
      call. = FALSE
    )
  }

  # 3. Connection management ----------------------------------------------
  own_con <- is.null(con)
  if (own_con) {
    con <- .open_transient_duckdb()
    on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  }

  # 4. Register / canonicalise the edge source ----------------------------
  view_name <- .register_edges(
    con = con,
    edges = edges,
    row_col = row_col,
    col_col = col_col,
    weight_col = weight_col
  )

  # 5. Count input edges (after canonicalisation, before filters) ---------
  n_edges <- as.integer(
    DBI::dbGetQuery(
      con,
      sprintf("SELECT COUNT(*) AS n FROM %s", .quote_ident(view_name))
    )$n
  )

  # 6. Small-cell reduction -----------------------------------------------
  if (small_cell_reduction > 0) {
    view_name <- .apply_scr(con, view_name, small_cell_reduction)
  }

  # 7. Main RR SQL --------------------------------------------------------
  sql_result <- .compute_rr_sql(
    con = con,
    view_name = view_name,
    zero_diagonal = zero_diagonal,
    min_count = min_count,
    bipartite = bipartite
  )

  if (verbose) {
    message(
      "[rr_from_duckdb] ",
      nrow(sql_result$edges),
      " non-zero RR cells from ",
      n_edges,
      " input edges."
    )
  }

  # 8. Sparse materialisation ---------------------------------------------
  out <- .result_to_sparse(sql_result, bipartite = bipartite)

  out$bipartite <- bipartite
  out$n_edges <- n_edges
  out$n_nonzero <- length(out$rr@x)
  class(out) <- "monecascale_rr"
  out
}
