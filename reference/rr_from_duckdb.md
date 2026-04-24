# Build a Sparse Relative-Risk Matrix Out-of-Core via DuckDB

Aggregates an edge list into a sparse relative-risk (\`RR = O / E\`)
matrix plus the companion raw-count matrix, using DuckDB as the
aggregation engine. The R side never allocates a dense \`n_rows x
n_cols\` matrix, so the function scales to person x employer or worker x
occupation mobility tables that would OOM under
\`moneca::weight.matrix()\`.

## Usage

``` r
rr_from_duckdb(
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
)
```

## Arguments

- edges:

  One of: (a) an in-memory \`data.frame\` / \`data.table\` with columns
  named by \`row_col\`, \`col_col\`, \`weight_col\`; (b) a single
  character file path ending in \`.parquet\`, \`.csv\`, or \`.csv.gz\`;
  or (c) a character scalar naming a pre-existing DuckDB table / view on
  the supplied \`con\`.

- row_col, col_col, weight_col:

  Column names inside the source. Defaults are \`"row"\`, \`"col"\`,
  \`"weight"\`.

- con:

  Optional \`DBIConnection\` to DuckDB. Required when \`edges\` is a
  bare table name. When \`NULL\`, a transient in-memory DuckDB
  connection is opened for the call and closed on exit.

- bipartite:

  Logical. If \`FALSE\` (default), row and column identifiers live in
  the same namespace and are unioned into a shared factor, yielding a
  square output with identical row and column names. If \`TRUE\`,
  identifiers are disjoint and the output is rectangular with separate
  row / column names.

- min_count:

  Numeric. Drop edges with \`w \< min_count\` before the RR computation.
  Defaults to \`0\` (no filter).

- small_cell_reduction:

  Numeric. Drop rows and columns whose marginal sum is below this
  threshold before RR computation. Mirrors the \`small.cell.reduction\`
  argument of \`moneca::moneca_fast()\`. Defaults to \`0\` (no filter).

- zero_diagonal:

  Logical. If \`TRUE\` (default) and \`bipartite = FALSE\`, self-loops
  (\`row_id == col_id\`) are dropped before the RR computation. Ignored
  when \`bipartite = TRUE\`.

- verbose:

  Logical. If \`TRUE\`, message the number of non-zero RR cells
  returned.

- ...:

  Reserved for forward compatibility.

## Value

A list of class \`"monecascale_rr"\` with elements \* \`rr\` - sparse
\`dgCMatrix\` of observed / expected; one entry per non-zero input cell.
\* \`counts\` - sparse \`dgCMatrix\` of raw observed counts; same
sparsity pattern as \`rr\`. \* \`margins\` - list with numeric vectors
\`row\` and \`col\`, plus the scalar \`total\`. \* \`node_names\` - list
with character vectors \`row\` and \`col\`. Identical when \`bipartite =
FALSE\`. \* \`bipartite\` - logical echoing the input. \* \`n_edges\` -
integer count of input rows read. \* \`n_nonzero\` - integer number of
non-zero cells in \`rr\`.

## Details

The RR aggregation is pushed into a single SQL query with the template

“\`sql WITH row_m AS (SELECT row_id, SUM(w) AS s FROM edges GROUP BY
row_id), col_m AS (SELECT col_id, SUM(w) AS s FROM edges GROUP BY
col_id), tot AS (SELECT SUM(w) AS total FROM edges) SELECT e.row_id,
e.col_id, e.w AS count, e.w \* tot.total / (r.s \* c.s) AS rr FROM edges
e JOIN row_m r USING (row_id) JOIN col_m c USING (col_id) CROSS JOIN tot
WHERE e.w \> 0; “\`

The join-and-divide pattern only emits rows whose edge weight is
non-zero, so sparsity is preserved end to end: the result set has one
row per non-zero input edge, and the R-side \`Matrix::sparseMatrix()\`
call materialises it directly as a \`dgCMatrix\`. When \`bipartite =
FALSE\`, row and column IDs are unioned into a single factor so that
both output matrices are square with identical dimnames. When
\`bipartite = TRUE\`, row and column namespaces stay disjoint and the
output is rectangular.

Preprocessing is applied in order: optional small-cell reduction (drops
rows/cols whose marginal sum falls below a threshold), optional
\`min_count\` filter (drops individual edges), optional zero-diagonal
filter (only active when \`bipartite = FALSE\`). All three filters are
applied in SQL.

## See also

\[moneca_bipartite()\] for the bipartite clustering path that consumes
\`\$rr\`; \[moneca_sbm()\] and \[moneca_flow()\] for the square
clustering paths that consume \`\$counts\`.

## Examples

``` r
if (FALSE) { # \dontrun{
# In-memory data.frame path
edges <- data.frame(
  row = c("a", "a", "b", "b", "c"),
  col = c("b", "c", "a", "c", "a"),
  weight = c(3, 1, 2, 4, 1)
)
res <- rr_from_duckdb(edges, weight_col = "weight")
res$rr

# Parquet file path
# rr_from_duckdb("mobility.parquet", row_col = "worker",
#                col_col = "employer", weight_col = "n")

# Pre-existing DuckDB connection
# con <- DBI::dbConnect(duckdb::duckdb(), "mobility.duckdb")
# rr_from_duckdb("mobility_edges", con = con,
#                row_col = "worker", col_col = "employer",
#                weight_col = "n")
# DBI::dbDisconnect(con, shutdown = TRUE)
} # }
```
