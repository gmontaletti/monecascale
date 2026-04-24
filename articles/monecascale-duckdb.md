# Out-of-core RR computation with DuckDB

## 1. Introduction: the RR bottleneck

MONECA clusters mobility data on the relative-risk matrix `RR = O / E`,
where `E_ij = r_i * c_j / T` is the expected count under row / column
independence. The reference implementation
[`moneca::weight.matrix()`](https://gmontaletti.github.io/MONECA/reference/weight.matrix.html)
densifies the input on its way to RR (`analytical_functions.R:581-595`):
the observed count matrix `O` is materialised in full before the ratio
is formed. That densification sets a hard R-side ceiling on problem
size. A mobility graph on 10^5+ nodes may fit comfortably as a sparse
edge list, yet its dense matrix representation does not fit in R memory
at all.

[`rr_from_duckdb()`](https://gmontaletti.github.io/monecascale/reference/rr_from_duckdb.md)
lifts the RR construction step out of R and into DuckDB SQL. The
margins, the grand total, and the cell-wise ratio are computed as an
aggregation over the edge list; only the non-zero cells are ever
materialised, and only once, on their way into a sparse `dgCMatrix`. The
densification is removed at its root rather than worked around.

The output composes directly with monecascale’s sparse-aware backends.
[`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md)
consumes the RR matrix via its rectangular argument path;
[`moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)
and
[`moneca_flow()`](https://gmontaletti.github.io/monecascale/reference/moneca_flow.md)
consume the raw counts exposed on `$counts`. The end result is that
mobility datasets an order of magnitude larger than
[`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html)
can handle become tractable, without giving up the moneca-class contract
that downstream analysis and plotting APIs rely on.

## 2. The SQL recipe

The RR computation is a three-CTE aggregation over the edge list. One
CTE sums the counts per row; one sums per column; one takes the grand
total. A final join multiplies the observed count by
`total / (row_sum * col_sum)` to yield `rr`.

``` sql
WITH row_sums AS (
  SELECT row_id, SUM(w) AS s FROM edges GROUP BY row_id
),
col_sums AS (
  SELECT col_id, SUM(w) AS s FROM edges GROUP BY col_id
),
totals AS (
  SELECT SUM(w) AS total FROM edges
)
SELECT
  e.row_id,
  e.col_id,
  e.w AS count,
  e.w * t.total / (r.s * c.s) AS rr
FROM edges e
JOIN row_sums r ON r.row_id = e.row_id
JOIN col_sums c ON c.col_id = e.col_id
CROSS JOIN totals t;
```

No ordering is required, and no dense intermediate is constructed.
DuckDB streams the aggregations over the edge list, then streams the
join result back to R in long form. The R side receives a
`(row_id, col_id, count, rr)` table whose length equals the number of
non-zero cells in `O`, and feeds that directly into
[`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).
`$margins$total` reflects the full input total rather than a post-filter
total, because that is what the RR denominator uses.

## 3. Three input modes

[`rr_from_duckdb()`](https://gmontaletti.github.io/monecascale/reference/rr_from_duckdb.md)
accepts three shapes of input. The examples below run on a small square
fixture so they evaluate in a vignette build; the same code shape scales
to arbitrarily large edge lists on disk.

### 3.1 In-memory data frame

``` r
library(monecascale)

set.seed(2026)
edges <- data.frame(
  row    = sample(paste0("n", 1:30), 200, replace = TRUE),
  col    = sample(paste0("n", 1:30), 200, replace = TRUE),
  weight = rpois(200, 3) + 1L,
  stringsAsFactors = FALSE
)

res <- rr_from_duckdb(edges, weight_col = "weight", verbose = TRUE)
str(res, max.level = 1)
#> List of 7
#>  $ counts    :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>  $ rr        :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>  $ margins   :List of 3
#>  $ node_names:List of 2
#>  $ bipartite : logi FALSE
#>  $ n_edges   : int 200
#>  $ n_nonzero : int 174
#>  - attr(*, "class")= chr "monecascale_rr"
```

The return value is a `"monecascale_rr"` list: `$rr` is the sparse RR
matrix, `$counts` the sparse observed count matrix,
`$margins$row / $col / $total` the raw margins (with `$total` being the
full pre-diagonal-filter grand total used as the RR denominator),
`$node_names$row / $col` the dimnames, `$bipartite` the mode flag, and
`$n_edges / $n_nonzero` integer summaries.

### 3.2 On-disk CSV or Parquet

When the edge list is already on disk, pass a path instead of a data
frame. DuckDB reads CSV natively via `read_csv_auto()`; R never holds
the edges in memory.

``` r
path <- tempfile(fileext = ".csv")
utils::write.csv(edges, path, row.names = FALSE)

res_csv <- rr_from_duckdb(path, weight_col = "weight")
identical(dim(res_csv$rr), dim(res$rr))
#> [1] TRUE
```

Parquet input works identically: DuckDB dispatches on the file extension
to `read_parquet()`. The Parquet path is not shown here because it would
introduce an `arrow` dependency for a demonstration, but the calling
convention is unchanged.

### 3.3 Pre-connected DuckDB table

The third mode is for users who already have a DuckDB connection with
data registered. When `edges` is a character scalar that does not
resolve to a file on disk,
[`rr_from_duckdb()`](https://gmontaletti.github.io/monecascale/reference/rr_from_duckdb.md)
treats it as a table or view name on the supplied connection. This is
the pattern for streaming Parquet datasets or remote tables that a user
attaches once and queries many times.

``` r
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
DBI::dbWriteTable(con, "edges_tbl", edges)

res_con <- rr_from_duckdb("edges_tbl", con = con, weight_col = "weight")
identical(dim(res_con$rr), dim(res$rr))
#> [1] TRUE

DBI::dbDisconnect(con, shutdown = TRUE)
```

## 4. Composition with `moneca_bipartite()`

For rectangular input (`bipartite = TRUE`) the RR matrix is the natural
input to
[`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md).
A 30 x 20 fixture is small enough to fit in a vignette and large enough
to illustrate the composition.

``` r
set.seed(2026)
edges_bip <- data.frame(
  row    = paste0("r", sample(1:30, 300, replace = TRUE)),
  col    = paste0("c", sample(1:20, 300, replace = TRUE)),
  weight = rpois(300, 3) + 1L
)

res_bip <- rr_from_duckdb(edges_bip, weight_col = "weight", bipartite = TRUE)
dim(res_bip$rr)
#> [1] 30 20
res_bip$bipartite
#> [1] TRUE
```

``` r
fit <- monecascale::moneca_bipartite(as.matrix(res_bip$rr))
print(fit)
#> <moneca_bipartite> rows: 30 x cols: 20, backend: greed/DcLbm, levels: 27
#>  level Krow Kcol joint_ICL
#>      1   15   14 -3403.192
#>      2   14   14 -3403.607
#>      3   13   14 -3405.961
#>      4   13   13 -3406.968
#>      5   12   13 -3410.892
#>      6   11   13 -3415.326
#>      7   11   12 -3417.523
#>      8   10   12 -3420.591
#>      9   10   11 -3421.500
#>     10    9   11 -3425.928
#>     11    9   10 -3427.863
#>     12    9    9 -3431.576
#>     13    9    8 -3438.874
#>     14    8    8 -3445.425
#>     15    8    7 -3454.015
#>     16    7    7 -3460.223
#>     17    7    6 -3467.302
#>     18    6    6 -3471.411
#>     19    5    6 -3476.887
#>     20    5    5 -3484.651
#>     21    4    5 -3492.419
#>     22    4    4 -3494.997
#>     23    3    4 -3505.551
#>     24    3    3 -3497.945
#>     25    2    3 -3520.937
#>     26    2    2 -3516.221
#>     27    1    2 -3538.827
```

The dense coercion `as.matrix(res_bip$rr)` is used here because the
current
[`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md)
contract is tested on dense input. At true 10^5-row scale the dense cast
defeats the purpose of the out-of-core build; that path is reserved for
a sparse-aware bipartite fit planned as a follow-up. The point of this
section is the composition pattern, not the scale.

## 5. Composition with `moneca_sbm()` and `moneca_flow()` via `$counts`

[`moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)
and
[`moneca_flow()`](https://gmontaletti.github.io/monecascale/reference/moneca_flow.md)
consume raw mobility counts, not RR: the SBM backend builds its own
independence null internally, and Infomap needs the unnormalised flow
magnitudes to weight the random walk. For those entry points the useful
slot is `$counts`, exposed alongside `$rr` specifically so that one
DuckDB pass feeds both downstream paths.

``` r
res <- rr_from_duckdb(edges, weight_col = "weight")
dense_counts <- as.matrix(res$counts)
dim(dense_counts)
#> [1] 30 30
```

``` r
fit_flow <- monecascale::moneca_flow(
  dense_counts,
  segment.levels = 2L,
  seed = 2026
)
print(fit_flow)
#> 
#> ================================================================================
#>                         moneca MOBILITY ANALYSIS RESULTS                        
#> ================================================================================
#> 
#> OVERALL MOBILITY PATTERNS
#> -------------------------------------------------------------------------------
#> Overall Population Mobility Rate:                   100.0%
#> Average Mobility Concentration (all levels):         99.9%
#> 
#> HIERARCHICAL SEGMENTATION ANALYSIS
#> -------------------------------------------------------------------------------
#> 
#> Internal Mobility Within Segments (%):
#> Level 1 Level 2 Level 3 
#>       0     100     100 
#> 
#> Mobility Concentration in Significant Pathways by Level (%):
#> Level 1 Level 2 Level 3 
#>    99.6   100.0   100.0 
#> 
#> Network Structure by Level:
#>                                    Level 1      Level 2      Level 3 
#> -------------------------------------------------------------------------------
#> Active Segments/Classes:                30            1            1 
#> Significant Edges:                     171            0            0 
#> Network Density:                     0.197          NaN          NaN 
#> Isolated Segments:                       0            1            1 
#> 
#> DETAILED WEIGHTED DEGREE DISTRIBUTIONS (STRENGTH)
#> -------------------------------------------------------------------------------
#> 
#> Total Weighted Connections (Strength In + Out):
#>           Min    Q1 Median  Mean    Q3   Max
#> Level 1 43.72 53.66   57.7 59.09 63.38 88.63
#> Level 2  0.00  0.00    0.0  0.00  0.00  0.00
#> Level 3  0.00  0.00    0.0  0.00  0.00  0.00
#> 
#> Outward Mobility Strength (Weighted Out-Degree):
#>           Min    Q1 Median  Mean    Q3  Max
#> Level 1 22.79 26.68  29.29 29.54 31.58 42.1
#> Level 2  0.00  0.00   0.00  0.00  0.00  0.0
#> Level 3  0.00  0.00   0.00  0.00  0.00  0.0
#> 
#> Inward Mobility Strength (Weighted In-Degree):
#>          Min    Q1 Median  Mean    Q3   Max
#> Level 1 19.7 25.06  27.78 29.54 33.01 51.52
#> Level 2  0.0  0.00   0.00  0.00  0.00  0.00
#> Level 3  0.0  0.00   0.00  0.00  0.00  0.00
#> 
#> Edge Weight Distribution (Relative Risk Values):
#>          Min   Q1 Median Mean   Q3   Max
#> Level 1 1.04 2.61    4.1 5.18 6.26 23.15
#> Level 2   NA   NA     NA  NaN   NA    NA
#> Level 3   NA   NA     NA  NaN   NA    NA
#> 
#> ================================================================================
```

``` r
fit_sbm <- monecascale::moneca_sbm(dense_counts, backend = "greed")
print(fit_sbm)
#> 
#> ================================================================================
#>                         moneca MOBILITY ANALYSIS RESULTS                        
#> ================================================================================
#> 
#> OVERALL MOBILITY PATTERNS
#> -------------------------------------------------------------------------------
#> Overall Population Mobility Rate:                   100.0%
#> Average Mobility Concentration (all levels):         83.4%
#> 
#> HIERARCHICAL SEGMENTATION ANALYSIS
#> -------------------------------------------------------------------------------
#> 
#> Internal Mobility Within Segments (%):
#>  Level 1  Level 2  Level 3  Level 4  Level 5  Level 6  Level 7  Level 8 
#>      0.0      0.6      0.6      1.7      1.7      1.7      1.7      2.0 
#>  Level 9 Level 10 Level 11 Level 12 Level 13 Level 14 Level 15 Level 16 
#>      2.0      2.0      3.0      3.0      3.9      5.3      7.1      7.1 
#> Level 17 Level 18 Level 19 Level 20 Level 21 Level 22 
#>      8.9     10.0     12.7     17.7     28.6     75.3 
#> 
#> Mobility Concentration in Significant Pathways by Level (%):
#>  Level 1  Level 2  Level 3  Level 4  Level 5  Level 6  Level 7  Level 8 
#>     99.6     95.0     94.2     93.6     92.8     92.9     92.4     91.6 
#>  Level 9 Level 10 Level 11 Level 12 Level 13 Level 14 Level 15 Level 16 
#>     90.1     88.9     89.6     86.8     84.1     79.9     77.8     74.5 
#> Level 17 Level 18 Level 19 Level 20 Level 21 Level 22 
#>     71.9     70.6     68.6     64.5     60.0     75.3 
#> 
#> Network Structure by Level:
#>                                    Level 1      Level 2      Level 3      Level 4      Level 5      Level 6      Level 7      Level 8      Level 9     Level 10     Level 11     Level 12     Level 13     Level 14     Level 15     Level 16     Level 17     Level 18     Level 19     Level 20     Level 21     Level 22 
#> -------------------------------------------------------------------------------
#> Active Segments/Classes:                30           22           21           20           19           18           17           16           15           14           13           12           11           10            9            8            7            6            5            4            3            2 
#> Significant Edges:                     171          121          113          105           97           93           89           85           76           69           62           52           43           36           28           22           17           13            9            5            2            0 
#> Network Density:                     0.197        0.262        0.269        0.276        0.284        0.304        0.327        0.354        0.362        0.379        0.397        0.394        0.391        0.400        0.389        0.393        0.405        0.433        0.450        0.417        0.333        0.000 
#> Isolated Segments:                       0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            0            1            2 
#> 
#> DETAILED WEIGHTED DEGREE DISTRIBUTIONS (STRENGTH)
#> -------------------------------------------------------------------------------
#> 
#> Total Weighted Connections (Strength In + Out):
#>            Min    Q1 Median  Mean    Q3   Max
#> Level 1  43.72 53.66  57.70 59.09 63.38 88.63
#> Level 2  29.56 37.92  42.78 43.49 47.31 61.53
#> Level 3  28.15 37.12  41.41 42.08 45.29 61.53
#> Level 4  22.88 34.60  40.85 40.56 44.82 61.53
#> Level 5  21.89 34.52  38.54 38.56 42.96 59.06
#> Level 6  21.54 31.93  36.76 36.81 43.53 49.59
#> Level 7  20.07 30.28  31.66 32.98 36.49 43.41
#> Level 8  20.07 27.58  30.22 31.10 34.94 41.81
#> Level 9  20.07 25.32  27.83 28.39 31.28 38.06
#> Level 10 18.92 22.27  25.77 25.97 29.92 34.59
#> Level 11 18.92 21.42  22.13 23.63 26.35 30.95
#> Level 12 17.51 18.58  20.52 21.34 21.42 30.95
#> Level 13 15.93 16.39  17.77 17.96 18.71 22.18
#> Level 14 12.39 14.28  16.48 16.18 17.65 20.34
#> Level 15 11.35 11.69  13.57 13.92 14.94 19.18
#> Level 16 10.16 11.08  11.53 11.39 11.84 12.34
#> Level 17  6.92  8.41   9.61  9.51 10.47 12.30
#> Level 18  6.80  7.30   8.56  8.30  9.24  9.54
#> Level 19  4.68  5.71   7.64  6.88  7.95  8.41
#> Level 20  3.41  4.18   5.68  5.68  7.17  7.95
#> Level 21  0.00  1.70   3.41  2.27  3.41  3.41
#> Level 22  0.00  0.00   0.00  0.00  0.00  0.00
#> 
#> Outward Mobility Strength (Weighted Out-Degree):
#>            Min    Q1 Median  Mean    Q3   Max
#> Level 1  22.79 26.68  29.29 29.54 31.58 42.10
#> Level 2  13.18 19.07  21.25 21.75 24.40 31.69
#> Level 3  12.88 18.61  20.42 21.04 22.79 30.93
#> Level 4  11.22 17.79  19.97 20.28 21.64 30.93
#> Level 5  11.22 16.77  18.61 19.28 21.28 29.57
#> Level 6  11.22 15.02  17.40 18.40 21.13 28.61
#> Level 7   9.76 14.15  14.98 16.49 18.92 26.02
#> Level 8   9.76 12.96  14.63 15.55 15.80 24.69
#> Level 9   9.76 11.45  13.92 14.20 15.03 21.55
#> Level 10  9.58 10.64  12.77 12.98 14.04 20.67
#> Level 11  9.57  9.79  10.89 11.81 12.68 19.07
#> Level 12  8.08  9.36   9.69 10.67 10.85 19.07
#> Level 13  7.50  8.05   8.47  8.98  9.75 11.69
#> Level 14  5.96  7.35   8.18  8.09  8.52 11.13
#> Level 15  5.64  5.88   6.60  6.96  6.79 11.13
#> Level 16  4.71  5.49   5.70  5.70  6.02  6.60
#> Level 17  3.60  4.16   4.84  4.76  5.40  5.72
#> Level 18  3.43  3.69   4.18  4.15  4.53  4.95
#> Level 19  1.71  3.57   3.93  3.44  3.96  4.01
#> Level 20  1.61  1.75   2.87  2.84  3.95  4.01
#> Level 21  0.00  0.80   1.61  1.14  1.70  1.80
#> Level 22  0.00  0.00   0.00  0.00  0.00  0.00
#> 
#> Inward Mobility Strength (Weighted In-Degree):
#>            Min    Q1 Median  Mean    Q3   Max
#> Level 1  19.70 25.06  27.78 29.54 33.01 51.52
#> Level 2  15.69 17.55  20.16 21.75 24.08 31.55
#> Level 3  14.39 16.99  19.19 21.04 24.45 31.55
#> Level 4  11.66 16.34  18.79 20.28 23.16 30.60
#> Level 5  10.67 16.14  18.65 19.28 21.51 29.49
#> Level 6  10.32 15.29  16.84 18.40 20.96 28.84
#> Level 7  10.32 14.51  16.53 16.49 18.65 21.72
#> Level 8  10.32 14.39  15.73 15.55 17.13 20.74
#> Level 9  10.32 13.34  14.42 14.20 15.52 17.15
#> Level 10  9.00 11.66  13.05 12.98 13.86 17.15
#> Level 11  9.00 10.52  11.88 11.81 12.33 16.84
#> Level 12  8.31  9.38  10.29 10.67 11.48 15.53
#> Level 13  7.35  8.41   9.13  8.98  9.48 10.48
#> Level 14  6.10  7.00   8.57  8.09  9.18  9.44
#> Level 15  4.96  6.05   6.88  6.96  8.05  9.13
#> Level 16  4.27  4.84   5.61  5.70  6.63  7.29
#> Level 17  3.32  3.95   4.77  4.76  5.37  6.58
#> Level 18  2.39  3.61   4.48  4.15  4.85  5.26
#> Level 19  1.75  2.97   3.93  3.44  4.07  4.48
#> Level 20  1.80  2.43   2.81  2.84  3.22  3.93
#> Level 21  0.00  0.80   1.61  1.14  1.70  1.80
#> Level 22  0.00  0.00   0.00  0.00  0.00  0.00
#> 
#> Edge Weight Distribution (Relative Risk Values):
#>           Min   Q1 Median Mean   Q3   Max
#> Level 1  1.04 2.61   4.10 5.18 6.26 23.15
#> Level 2  1.07 2.08   3.11 3.95 4.47 15.90
#> Level 3  1.07 2.08   2.98 3.91 4.53 15.90
#> Level 4  1.07 1.99   2.88 3.86 4.47 15.90
#> Level 5  1.08 1.80   2.88 3.78 4.40 15.90
#> Level 6  1.01 1.71   2.67 3.56 4.20 14.93
#> Level 7  1.01 1.51   2.50 3.15 4.04 13.66
#> Level 8  1.01 1.49   2.41 2.93 3.71 13.66
#> Level 9  1.01 1.51   2.29 2.80 3.56 13.46
#> Level 10 1.01 1.48   2.14 2.63 3.26 13.46
#> Level 11 1.01 1.43   2.11 2.48 2.79 13.46
#> Level 12 1.01 1.47   1.98 2.46 2.72 13.46
#> Level 13 1.01 1.46   2.07 2.30 2.58  7.06
#> Level 14 1.04 1.49   2.03 2.25 2.45  7.06
#> Level 15 1.09 1.52   1.96 2.24 2.46  7.06
#> Level 16 1.03 1.44   1.84 2.07 2.50  3.93
#> Level 17 1.03 1.59   1.72 1.96 2.35  3.93
#> Level 18 1.01 1.50   1.71 1.92 2.39  3.93
#> Level 19 1.03 1.32   1.71 1.91 1.83  3.93
#> Level 20 1.03 1.61   1.80 2.27 2.98  3.93
#> Level 21 1.61 1.66   1.70 1.70 1.75  1.80
#> Level 22   NA   NA     NA  NaN   NA    NA
#> 
#> ================================================================================
```

`$rr` is the input for paths that use the relative-risk matrix directly,
such as
[`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md)
and the planned `moneca_localclique()`. `$counts` is the input for paths
that expect raw mobility. A single DuckDB aggregation yields both.

## 6. Limitations

- **No Arrow dataset input yet.** Partitioned Parquet datasets must be
  materialised through DuckDB’s own `read_parquet()` rather than
  registered as Arrow datasets. A direct Arrow path is planned as a
  follow-up.
- **No streaming RR persistence.** The sparse RR is returned to R in
  memory. Users whose RR does not fit in R memory need to persist the
  long-format cell table through DuckDB directly; a streaming writer is
  not exposed.
- **Downstream ceiling from moneca metadata.**
  [`moneca::moneca_segments()`](https://gmontaletti.github.io/MONECA/reference/moneca_segments.html),
  called internally by the SBM and flow adapters, still densifies
  `mat.list[[1]]` when it builds per-segment summaries. The full
  pipeline therefore retains an R-side ceiling around 10^5 nodes until
  upstream moneca adds a sparse metadata path.
  [`rr_from_duckdb()`](https://gmontaletti.github.io/monecascale/reference/rr_from_duckdb.md)
  removes the RR-build bottleneck; it does not remove the metadata-build
  one.
- **Optional dependency.** `duckdb` is a `Suggests`. Install with
  `install.packages("duckdb")`.
  [`rr_from_duckdb()`](https://gmontaletti.github.io/monecascale/reference/rr_from_duckdb.md)
  errors informatively when the package is not available.

## 7. References

- DuckDB: <https://duckdb.org/>.
- Raasveldt, M., & Muhleisen, H. (2019). DuckDB: an embeddable
  analytical database. *Proceedings of the 2019 International Conference
  on Management of Data*, 1981-1984.
- The monecascale scaling roadmap places
  [`rr_from_duckdb()`](https://gmontaletti.github.io/monecascale/reference/rr_from_duckdb.md)
  under direction D7 (out-of-core / streaming backend), shipped in
  version 0.5.0 alongside the SBM, bipartite, flow, and auto-level
  backends documented in the other vignettes.
