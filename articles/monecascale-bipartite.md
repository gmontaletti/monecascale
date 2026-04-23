# Bipartite Segmentation with moneca_bipartite()

## 1. Introduction

[`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md)
extends the MONECA analysis stack to rectangular (two-mode) mobility
data, where rows and columns are distinct entity types. Typical use
cases are person x employer, worker x occupation, and student x school
tables. The square mobility matrix that
[`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html)
and
[`monecascale::moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)
operate on is already a one-mode aggregation; bipartite input preserves
individual row and column identity and lets the block model infer groups
on both sides jointly.

This vignette is the rectangular counterpart of
[`vignette("monecascale-sbm", package = "monecascale")`](https://gmontaletti.github.io/monecascale/articles/monecascale-sbm.md):
the square case motivates why relative risk is the DC-SBM residual, this
one motivates why the same argument carries over to the Latent Block
Model (LBM) and how the joint row / column clustering is exposed to
moneca’s downstream `segment.*()` and `plot_moneca_*()` APIs.

## 2. Theoretical bridge: rectangular RR as DcLbm residual

MONECA’s weight is the relative-risk deviation of observed mobility from
the independence null. For a rectangular core `O` with row margins `r`,
column margins `c`, and grand total `N`:

$$RR_{ij} = \frac{O_{ij}}{E_{ij}},\qquad E_{ij} = \frac{r_{i} \cdot c_{j}}{N}.$$

The degree-corrected Latent Block Model (DcLbm, Come et al. 2021) models
the expected edge count between row `i` in row-block `g` and column `j`
in column-block `h` as `theta_i^row * theta_j^col * omega_{g,h}`, with
`theta` degree parameters and `omega` the row-block by column-block
interaction matrix. Setting `omega` to the all-ones matrix recovers
exactly `(r_i c_j) / N`, i.e. MONECA’s independence null. DcLbm
inference on the rectangular mobility is therefore the principled
scalable drop-in for clique enumeration on the bipartite case, in the
same way that DcSbm is for the unipartite case covered in
`monecascale-sbm`.

The current backend is
[`greed::DcLbm`](https://comeetie.github.io/greed/reference/DcLbm.html),
an R-native hierarchical implementation built on Integrated
Classification Likelihood with merge-split local search. A `"graphtool"`
backend via `graph_tool.inference.minimize_blockmodel_dl(pclabel = ...)`
is reserved for scale beyond 10^5 per side.

## 3. Quick start

``` r
library(monecascale)

# 1. synthetic rectangular fixture -----
set.seed(2026)
n_row <- 60
n_col <- 40
row_block <- rep(1:4, length.out = n_row)
col_block <- rep(1:3, length.out = n_col)
# block-constant Poisson intensities: diagonal-ish structure
lambda <- outer(row_block, col_block, function(g, h) {
  ifelse((g - 1) %% 3 + 1 == h, 6, 1)
})
mx <- matrix(stats::rpois(n_row * n_col, lambda = lambda), n_row, n_col)
rownames(mx) <- paste0("R", seq_len(n_row))
colnames(mx) <- paste0("C", seq_len(n_col))

# 2. fit -----
out <- monecascale::moneca_bipartite(mx, seed = 2026, verbose = FALSE)

# 3. inspect -----
out
#> <moneca_bipartite> rows: 60 x cols: 40, backend: greed/DcLbm, levels: 4
#>  level Krow Kcol joint_ICL
#>      1    3    3 -62529.90
#>      2    2    3 -63607.83
#>      3    2    2 -63588.51
#>      4    1    2 -64976.75
summary(out)
#> <moneca_bipartite> rows: 60 x cols: 40, backend: greed/DcLbm, levels: 4
#>   backend: greed   model: DcLbm   n_init: 1   rr_scale: 3
#> Hierarchy levels:
#>  level Krow Kcol joint_ICL joint_MDL
#>      1    3    3 -62529.90  62529.90
#>      2    2    3 -63607.83  63607.83
#>      3    2    2 -63588.51  63588.51
#>      4    1    2 -64976.75  64976.75
out$bipartite_diagnostics$n_blocks_row_per_level
#> [1] 3 2 2 1
out$bipartite_diagnostics$n_blocks_col_per_level
#> [1] 3 3 2 2
```

The
[`print.moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/print.moneca_bipartite.md)
method shows the per-level table of row block count, column block count,
and joint ICL. [`summary()`](https://rdrr.io/r/base/summary.html) adds
the backend metadata block (`n_init`, `rr_scale`, margin status).

## 4. Output anatomy

[`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md)
returns an object of class `"moneca_bipartite"` (disjoint from
`"moneca"`) with three slots: `$rows`, `$cols`, and
`$bipartite_diagnostics`. The first two are full moneca-class objects
built from one-mode projections of the rectangular RR. The third retains
the joint bipartite information that a single square object cannot
represent.

``` r
class(out)
#> [1] "moneca_bipartite"
class(out$rows)
#> [1] "moneca"
class(out$cols)
#> [1] "moneca"
```

The row view is a square, symmetric, zero-diagonal matrix on `n_row`
entities; the column view is square on `n_col`. Margins are appended as
the last row and column, matching moneca’s convention.

``` r
dim(out$rows$mat.list[[1]])
#> [1] 61 61
dim(out$cols$mat.list[[1]])
#> [1] 41 41
isSymmetric(as.matrix(out$rows$mat.list[[1]]))
#> [1] TRUE
```

The diagnostics slot exposes the joint structure:

``` r
str(out$bipartite_diagnostics, max.level = 1)
#> List of 18
#>  $ backend                           : chr "greed"
#>  $ model                             : chr "DcLbm"
#>  $ n_row                             : int 60
#>  $ n_col                             : int 40
#>  $ n_blocks_row_per_level            : int [1:4] 3 2 2 1
#>  $ n_blocks_col_per_level            : int [1:4] 3 3 2 2
#>  $ joint_icl_per_level               : num [1:4] -62530 -63608 -63589 -64977
#>  $ joint_mdl_per_level               : num [1:4] 62530 63608 63589 64977
#>  $ block_interaction_matrix_per_level:List of 4
#>  $ rr_rect                           :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
#>  $ row_margin                        : Named num [1:60] 106 110 108 96 113 95 93 121 111 110 ...
#>   ..- attr(*, "names")= chr [1:60] "R1" "R2" "R3" "R4" ...
#>  $ col_margin                        : Named num [1:40] 184 137 135 187 147 138 200 124 138 206 ...
#>   ..- attr(*, "names")= chr [1:40] "C1" "C2" "C3" "C4" ...
#>  $ grand_total                       : int 6443
#>  $ rr_scale                          : int 3
#>  $ n_init                            : int 1
#>  $ margins_added                     : logi TRUE
#>  $ row_names                         : chr [1:60] "R1" "R2" "R3" "R4" ...
#>  $ col_names                         : chr [1:40] "C1" "C2" "C3" "C4" ...
out$bipartite_diagnostics$n_blocks_row_per_level
#> [1] 3 2 2 1
out$bipartite_diagnostics$n_blocks_col_per_level
#> [1] 3 3 2 2
round(out$bipartite_diagnostics$joint_icl_per_level, 2)
#> [1] -62529.90 -63607.83 -63588.51 -64976.75
```

`block_interaction_matrix_per_level` is a list of `Krow x Kcol`
matrices, one per hierarchy level, cross-tabulating the rectangular core
by row and column memberships. Level 1 is finest.

``` r
bim <- out$bipartite_diagnostics$block_interaction_matrix_per_level
length(bim)
#> [1] 4
dim(bim[[1]])
#> [1] 3 3
round(bim[[1]], 1)
#>      C1   C2   C3
#> R1 1212  186  185
#> R2  221 1122  238
#> R3  371  395 2513
```

## 5. Downstream compatibility

Because `$rows` and `$cols` are full moneca-class objects, the entire
moneca consumer API runs on them unchanged. The two views are queried
independently: the row view answers questions about how row entities
cluster given their column-mediated co-movement, and symmetrically for
the column view.

``` r
head(moneca::segment.membership(out$rows))
#>   name membership
#> 1   R1        5.1
#> 2   R2        5.1
#> 3   R3        5.1
#> 4   R4        5.1
#> 5   R5        5.1
#> 6   R6        5.1
head(moneca::segment.membership(out$cols))
#>   name membership
#> 1   C1        5.2
#> 2   C2        5.1
#> 3   C3        5.1
#> 4   C4        5.2
#> 5   C5        5.1
#> 6   C6        5.1
```

``` r
moneca::segment.quality(out$rows)
#>     Membership 1: Segment 1: within.mobility 1: share.of.mobility 1: Density
#> R3         5.1          3                 NA                   NA        NaN
#> R11        5.1         11                 NA                   NA        NaN
#> R23        5.1         23                 NA                   NA        NaN
#> R27        5.1         27                 NA                   NA        NaN
#> R39        5.1         39                 NA                   NA        NaN
#> R7         5.1          7                 NA                   NA        NaN
#> R15        5.1         15                 NA                   NA        NaN
#> R19        5.1         19                 NA                   NA        NaN
#> R31        5.1         31                 NA                   NA        NaN
#> R35        5.1         35                 NA                   NA        NaN
#> R43        5.1         43                 NA                   NA        NaN
#> R47        5.1         47                 NA                   NA        NaN
#> R51        5.1         51                 NA                   NA        NaN
#> R55        5.1         55                 NA                   NA        NaN
#> R59        5.1         59                 NA                   NA        NaN
#> R6         5.1          6                 NA                   NA        NaN
#> R10        5.1         10                 NA                   NA        NaN
#> R14        5.1         14                 NA                   NA        NaN
#> R18        5.1         18                 NA                   NA        NaN
#> R38        5.1         38                 NA                   NA        NaN
#> R42        5.1         42                 NA                   NA        NaN
#> R54        5.1         54                 NA                   NA        NaN
#> R2         5.1          2                 NA                   NA        NaN
#> R22        5.1         22                 NA                   NA        NaN
#> R26        5.1         26                 NA                   NA        NaN
#> R34        5.1         34                 NA                   NA        NaN
#> R50        5.1         50                 NA                   NA        NaN
#> R58        5.1         58                 NA                   NA        NaN
#> R30        5.1         30                 NA                   NA        NaN
#> R46        5.1         46                 NA                   NA        NaN
#> R41        5.1         41                 NA                   NA        NaN
#> R1         5.1          1                 NA                   NA        NaN
#> R4         5.1          4                 NA                   NA        NaN
#> R13        5.1         13                 NA                   NA        NaN
#> R16        5.1         16                 NA                   NA        NaN
#> R28        5.1         28                 NA                   NA        NaN
#> R32        5.1         32                 NA                   NA        NaN
#> R33        5.1         33                 NA                   NA        NaN
#> R45        5.1         45                 NA                   NA        NaN
#> R49        5.1         49                 NA                   NA        NaN
#> R52        5.1         52                 NA                   NA        NaN
#> R53        5.1         53                 NA                   NA        NaN
#> R56        5.1         56                 NA                   NA        NaN
#> R57        5.1         57                 NA                   NA        NaN
#> R60        5.1         60                 NA                   NA        NaN
#> R5         5.1          5                 NA                   NA        NaN
#> R8         5.1          8                 NA                   NA        NaN
#> R9         5.1          9                 NA                   NA        NaN
#> R12        5.1         12                 NA                   NA        NaN
#> R17        5.1         17                 NA                   NA        NaN
#> R20        5.1         20                 NA                   NA        NaN
#> R21        5.1         21                 NA                   NA        NaN
#> R24        5.1         24                 NA                   NA        NaN
#> R25        5.1         25                 NA                   NA        NaN
#> R29        5.1         29                 NA                   NA        NaN
#> R36        5.1         36                 NA                   NA        NaN
#> R37        5.1         37                 NA                   NA        NaN
#> R40        5.1         40                 NA                   NA        NaN
#> R48        5.1         48                 NA                   NA        NaN
#> R44        5.1         44                 NA                   NA        NaN
#>     1: Nodes 1: Max.path 1: share.of.total 2: Segment 2: within.mobility
#> R3         1           0             0.022          2                 NA
#> R11        1           0             0.022          2                 NA
#> R23        1           0             0.022          2                 NA
#> R27        1           0             0.022          2                 NA
#> R39        1           0             0.022          2                 NA
#> R7         1           0             0.021          2                 NA
#> R15        1           0             0.021          2                 NA
#> R19        1           0             0.021          2                 NA
#> R31        1           0             0.021          2                 NA
#> R35        1           0             0.021          2                 NA
#> R43        1           0             0.021          2                 NA
#> R47        1           0             0.021          2                 NA
#> R51        1           0             0.021          2                 NA
#> R55        1           0             0.021          2                 NA
#> R59        1           0             0.021          2                 NA
#> R6         1           0             0.021          1                 NA
#> R10        1           0             0.021          1                 NA
#> R14        1           0             0.021          1                 NA
#> R18        1           0             0.021          1                 NA
#> R38        1           0             0.021          1                 NA
#> R42        1           0             0.021          1                 NA
#> R54        1           0             0.021          1                 NA
#> R2         1           0             0.020          1                 NA
#> R22        1           0             0.020          1                 NA
#> R26        1           0             0.020          1                 NA
#> R34        1           0             0.020          1                 NA
#> R50        1           0             0.020          1                 NA
#> R58        1           0             0.020          1                 NA
#> R30        1           0             0.019          1                 NA
#> R46        1           0             0.019          1                 NA
#> R41        1           0             0.014          3                 NA
#> R1         1           0             0.013          3                 NA
#> R4         1           0             0.013          3                 NA
#> R13        1           0             0.013          3                 NA
#> R16        1           0             0.013          3                 NA
#> R28        1           0             0.013          3                 NA
#> R32        1           0             0.013          3                 NA
#> R33        1           0             0.013          3                 NA
#> R45        1           0             0.013          3                 NA
#> R49        1           0             0.013          3                 NA
#> R52        1           0             0.013          3                 NA
#> R53        1           0             0.013          3                 NA
#> R56        1           0             0.013          3                 NA
#> R57        1           0             0.013          3                 NA
#> R60        1           0             0.013          3                 NA
#> R5         1           0             0.012          3                 NA
#> R8         1           0             0.012          3                 NA
#> R9         1           0             0.012          3                 NA
#> R12        1           0             0.012          3                 NA
#> R17        1           0             0.012          3                 NA
#> R20        1           0             0.012          3                 NA
#> R21        1           0             0.012          3                 NA
#> R24        1           0             0.012          3                 NA
#> R25        1           0             0.012          3                 NA
#> R29        1           0             0.012          3                 NA
#> R36        1           0             0.012          3                 NA
#> R37        1           0             0.012          3                 NA
#> R40        1           0             0.012          3                 NA
#> R48        1           0             0.012          3                 NA
#> R44        1           0             0.011          3                 NA
#>     2: share.of.mobility 2: Density 2: Nodes 2: Max.path 2: share.of.total
#> R3                    NA          1       15           1             0.319
#> R11                   NA          1       15           1             0.319
#> R23                   NA          1       15           1             0.319
#> R27                   NA          1       15           1             0.319
#> R39                   NA          1       15           1             0.319
#> R7                    NA          1       15           1             0.319
#> R15                   NA          1       15           1             0.319
#> R19                   NA          1       15           1             0.319
#> R31                   NA          1       15           1             0.319
#> R35                   NA          1       15           1             0.319
#> R43                   NA          1       15           1             0.319
#> R47                   NA          1       15           1             0.319
#> R51                   NA          1       15           1             0.319
#> R55                   NA          1       15           1             0.319
#> R59                   NA          1       15           1             0.319
#> R6                    NA          1       15           1             0.304
#> R10                   NA          1       15           1             0.304
#> R14                   NA          1       15           1             0.304
#> R18                   NA          1       15           1             0.304
#> R38                   NA          1       15           1             0.304
#> R42                   NA          1       15           1             0.304
#> R54                   NA          1       15           1             0.304
#> R2                    NA          1       15           1             0.304
#> R22                   NA          1       15           1             0.304
#> R26                   NA          1       15           1             0.304
#> R34                   NA          1       15           1             0.304
#> R50                   NA          1       15           1             0.304
#> R58                   NA          1       15           1             0.304
#> R30                   NA          1       15           1             0.304
#> R46                   NA          1       15           1             0.304
#> R41                   NA          1       30           1             0.376
#> R1                    NA          1       30           1             0.376
#> R4                    NA          1       30           1             0.376
#> R13                   NA          1       30           1             0.376
#> R16                   NA          1       30           1             0.376
#> R28                   NA          1       30           1             0.376
#> R32                   NA          1       30           1             0.376
#> R33                   NA          1       30           1             0.376
#> R45                   NA          1       30           1             0.376
#> R49                   NA          1       30           1             0.376
#> R52                   NA          1       30           1             0.376
#> R53                   NA          1       30           1             0.376
#> R56                   NA          1       30           1             0.376
#> R57                   NA          1       30           1             0.376
#> R60                   NA          1       30           1             0.376
#> R5                    NA          1       30           1             0.376
#> R8                    NA          1       30           1             0.376
#> R9                    NA          1       30           1             0.376
#> R12                   NA          1       30           1             0.376
#> R17                   NA          1       30           1             0.376
#> R20                   NA          1       30           1             0.376
#> R21                   NA          1       30           1             0.376
#> R24                   NA          1       30           1             0.376
#> R25                   NA          1       30           1             0.376
#> R29                   NA          1       30           1             0.376
#> R36                   NA          1       30           1             0.376
#> R37                   NA          1       30           1             0.376
#> R40                   NA          1       30           1             0.376
#> R48                   NA          1       30           1             0.376
#> R44                   NA          1       30           1             0.376
#>     3: Segment 3: within.mobility 3: share.of.mobility 3: Density 3: Nodes
#> R3           1                 NA                   NA  0.4827586       30
#> R11          1                 NA                   NA  0.4827586       30
#> R23          1                 NA                   NA  0.4827586       30
#> R27          1                 NA                   NA  0.4827586       30
#> R39          1                 NA                   NA  0.4827586       30
#> R7           1                 NA                   NA  0.4827586       30
#> R15          1                 NA                   NA  0.4827586       30
#> R19          1                 NA                   NA  0.4827586       30
#> R31          1                 NA                   NA  0.4827586       30
#> R35          1                 NA                   NA  0.4827586       30
#> R43          1                 NA                   NA  0.4827586       30
#> R47          1                 NA                   NA  0.4827586       30
#> R51          1                 NA                   NA  0.4827586       30
#> R55          1                 NA                   NA  0.4827586       30
#> R59          1                 NA                   NA  0.4827586       30
#> R6           1                 NA                   NA  0.4827586       30
#> R10          1                 NA                   NA  0.4827586       30
#> R14          1                 NA                   NA  0.4827586       30
#> R18          1                 NA                   NA  0.4827586       30
#> R38          1                 NA                   NA  0.4827586       30
#> R42          1                 NA                   NA  0.4827586       30
#> R54          1                 NA                   NA  0.4827586       30
#> R2           1                 NA                   NA  0.4827586       30
#> R22          1                 NA                   NA  0.4827586       30
#> R26          1                 NA                   NA  0.4827586       30
#> R34          1                 NA                   NA  0.4827586       30
#> R50          1                 NA                   NA  0.4827586       30
#> R58          1                 NA                   NA  0.4827586       30
#> R30          1                 NA                   NA  0.4827586       30
#> R46          1                 NA                   NA  0.4827586       30
#> R41          2                 NA                   NA  1.0000000       30
#> R1           2                 NA                   NA  1.0000000       30
#> R4           2                 NA                   NA  1.0000000       30
#> R13          2                 NA                   NA  1.0000000       30
#> R16          2                 NA                   NA  1.0000000       30
#> R28          2                 NA                   NA  1.0000000       30
#> R32          2                 NA                   NA  1.0000000       30
#> R33          2                 NA                   NA  1.0000000       30
#> R45          2                 NA                   NA  1.0000000       30
#> R49          2                 NA                   NA  1.0000000       30
#> R52          2                 NA                   NA  1.0000000       30
#> R53          2                 NA                   NA  1.0000000       30
#> R56          2                 NA                   NA  1.0000000       30
#> R57          2                 NA                   NA  1.0000000       30
#> R60          2                 NA                   NA  1.0000000       30
#> R5           2                 NA                   NA  1.0000000       30
#> R8           2                 NA                   NA  1.0000000       30
#> R9           2                 NA                   NA  1.0000000       30
#> R12          2                 NA                   NA  1.0000000       30
#> R17          2                 NA                   NA  1.0000000       30
#> R20          2                 NA                   NA  1.0000000       30
#> R21          2                 NA                   NA  1.0000000       30
#> R24          2                 NA                   NA  1.0000000       30
#> R25          2                 NA                   NA  1.0000000       30
#> R29          2                 NA                   NA  1.0000000       30
#> R36          2                 NA                   NA  1.0000000       30
#> R37          2                 NA                   NA  1.0000000       30
#> R40          2                 NA                   NA  1.0000000       30
#> R48          2                 NA                   NA  1.0000000       30
#> R44          2                 NA                   NA  1.0000000       30
#>     3: Max.path 3: share.of.total 4: Segment 4: within.mobility
#> R3            1             0.624          1                 NA
#> R11           1             0.624          1                 NA
#> R23           1             0.624          1                 NA
#> R27           1             0.624          1                 NA
#> R39           1             0.624          1                 NA
#> R7            1             0.624          1                 NA
#> R15           1             0.624          1                 NA
#> R19           1             0.624          1                 NA
#> R31           1             0.624          1                 NA
#> R35           1             0.624          1                 NA
#> R43           1             0.624          1                 NA
#> R47           1             0.624          1                 NA
#> R51           1             0.624          1                 NA
#> R55           1             0.624          1                 NA
#> R59           1             0.624          1                 NA
#> R6            1             0.624          1                 NA
#> R10           1             0.624          1                 NA
#> R14           1             0.624          1                 NA
#> R18           1             0.624          1                 NA
#> R38           1             0.624          1                 NA
#> R42           1             0.624          1                 NA
#> R54           1             0.624          1                 NA
#> R2            1             0.624          1                 NA
#> R22           1             0.624          1                 NA
#> R26           1             0.624          1                 NA
#> R34           1             0.624          1                 NA
#> R50           1             0.624          1                 NA
#> R58           1             0.624          1                 NA
#> R30           1             0.624          1                 NA
#> R46           1             0.624          1                 NA
#> R41           1             0.376          2                 NA
#> R1            1             0.376          2                 NA
#> R4            1             0.376          2                 NA
#> R13           1             0.376          2                 NA
#> R16           1             0.376          2                 NA
#> R28           1             0.376          2                 NA
#> R32           1             0.376          2                 NA
#> R33           1             0.376          2                 NA
#> R45           1             0.376          2                 NA
#> R49           1             0.376          2                 NA
#> R52           1             0.376          2                 NA
#> R53           1             0.376          2                 NA
#> R56           1             0.376          2                 NA
#> R57           1             0.376          2                 NA
#> R60           1             0.376          2                 NA
#> R5            1             0.376          2                 NA
#> R8            1             0.376          2                 NA
#> R9            1             0.376          2                 NA
#> R12           1             0.376          2                 NA
#> R17           1             0.376          2                 NA
#> R20           1             0.376          2                 NA
#> R21           1             0.376          2                 NA
#> R24           1             0.376          2                 NA
#> R25           1             0.376          2                 NA
#> R29           1             0.376          2                 NA
#> R36           1             0.376          2                 NA
#> R37           1             0.376          2                 NA
#> R40           1             0.376          2                 NA
#> R48           1             0.376          2                 NA
#> R44           1             0.376          2                 NA
#>     4: share.of.mobility 4: Density 4: Nodes 4: Max.path 4: share.of.total
#> R3                    NA  0.4827586       30           1             0.624
#> R11                   NA  0.4827586       30           1             0.624
#> R23                   NA  0.4827586       30           1             0.624
#> R27                   NA  0.4827586       30           1             0.624
#> R39                   NA  0.4827586       30           1             0.624
#> R7                    NA  0.4827586       30           1             0.624
#> R15                   NA  0.4827586       30           1             0.624
#> R19                   NA  0.4827586       30           1             0.624
#> R31                   NA  0.4827586       30           1             0.624
#> R35                   NA  0.4827586       30           1             0.624
#> R43                   NA  0.4827586       30           1             0.624
#> R47                   NA  0.4827586       30           1             0.624
#> R51                   NA  0.4827586       30           1             0.624
#> R55                   NA  0.4827586       30           1             0.624
#> R59                   NA  0.4827586       30           1             0.624
#> R6                    NA  0.4827586       30           1             0.624
#> R10                   NA  0.4827586       30           1             0.624
#> R14                   NA  0.4827586       30           1             0.624
#> R18                   NA  0.4827586       30           1             0.624
#> R38                   NA  0.4827586       30           1             0.624
#> R42                   NA  0.4827586       30           1             0.624
#> R54                   NA  0.4827586       30           1             0.624
#> R2                    NA  0.4827586       30           1             0.624
#> R22                   NA  0.4827586       30           1             0.624
#> R26                   NA  0.4827586       30           1             0.624
#> R34                   NA  0.4827586       30           1             0.624
#> R50                   NA  0.4827586       30           1             0.624
#> R58                   NA  0.4827586       30           1             0.624
#> R30                   NA  0.4827586       30           1             0.624
#> R46                   NA  0.4827586       30           1             0.624
#> R41                   NA  1.0000000       30           1             0.376
#> R1                    NA  1.0000000       30           1             0.376
#> R4                    NA  1.0000000       30           1             0.376
#> R13                   NA  1.0000000       30           1             0.376
#> R16                   NA  1.0000000       30           1             0.376
#> R28                   NA  1.0000000       30           1             0.376
#> R32                   NA  1.0000000       30           1             0.376
#> R33                   NA  1.0000000       30           1             0.376
#> R45                   NA  1.0000000       30           1             0.376
#> R49                   NA  1.0000000       30           1             0.376
#> R52                   NA  1.0000000       30           1             0.376
#> R53                   NA  1.0000000       30           1             0.376
#> R56                   NA  1.0000000       30           1             0.376
#> R57                   NA  1.0000000       30           1             0.376
#> R60                   NA  1.0000000       30           1             0.376
#> R5                    NA  1.0000000       30           1             0.376
#> R8                    NA  1.0000000       30           1             0.376
#> R9                    NA  1.0000000       30           1             0.376
#> R12                   NA  1.0000000       30           1             0.376
#> R17                   NA  1.0000000       30           1             0.376
#> R20                   NA  1.0000000       30           1             0.376
#> R21                   NA  1.0000000       30           1             0.376
#> R24                   NA  1.0000000       30           1             0.376
#> R25                   NA  1.0000000       30           1             0.376
#> R29                   NA  1.0000000       30           1             0.376
#> R36                   NA  1.0000000       30           1             0.376
#> R37                   NA  1.0000000       30           1             0.376
#> R40                   NA  1.0000000       30           1             0.376
#> R48                   NA  1.0000000       30           1             0.376
#> R44                   NA  1.0000000       30           1             0.376
#>     5: Segment 5: within.mobility 5: share.of.mobility 5: Density 5: Nodes
#> R3           1                  1                    1  0.3644068       60
#> R11          1                  1                    1  0.3644068       60
#> R23          1                  1                    1  0.3644068       60
#> R27          1                  1                    1  0.3644068       60
#> R39          1                  1                    1  0.3644068       60
#> R7           1                  1                    1  0.3644068       60
#> R15          1                  1                    1  0.3644068       60
#> R19          1                  1                    1  0.3644068       60
#> R31          1                  1                    1  0.3644068       60
#> R35          1                  1                    1  0.3644068       60
#> R43          1                  1                    1  0.3644068       60
#> R47          1                  1                    1  0.3644068       60
#> R51          1                  1                    1  0.3644068       60
#> R55          1                  1                    1  0.3644068       60
#> R59          1                  1                    1  0.3644068       60
#> R6           1                  1                    1  0.3644068       60
#> R10          1                  1                    1  0.3644068       60
#> R14          1                  1                    1  0.3644068       60
#> R18          1                  1                    1  0.3644068       60
#> R38          1                  1                    1  0.3644068       60
#> R42          1                  1                    1  0.3644068       60
#> R54          1                  1                    1  0.3644068       60
#> R2           1                  1                    1  0.3644068       60
#> R22          1                  1                    1  0.3644068       60
#> R26          1                  1                    1  0.3644068       60
#> R34          1                  1                    1  0.3644068       60
#> R50          1                  1                    1  0.3644068       60
#> R58          1                  1                    1  0.3644068       60
#> R30          1                  1                    1  0.3644068       60
#> R46          1                  1                    1  0.3644068       60
#> R41          1                  1                    1  0.3644068       60
#> R1           1                  1                    1  0.3644068       60
#> R4           1                  1                    1  0.3644068       60
#> R13          1                  1                    1  0.3644068       60
#> R16          1                  1                    1  0.3644068       60
#> R28          1                  1                    1  0.3644068       60
#> R32          1                  1                    1  0.3644068       60
#> R33          1                  1                    1  0.3644068       60
#> R45          1                  1                    1  0.3644068       60
#> R49          1                  1                    1  0.3644068       60
#> R52          1                  1                    1  0.3644068       60
#> R53          1                  1                    1  0.3644068       60
#> R56          1                  1                    1  0.3644068       60
#> R57          1                  1                    1  0.3644068       60
#> R60          1                  1                    1  0.3644068       60
#> R5           1                  1                    1  0.3644068       60
#> R8           1                  1                    1  0.3644068       60
#> R9           1                  1                    1  0.3644068       60
#> R12          1                  1                    1  0.3644068       60
#> R17          1                  1                    1  0.3644068       60
#> R20          1                  1                    1  0.3644068       60
#> R21          1                  1                    1  0.3644068       60
#> R24          1                  1                    1  0.3644068       60
#> R25          1                  1                    1  0.3644068       60
#> R29          1                  1                    1  0.3644068       60
#> R36          1                  1                    1  0.3644068       60
#> R37          1                  1                    1  0.3644068       60
#> R40          1                  1                    1  0.3644068       60
#> R48          1                  1                    1  0.3644068       60
#> R44          1                  1                    1  0.3644068       60
#>     5: Max.path 5: share.of.total
#> R3            1                 1
#> R11           1                 1
#> R23           1                 1
#> R27           1                 1
#> R39           1                 1
#> R7            1                 1
#> R15           1                 1
#> R19           1                 1
#> R31           1                 1
#> R35           1                 1
#> R43           1                 1
#> R47           1                 1
#> R51           1                 1
#> R55           1                 1
#> R59           1                 1
#> R6            1                 1
#> R10           1                 1
#> R14           1                 1
#> R18           1                 1
#> R38           1                 1
#> R42           1                 1
#> R54           1                 1
#> R2            1                 1
#> R22           1                 1
#> R26           1                 1
#> R34           1                 1
#> R50           1                 1
#> R58           1                 1
#> R30           1                 1
#> R46           1                 1
#> R41           1                 1
#> R1            1                 1
#> R4            1                 1
#> R13           1                 1
#> R16           1                 1
#> R28           1                 1
#> R32           1                 1
#> R33           1                 1
#> R45           1                 1
#> R49           1                 1
#> R52           1                 1
#> R53           1                 1
#> R56           1                 1
#> R57           1                 1
#> R60           1                 1
#> R5            1                 1
#> R8            1                 1
#> R9            1                 1
#> R12           1                 1
#> R17           1                 1
#> R20           1                 1
#> R21           1                 1
#> R24           1                 1
#> R25           1                 1
#> R29           1                 1
#> R36           1                 1
#> R37           1                 1
#> R40           1                 1
#> R48           1                 1
#> R44           1                 1
moneca::segment.quality(out$cols)
#>     Membership 1: Segment 1: within.mobility 1: share.of.mobility 1: Density
#> C5         5.1          5                 NA                   NA        NaN
#> C8         5.1          8                 NA                   NA        NaN
#> C2         5.1          2                 NA                   NA        NaN
#> C14        5.1         14                 NA                   NA        NaN
#> C17        5.1         17                 NA                   NA        NaN
#> C20        5.1         20                 NA                   NA        NaN
#> C23        5.1         23                 NA                   NA        NaN
#> C32        5.1         32                 NA                   NA        NaN
#> C35        5.1         35                 NA                   NA        NaN
#> C38        5.1         38                 NA                   NA        NaN
#> C11        5.1         11                 NA                   NA        NaN
#> C26        5.1         26                 NA                   NA        NaN
#> C29        5.1         29                 NA                   NA        NaN
#> C9         5.1          9                 NA                   NA        NaN
#> C3         5.1          3                 NA                   NA        NaN
#> C6         5.1          6                 NA                   NA        NaN
#> C12        5.1         12                 NA                   NA        NaN
#> C18        5.1         18                 NA                   NA        NaN
#> C21        5.1         21                 NA                   NA        NaN
#> C30        5.1         30                 NA                   NA        NaN
#> C33        5.1         33                 NA                   NA        NaN
#> C36        5.1         36                 NA                   NA        NaN
#> C39        5.1         39                 NA                   NA        NaN
#> C15        5.1         15                 NA                   NA        NaN
#> C24        5.1         24                 NA                   NA        NaN
#> C27        5.1         27                 NA                   NA        NaN
#> C1         5.2          1                 NA                   NA        NaN
#> C4         5.2          4                 NA                   NA        NaN
#> C7         5.2          7                 NA                   NA        NaN
#> C10        5.2         10                 NA                   NA        NaN
#> C16        5.2         16                 NA                   NA        NaN
#> C19        5.2         19                 NA                   NA        NaN
#> C37        5.2         37                 NA                   NA        NaN
#> C40        5.2         40                 NA                   NA        NaN
#> C13        5.2         13                 NA                   NA        NaN
#> C22        5.2         22                 NA                   NA        NaN
#> C25        5.2         25                 NA                   NA        NaN
#> C28        5.2         28                 NA                   NA        NaN
#> C31        5.2         31                 NA                   NA        NaN
#> C34        5.2         34                 NA                   NA        NaN
#>     1: Nodes 1: Max.path 1: share.of.total 2: Segment 2: within.mobility
#> C5         1           0             0.028          1                 NA
#> C8         1           0             0.028          1                 NA
#> C2         1           0             0.027          1                 NA
#> C14        1           0             0.027          1                 NA
#> C17        1           0             0.027          1                 NA
#> C20        1           0             0.027          1                 NA
#> C23        1           0             0.027          1                 NA
#> C32        1           0             0.027          1                 NA
#> C35        1           0             0.027          1                 NA
#> C38        1           0             0.027          1                 NA
#> C11        1           0             0.026          1                 NA
#> C26        1           0             0.026          1                 NA
#> C29        1           0             0.026          1                 NA
#> C9         1           0             0.028          2                 NA
#> C3         1           0             0.027          2                 NA
#> C6         1           0             0.027          2                 NA
#> C12        1           0             0.027          2                 NA
#> C18        1           0             0.027          2                 NA
#> C21        1           0             0.027          2                 NA
#> C30        1           0             0.027          2                 NA
#> C33        1           0             0.027          2                 NA
#> C36        1           0             0.027          2                 NA
#> C39        1           0             0.027          2                 NA
#> C15        1           0             0.026          2                 NA
#> C24        1           0             0.026          2                 NA
#> C27        1           0             0.026          2                 NA
#> C1         1           0             0.022          3                 NA
#> C4         1           0             0.022          3                 NA
#> C7         1           0             0.022          3                 NA
#> C10        1           0             0.022          3                 NA
#> C16        1           0             0.022          3                 NA
#> C19        1           0             0.022          3                 NA
#> C37        1           0             0.022          3                 NA
#> C40        1           0             0.022          3                 NA
#> C13        1           0             0.021          3                 NA
#> C22        1           0             0.021          3                 NA
#> C25        1           0             0.021          3                 NA
#> C28        1           0             0.021          3                 NA
#> C31        1           0             0.021          3                 NA
#> C34        1           0             0.021          3                 NA
#>     2: share.of.mobility 2: Density 2: Nodes 2: Max.path 2: share.of.total
#> C5                    NA          1       13           1             0.349
#> C8                    NA          1       13           1             0.349
#> C2                    NA          1       13           1             0.349
#> C14                   NA          1       13           1             0.349
#> C17                   NA          1       13           1             0.349
#> C20                   NA          1       13           1             0.349
#> C23                   NA          1       13           1             0.349
#> C32                   NA          1       13           1             0.349
#> C35                   NA          1       13           1             0.349
#> C38                   NA          1       13           1             0.349
#> C11                   NA          1       13           1             0.349
#> C26                   NA          1       13           1             0.349
#> C29                   NA          1       13           1             0.349
#> C9                    NA          1       13           1             0.348
#> C3                    NA          1       13           1             0.348
#> C6                    NA          1       13           1             0.348
#> C12                   NA          1       13           1             0.348
#> C18                   NA          1       13           1             0.348
#> C21                   NA          1       13           1             0.348
#> C30                   NA          1       13           1             0.348
#> C33                   NA          1       13           1             0.348
#> C36                   NA          1       13           1             0.348
#> C39                   NA          1       13           1             0.348
#> C15                   NA          1       13           1             0.348
#> C24                   NA          1       13           1             0.348
#> C27                   NA          1       13           1             0.348
#> C1                    NA          1       14           1             0.304
#> C4                    NA          1       14           1             0.304
#> C7                    NA          1       14           1             0.304
#> C10                   NA          1       14           1             0.304
#> C16                   NA          1       14           1             0.304
#> C19                   NA          1       14           1             0.304
#> C37                   NA          1       14           1             0.304
#> C40                   NA          1       14           1             0.304
#> C13                   NA          1       14           1             0.304
#> C22                   NA          1       14           1             0.304
#> C25                   NA          1       14           1             0.304
#> C28                   NA          1       14           1             0.304
#> C31                   NA          1       14           1             0.304
#> C34                   NA          1       14           1             0.304
#>     3: Segment 3: within.mobility 3: share.of.mobility 3: Density 3: Nodes
#> C5           1                 NA                   NA          1       13
#> C8           1                 NA                   NA          1       13
#> C2           1                 NA                   NA          1       13
#> C14          1                 NA                   NA          1       13
#> C17          1                 NA                   NA          1       13
#> C20          1                 NA                   NA          1       13
#> C23          1                 NA                   NA          1       13
#> C32          1                 NA                   NA          1       13
#> C35          1                 NA                   NA          1       13
#> C38          1                 NA                   NA          1       13
#> C11          1                 NA                   NA          1       13
#> C26          1                 NA                   NA          1       13
#> C29          1                 NA                   NA          1       13
#> C9           2                 NA                   NA          1       13
#> C3           2                 NA                   NA          1       13
#> C6           2                 NA                   NA          1       13
#> C12          2                 NA                   NA          1       13
#> C18          2                 NA                   NA          1       13
#> C21          2                 NA                   NA          1       13
#> C30          2                 NA                   NA          1       13
#> C33          2                 NA                   NA          1       13
#> C36          2                 NA                   NA          1       13
#> C39          2                 NA                   NA          1       13
#> C15          2                 NA                   NA          1       13
#> C24          2                 NA                   NA          1       13
#> C27          2                 NA                   NA          1       13
#> C1           3                 NA                   NA          1       14
#> C4           3                 NA                   NA          1       14
#> C7           3                 NA                   NA          1       14
#> C10          3                 NA                   NA          1       14
#> C16          3                 NA                   NA          1       14
#> C19          3                 NA                   NA          1       14
#> C37          3                 NA                   NA          1       14
#> C40          3                 NA                   NA          1       14
#> C13          3                 NA                   NA          1       14
#> C22          3                 NA                   NA          1       14
#> C25          3                 NA                   NA          1       14
#> C28          3                 NA                   NA          1       14
#> C31          3                 NA                   NA          1       14
#> C34          3                 NA                   NA          1       14
#>     3: Max.path 3: share.of.total 4: Segment 4: within.mobility
#> C5            1             0.349          1                 NA
#> C8            1             0.349          1                 NA
#> C2            1             0.349          1                 NA
#> C14           1             0.349          1                 NA
#> C17           1             0.349          1                 NA
#> C20           1             0.349          1                 NA
#> C23           1             0.349          1                 NA
#> C32           1             0.349          1                 NA
#> C35           1             0.349          1                 NA
#> C38           1             0.349          1                 NA
#> C11           1             0.349          1                 NA
#> C26           1             0.349          1                 NA
#> C29           1             0.349          1                 NA
#> C9            1             0.348          1                 NA
#> C3            1             0.348          1                 NA
#> C6            1             0.348          1                 NA
#> C12           1             0.348          1                 NA
#> C18           1             0.348          1                 NA
#> C21           1             0.348          1                 NA
#> C30           1             0.348          1                 NA
#> C33           1             0.348          1                 NA
#> C36           1             0.348          1                 NA
#> C39           1             0.348          1                 NA
#> C15           1             0.348          1                 NA
#> C24           1             0.348          1                 NA
#> C27           1             0.348          1                 NA
#> C1            1             0.304          2                 NA
#> C4            1             0.304          2                 NA
#> C7            1             0.304          2                 NA
#> C10           1             0.304          2                 NA
#> C16           1             0.304          2                 NA
#> C19           1             0.304          2                 NA
#> C37           1             0.304          2                 NA
#> C40           1             0.304          2                 NA
#> C13           1             0.304          2                 NA
#> C22           1             0.304          2                 NA
#> C25           1             0.304          2                 NA
#> C28           1             0.304          2                 NA
#> C31           1             0.304          2                 NA
#> C34           1             0.304          2                 NA
#>     4: share.of.mobility 4: Density 4: Nodes 4: Max.path 4: share.of.total
#> C5                    NA       0.48       26           1             0.696
#> C8                    NA       0.48       26           1             0.696
#> C2                    NA       0.48       26           1             0.696
#> C14                   NA       0.48       26           1             0.696
#> C17                   NA       0.48       26           1             0.696
#> C20                   NA       0.48       26           1             0.696
#> C23                   NA       0.48       26           1             0.696
#> C32                   NA       0.48       26           1             0.696
#> C35                   NA       0.48       26           1             0.696
#> C38                   NA       0.48       26           1             0.696
#> C11                   NA       0.48       26           1             0.696
#> C26                   NA       0.48       26           1             0.696
#> C29                   NA       0.48       26           1             0.696
#> C9                    NA       0.48       26           1             0.696
#> C3                    NA       0.48       26           1             0.696
#> C6                    NA       0.48       26           1             0.696
#> C12                   NA       0.48       26           1             0.696
#> C18                   NA       0.48       26           1             0.696
#> C21                   NA       0.48       26           1             0.696
#> C30                   NA       0.48       26           1             0.696
#> C33                   NA       0.48       26           1             0.696
#> C36                   NA       0.48       26           1             0.696
#> C39                   NA       0.48       26           1             0.696
#> C15                   NA       0.48       26           1             0.696
#> C24                   NA       0.48       26           1             0.696
#> C27                   NA       0.48       26           1             0.696
#> C1                    NA       1.00       14           1             0.304
#> C4                    NA       1.00       14           1             0.304
#> C7                    NA       1.00       14           1             0.304
#> C10                   NA       1.00       14           1             0.304
#> C16                   NA       1.00       14           1             0.304
#> C19                   NA       1.00       14           1             0.304
#> C37                   NA       1.00       14           1             0.304
#> C40                   NA       1.00       14           1             0.304
#> C13                   NA       1.00       14           1             0.304
#> C22                   NA       1.00       14           1             0.304
#> C25                   NA       1.00       14           1             0.304
#> C28                   NA       1.00       14           1             0.304
#> C31                   NA       1.00       14           1             0.304
#> C34                   NA       1.00       14           1             0.304
#>     5: Segment 5: within.mobility 5: share.of.mobility 5: Density 5: Nodes
#> C5           1                 NA                   NA       0.48       26
#> C8           1                 NA                   NA       0.48       26
#> C2           1                 NA                   NA       0.48       26
#> C14          1                 NA                   NA       0.48       26
#> C17          1                 NA                   NA       0.48       26
#> C20          1                 NA                   NA       0.48       26
#> C23          1                 NA                   NA       0.48       26
#> C32          1                 NA                   NA       0.48       26
#> C35          1                 NA                   NA       0.48       26
#> C38          1                 NA                   NA       0.48       26
#> C11          1                 NA                   NA       0.48       26
#> C26          1                 NA                   NA       0.48       26
#> C29          1                 NA                   NA       0.48       26
#> C9           1                 NA                   NA       0.48       26
#> C3           1                 NA                   NA       0.48       26
#> C6           1                 NA                   NA       0.48       26
#> C12          1                 NA                   NA       0.48       26
#> C18          1                 NA                   NA       0.48       26
#> C21          1                 NA                   NA       0.48       26
#> C30          1                 NA                   NA       0.48       26
#> C33          1                 NA                   NA       0.48       26
#> C36          1                 NA                   NA       0.48       26
#> C39          1                 NA                   NA       0.48       26
#> C15          1                 NA                   NA       0.48       26
#> C24          1                 NA                   NA       0.48       26
#> C27          1                 NA                   NA       0.48       26
#> C1           2                 NA                   NA       1.00       14
#> C4           2                 NA                   NA       1.00       14
#> C7           2                 NA                   NA       1.00       14
#> C10          2                 NA                   NA       1.00       14
#> C16          2                 NA                   NA       1.00       14
#> C19          2                 NA                   NA       1.00       14
#> C37          2                 NA                   NA       1.00       14
#> C40          2                 NA                   NA       1.00       14
#> C13          2                 NA                   NA       1.00       14
#> C22          2                 NA                   NA       1.00       14
#> C25          2                 NA                   NA       1.00       14
#> C28          2                 NA                   NA       1.00       14
#> C31          2                 NA                   NA       1.00       14
#> C34          2                 NA                   NA       1.00       14
#>     5: Max.path 5: share.of.total
#> C5            1             0.696
#> C8            1             0.696
#> C2            1             0.696
#> C14           1             0.696
#> C17           1             0.696
#> C20           1             0.696
#> C23           1             0.696
#> C32           1             0.696
#> C35           1             0.696
#> C38           1             0.696
#> C11           1             0.696
#> C26           1             0.696
#> C29           1             0.696
#> C9            1             0.696
#> C3            1             0.696
#> C6            1             0.696
#> C12           1             0.696
#> C18           1             0.696
#> C21           1             0.696
#> C30           1             0.696
#> C33           1             0.696
#> C36           1             0.696
#> C39           1             0.696
#> C15           1             0.696
#> C24           1             0.696
#> C27           1             0.696
#> C1            1             0.304
#> C4            1             0.304
#> C7            1             0.304
#> C10           1             0.304
#> C16           1             0.304
#> C19           1             0.304
#> C37           1             0.304
#> C40           1             0.304
#> C13           1             0.304
#> C22           1             0.304
#> C25           1             0.304
#> C28           1             0.304
#> C31           1             0.304
#> C34           1             0.304
```

``` r
if (requireNamespace("ggraph", quietly = TRUE)) {
  moneca::plot_moneca_ggraph(
    out$rows,
    level = 2,
    node_color = "segment",
    title = "moneca_bipartite - row view"
  )
}
```

![Row-view network
plot](monecascale-bipartite_files/figure-html/downstream-ggraph-rows-1.png)

``` r
if (requireNamespace("ggraph", quietly = TRUE)) {
  moneca::plot_moneca_hierarchical(
    out$cols,
    title = "moneca_bipartite - column view"
  )
}
```

![Column-view hierarchical
plot](monecascale-bipartite_files/figure-html/downstream-hier-cols-1.png)

The joint row-block by column-block view is not a moneca-class object
and is plotted directly from the diagnostics slot:

``` r
image(
  t(bim[[1]])[, nrow(bim[[1]]):1],
  axes = FALSE,
  main = "Joint block interaction (finest level)",
  xlab = "column block",
  ylab = "row block"
)
```

![Joint block-interaction
heatmap](monecascale-bipartite_files/figure-html/joint-heatmap-1.png)

## 6. Latent-block recovery experiment

A block-constant generator with known row and column memberships lets us
measure how well DcLbm recovers the planted structure. We compare the
row-side and column-side memberships at the level closest to the planted
sizes against `row_block` and `col_block` with the adjusted Rand index.
The fixture is independent of the quick-start block to keep this section
self-contained.

``` r
# 1. planted design -----
set.seed(42)
n_row <- 80
n_col <- 60
K_row_true <- 4
K_col_true <- 3
row_block_true <- rep(seq_len(K_row_true), length.out = n_row)
col_block_true <- rep(seq_len(K_col_true), length.out = n_col)

# 2. block-constant Poisson with strong diagonal signal -----
omega <- matrix(0.8, K_row_true, K_col_true)
for (g in seq_len(K_row_true)) {
  h <- ((g - 1) %% K_col_true) + 1L
  omega[g, h] <- 8
}
lambda <- omega[row_block_true, col_block_true]
mx_parity <- matrix(
  stats::rpois(n_row * n_col, lambda = lambda),
  n_row,
  n_col
)
rownames(mx_parity) <- paste0("R", seq_len(n_row))
colnames(mx_parity) <- paste0("C", seq_len(n_col))

# 3. fit -----
fit_parity <- monecascale::moneca_bipartite(
  mx_parity,
  seed = 42,
  verbose = FALSE
)
```

``` r
# Minimal ARI. Returns 1 when labelings agree up to permutation.
ari <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  sc <- function(x) sum(choose(x, 2))
  a_sum <- sc(rowSums(tab))
  b_sum <- sc(colSums(tab))
  tab_sum <- sc(as.vector(tab))
  expected <- a_sum * b_sum / choose(n, 2)
  max_idx <- (a_sum + b_sum) / 2
  if (max_idx == expected) return(1)
  (tab_sum - expected) / (max_idx - expected)
}

level_labels <- function(mon, level) {
  n_core <- nrow(mon$mat.list[[1]]) - 1L
  labels <- rep(NA_integer_, n_core)
  cliques <- mon$segment.list[[level]]
  for (b in seq_along(cliques)) labels[cliques[[b]]] <- b
  na_idx <- is.na(labels)
  if (any(na_idx)) {
    labels[na_idx] <- max(0L, labels, na.rm = TRUE) + seq_len(sum(na_idx))
  }
  as.integer(labels)
}
```

``` r
# 1. pick the level whose row-block count matches the planted K_row -----
krow_per_level <- fit_parity$bipartite_diagnostics$n_blocks_row_per_level
kcol_per_level <- fit_parity$bipartite_diagnostics$n_blocks_col_per_level
row_level <- which.min(abs(krow_per_level - K_row_true)) + 1L
col_level <- which.min(abs(kcol_per_level - K_col_true)) + 1L

# 2. extract moneca-style node labels at those levels -----
row_hat <- level_labels(fit_parity$rows, row_level)
col_hat <- level_labels(fit_parity$cols, col_level)

# 3. score against planted memberships -----
ari_row <- ari(row_block_true, row_hat)
ari_col <- ari(col_block_true, col_hat)
data.frame(
  side = c("rows", "cols"),
  K_true = c(K_row_true, K_col_true),
  K_hat = c(krow_per_level[row_level - 1L], kcol_per_level[col_level - 1L]),
  ARI = round(c(ari_row, ari_col), 3)
)
#>   side K_true K_hat   ARI
#> 1 rows      4     3 0.706
#> 2 cols      3     3 1.000
```

On this fixture DcLbm typically recovers both sides with ARI \> 0.5, and
often well above. Exact values vary with the seed and with the
signal-to-noise contrast between diagonal and off-diagonal block
intensities. The load-bearing claim is directional: the rectangular RR =
DcLbm residual bridge is empirically supported on data where the planted
structure is known.

## 7. Limitations and future directions

- **One-mode projections, not joint clustering.** Exposing the fit to
  moneca requires two square objects, so the joint `Krow x Kcol`
  co-clustering is not representable inside `$rows` or `$cols`. The full
  block-interaction matrices are retained in
  `$bipartite_diagnostics$block_interaction_matrix_per_level` for
  consumers that need the joint view.
- **No isolates summary on the bipartite side yet.** Row- and
  column-side isolates surfaces are deferred to release 0.2.1.
- **Single backend.**
  [`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md)
  currently exposes the `greed` DcLbm backend only. A `graph-tool`
  bipartite backend via `pclabel` is reserved for input sizes beyond
  ~10^5 per side.
- **Auto-level composition pending.** Until direction **D6**
  (`auto_segment_levels()`) lands, cap the hierarchy explicitly via
  `segment.levels` when downstream consumers expect a fixed depth.
- **Square input is rejected.** For square mobility matrices use
  [`monecascale::moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)
  or
  [`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html).

## 8. References

- Come, E., Jouvin, N., Latouche, P., Bouveyron, C., Bataillon, E., &
  Chen, A. (2021). Hierarchical clustering with discrete latent variable
  models and the integrated classification likelihood. *Advances in Data
  Analysis and Classification*, 15, 957-986.
- Larremore, D. B., Clauset, A., & Jacobs, A. Z. (2014). Efficiently
  inferring community structure in bipartite networks. *Physical Review
  E*, 90(1), 012805.
- Karrer, B., & Newman, M. E. J. (2011). Stochastic blockmodels and
  community structure in networks. *Physical Review E*, 83, 016107.
- Peixoto, T. P. (2014). Hierarchical block structures and
  high-resolution model selection in large networks. *Physical Review
  X*, 4, 011047.
