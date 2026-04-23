# Bipartite (Rectangular) DC-LBM Backend for MONECA Segmentation

Rectangular counterpart of \[moneca_sbm()\]. Accepts a two-mode mobility
matrix (rows and columns are distinct entity types, e.g. persons x
employers, workers x occupations) and fits a hierarchical
degree-corrected Latent Block Model (\`greed::DcLbm\`). The rectangular
relative-risk matrix \`RR = O / E\` acts as the DcLbm residual under
identity block interaction, mirroring the square RR = DC-SBM bridge that
underpins \`moneca_sbm()\`.

## Usage

``` r
moneca_bipartite(
  mx,
  backend = "greed",
  deg_corr = TRUE,
  edge_model = "poisson",
  small.cell.reduction = 0,
  margin_policy = c("row_col", "none"),
  segment.levels = NULL,
  max_K = NULL,
  seed = NULL,
  n_init = 1L,
  has_margins = "auto",
  verbose = FALSE,
  ...
)
```

## Arguments

- mx:

  A rectangular mobility matrix of shape \`n_row x n_col\`. Optional
  moneca-style margins in the last row and column are auto-detected when
  \`has_margins = "auto"\`. Sparse \`dgCMatrix\` input is supported end
  to end.

- backend:

  Currently only \`"greed"\` (DcLbm via the \`greed\` package). Argument
  kept for forward compatibility with future bipartite backends.

- deg_corr:

  Logical. Reserved. \`greed::DcLbm\` is always degree corrected; a
  non-default value is silently accepted but has no effect.

- edge_model:

  One of \`"poisson"\` (default, appropriate for count-valued mobility)
  or \`"bernoulli"\`. Reserved. DcLbm uses the Poisson edge model;
  non-default values emit a warning and are otherwise ignored.

- small.cell.reduction:

  Numeric. Cells in the core rectangle with values below this threshold
  are zeroed before RR computation. Mirrors the
  \`moneca::moneca_fast()\` parameter.

- margin_policy:

  One of \`"row_col"\` (default; append margins to each one-mode
  projection, matching moneca's convention) or \`"none"\` (skip
  appending margins). \`"row_col"\` is the load-bearing choice for
  downstream moneca compatibility.

- segment.levels:

  Integer or \`NULL\`. If \`NULL\` (default), all hierarchy levels
  produced by the backend (from the finest joint \`K\` down to \`2\`)
  are returned. If an integer, the hierarchy is truncated to that many
  approximately log-spaced levels.

- max_K:

  Integer or \`NULL\`. Reserved for backend-specific caps; the \`greed\`
  backend currently uses its internal default.

- seed:

  Integer seed for reproducible inference.

- n_init:

  Integer. Number of independent backend runs; the best-scoring fit
  (highest ICL) is retained. Defaults to \`1L\`.

- has_margins:

  One of \`"auto"\` (default), \`TRUE\`, or \`FALSE\`. Controls whether
  to detect and/or append row/column margins to the input matrix.

- verbose:

  Logical. If \`TRUE\`, prints backend progress.

- ...:

  Additional arguments forwarded to the backend.

## Value

An object of class \`"moneca_bipartite"\` with slots: \* \`rows\` - full
moneca-class object for the row view (one-mode projection on rows). \*
\`cols\` - full moneca-class object for the column view. \*
\`bipartite_diagnostics\` - list with backend metadata, per-level row /
column block counts, joint ICL / MDL traces, and per-level \`Krow x
Kcol\` block-interaction matrices computed directly from the rectangular
core.

## Details

Because \`moneca::weight.matrix()\` and the downstream \`segment.\*()\`
stack assume square inputs, the bipartite fit is exposed to moneca
through two linked moneca-class objects built from one-mode projections
of the rectangular RR. The row view uses \`RR column view the transpose,
with \`D_c = diag(col_margin / grand_total)\` (analogously for rows).
Each projection is zero-diagonal, symmetric, sparse, and carries
moneca-style margins appended as the last row and column.

The row view and column view each expose the full moneca output contract
(\`segment.list\`, \`mat.list\`, \`segment_metadata\`, and companion
scalars). \`moneca::segment.membership()\`,
\`moneca::segment.quality()\`, \`moneca::plot_moneca_ggraph()\`, and
\`moneca::plot_moneca_hierarchical()\` run unchanged on \`out\$rows\` or
\`out\$cols\`.

The joint \`Krow x Kcol\` co-clustering is not representable as a single
square moneca object. It is preserved in
\`bipartite_diagnostics\$block_interaction_matrix_per_level\`, which
cross-tabulates the rectangular core by row and column memberships at
each level.

## See also

\[moneca_sbm()\] for the square (one-mode) counterpart;
\`moneca::moneca_fast()\` for the reference clique-based algorithm on
square matrices.

## Examples

``` r
if (requireNamespace("greed", quietly = TRUE)) {
  # Toy rectangular mobility fixture
  set.seed(1)
  rect <- matrix(stats::rpois(40 * 25, lambda = 2), 40, 25)
  rownames(rect) <- paste0("R", seq_len(40))
  colnames(rect) <- paste0("C", seq_len(25))
  out <- moneca_bipartite(rect, seed = 1L, verbose = FALSE)
  moneca::segment.membership(out$rows)
  moneca::segment.membership(out$cols)
}
#>    name membership
#> 1    C1        1.1
#> 2    C2        1.2
#> 3    C3        1.3
#> 4    C4        1.4
#> 5    C5        1.5
#> 6    C6        1.6
#> 7    C7        1.7
#> 8    C8        1.8
#> 9    C9        1.9
#> 10  C10       1.10
#> 11  C11       1.11
#> 12  C12       1.12
#> 13  C13       1.13
#> 14  C14       1.14
#> 15  C15       1.15
#> 16  C16       1.16
#> 17  C17       1.17
#> 18  C18       1.18
#> 19  C19       1.19
#> 20  C20       1.20
#> 21  C21       1.21
#> 22  C22       1.22
#> 23  C23       1.23
#> 24  C24       1.24
#> 25  C25       1.25
```
