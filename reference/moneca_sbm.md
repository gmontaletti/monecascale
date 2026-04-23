# Hierarchical DC-SBM Backend for MONECA Segmentation

Alternative to
[`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html)
that replaces clique enumeration with hierarchical degree-corrected
stochastic block model (DC-SBM) inference. MONECA's weight matrix
`RR = O / E` is the expected-value residual under a DC-SBM with identity
block interaction, so DC-SBM clustering is a principled scalable
replacement for the clique step.

## Usage

``` r
moneca_sbm(
  mx,
  backend = c("auto", "greed", "graphtool"),
  deg_corr = TRUE,
  edge_model = c("poisson", "bernoulli"),
  small.cell.reduction = 0,
  symmetric_method = c("sum", "min", "none"),
  segment.levels = NULL,
  auto_method = c("mdl", "mi_plateau"),
  max_K = NULL,
  isolates = FALSE,
  has_margins = "auto",
  seed = NULL,
  n_init = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- mx:

  A square mobility matrix, optionally with margins as the last row and
  column (auto-detected via `has_margins = "auto"`). Sparse `dgCMatrix`
  input is supported end-to-end through the `greed` backend.

- backend:

  One of `"auto"`, `"greed"`, `"graphtool"`. `"auto"` prefers
  `graphtool` when available and the matrix exceeds 5000 nodes;
  otherwise uses `greed`.

- deg_corr:

  Logical. If `TRUE` (default), fits a degree-corrected SBM. Set to
  `FALSE` for plain SBM.

- edge_model:

  One of `"poisson"` (default, appropriate for count-valued mobility) or
  `"bernoulli"` (binary edges).

- small.cell.reduction:

  Numeric. Cells in the core matrix below this value are zeroed before
  SBM inference. Mirrors the
  [`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html)
  parameter.

- symmetric_method:

  One of `"sum"` (default; `mx + t(mx)`), `"min"`
  (`2 * pmin(mx, t(mx))`, down-weights one-way flows), or `"none"` (keep
  directed). Applied to the core before SBM.

- segment.levels:

  Integer, `NULL`, or the string `"auto"`. If `NULL` (default), all
  hierarchy levels produced by the backend (from K_final down to 2) are
  returned. If an integer, the hierarchy is truncated to that many
  approximately log-spaced levels. If `"auto"`, the full hierarchy is
  fit and then
  [`auto_segment_levels`](https://gmontaletti.github.io/monecascale/reference/auto_segment_levels.md)
  picks a single preferred level via `auto_method`; the object is
  trimmed to that level and carries the picker result under
  `$auto_level`.

- auto_method:

  One of `"mdl"` (default) or `"mi_plateau"`. Ignored unless
  `segment.levels = "auto"`.

- max_K:

  Integer or `NULL`. Upper bound on number of blocks at the finest level
  (greed backend only). Defaults to `max(20, n_core / 5)`.

- isolates:

  Logical. If `TRUE`, attaches an `$isolates_summary` slot mirroring
  [`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html).

- has_margins:

  One of `"auto"` (default), `TRUE`, or `FALSE`. Controls whether to
  detect and/or append row/column margins to the input matrix.

- seed:

  Integer seed for reproducible inference.

- n_init:

  Integer. Number of independent backend runs; the best-scoring fit is
  retained. Defaults to 1 for `greed` (which is quasi-deterministic) and
  5 is typical for `graphtool` (MCMC).

- verbose:

  Logical. If `TRUE`, prints backend progress.

- ...:

  Additional arguments forwarded to the backend.

## Value

An object of class `"moneca"` with the same structure as
[`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html),
plus an extra slot `$sbm_diagnostics` describing the backend used, the
number of blocks per level, and the MDL/ICL trace.

## Details

The number of hierarchical levels is selected automatically by the
backend (ICL for `greed`; MDL for `graphtool`), which eliminates
MONECA's manual `segment.levels` requirement. Passing `segment.levels`
truncates the backend's full hierarchy to that many approximately
log-spaced levels, preserving the finest and coarsest.

## See also

[`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html)
for the clique-based reference implementation;
[`moneca_sbm_install_graphtool`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm_install_graphtool.md)
to set up the `graphtool` backend.

## Examples

``` r
if (FALSE) { # \dontrun{
mob <- moneca::generate_mobility_data(
  n_classes = 20, n_total = 2000, seed = 42
)
seg <- moneca_sbm(mob, backend = "greed", seed = 1L)
moneca::segment.membership(seg)
moneca::plot_moneca_hierarchical(seg)
} # }
```
