# monecascale: Scalable Backends for MONECA

**Author: Giampaolo Montaletti** ORCID:
<https://orcid.org/0009-0002-5327-1122> Email:
<giampaolo.montaletti@gmail.com>

`monecascale` is a sibling package to
[**moneca**](https://github.com/gmontaletti/MONECA) that provides
scalable clustering backends for mobility network analysis. It is
designed to extend moneca beyond classification-scale matrices (dozens
to a few hundred categories) toward settings such as person × employer
mobility (10^6 × 10^5).

## Relationship to moneca

- `moneca` remains the home of the reference algorithm
  ([`moneca()`](https://gmontaletti.github.io/MONECA/reference/moneca.html),
  [`moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html)),
  weight-matrix construction, segment metadata, and the full
  visualization stack.
- `monecascale` adds alternative clustering backends whose output is
  consumed unchanged by moneca’s analysis and plotting functions. Every
  `monecascale` entry point returns a `moneca`-class S3 object with an
  extra diagnostics slot.

## Current backends

### `moneca_sbm()`

Hierarchical degree-corrected stochastic block model (DC-SBM)
clustering. Two dispatchable backends:

- `"greed"` (default, R-native via
  [`greed::DcSbm`](https://comeetie.github.io/greed/reference/DcSbm.html))
  — sparse `dgCMatrix` end-to-end, hierarchical dendrogram via
  `cut(sol, k)`, portable parallel via `future`. Scales to ~10^(4–10)5
  nodes.
- `"graphtool"` (opt-in, Python via `reticulate`) — Peixoto’s
  `minimize_nested_blockmodel_dl()`, MDL-based level selection. Scales
  to 10^6+ nodes.

**Theoretical rationale.** moneca’s `RR = O / E` is the DC-SBM residual
under identity block interaction; DC-SBM inference is therefore a
principled scalable drop-in for moneca’s clique step. Empirical parity
on a 20-class synthetic fixture: `ARI = 0.671`, `NMI = 0.905` vs
[`moneca::moneca_fast()`](https://gmontaletti.github.io/MONECA/reference/moneca_fast.html).

### `moneca_bipartite()`

Degree-corrected Latent Block Model (DcLbm) for rectangular person ×
employer (or worker × occupation) mobility matrices. Output comprises
two linked `moneca`-class objects (`$rows`, `$cols`) consumable by
moneca’s full analysis and plotting stack. Scales to ~10^(4–10)5
rows/columns per side.

### `moneca_flow()`

Flow-based hierarchical clustering via Infomap and the Map Equation.
Backend via
[`igraph::cluster_infomap()`](https://r.igraph.org/reference/cluster_infomap.html)
(R-native, scales to ~10^5 nodes). A recursive-flat wrapper synthesises
a 2–3 level hierarchy. Returns a `moneca`-class object with
`$flow_diagnostics` (codelength per level, node assignments).

### `auto_segment_levels()`

Post-hoc hierarchy-level auto-selection from a full hierarchical fit.
Two criteria:

- `method = "mdl"` — MDL elbow detection via kneedle algorithm (Satopää
  et al. 2011), with `max_second_diff` fallback.
- `method = "mi_plateau"` — mutual-information plateau: picks the level
  where marginal MI gain drops below a threshold of the maximum MI.

Available as `segment.levels = "auto"` on
[`moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md),
[`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md),
and
[`moneca_flow()`](https://gmontaletti.github.io/monecascale/reference/moneca_flow.md)
for one-shot hierarchy fitting and trimming. Recognises backend-specific
diagnostics (MDL for SBM and flow; mutual information for bipartite).

### `rr_from_duckdb()`

Out-of-core sparse relative-risk matrix construction via DuckDB SQL.
Accepts edge-list input as a data.frame, a file path (`.csv` /
`.parquet`), or a table name on a pre-connected DuckDB instance. Returns
a `monecascale_rr` list with sparse `$rr`, sparse `$counts`, and margin
metadata. Composes seamlessly with
[`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md)
(via `$rr`) and
[`moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)
/
[`moneca_flow()`](https://gmontaletti.github.io/monecascale/reference/moneca_flow.md)
(via `$counts`).

## Installation

``` r
# Install moneca first (required)
devtools::install_github("gmontaletti/MONECA")

# Then monecascale
devtools::install_github("gmontaletti/monecascale")

# Optional backend dependencies
install.packages("greed")                              # default R-native backend
install.packages("reticulate")                         # for graph-tool backend
monecascale::moneca_sbm_install_graphtool()            # installs graph-tool via conda
```

## Quick start

``` r
library(monecascale)
library(moneca)

mob <- generate_mobility_data(
  n_classes = 20, n_total = 5000,
  immobility_strength = 0.5, class_clustering = 0.75,
  seed = 2026
)

seg <- moneca_sbm(mob, backend = "greed", seed = 2026)

seg$sbm_diagnostics$n_blocks_per_level
segment.membership(seg)
plot_moneca_hierarchical(seg)
```

See
[`vignette("monecascale-sbm")`](https://gmontaletti.github.io/monecascale/articles/monecascale-sbm.md)
for the RR ≡ DC-SBM residual walkthrough, the parity experiment, and
guidance on choosing a backend.

## Roadmap

`monecascale` is the home for the following scaling directions from the
[MONECA Scaling
Roadmap](https://github.com/gmontaletti/MONECA/blob/master/reference/moneca/SCALING_ROADMAP.md):

- **D2** — hierarchical DC-SBM (`moneca_sbm`, *shipped*).
- **D1** — bipartite / two-mode (`moneca_bipartite`, *shipped*).
- **D6** — auto-level detection (`auto_segment_levels`, *shipped*).
- **D3** — flow-based (`moneca_flow` via Infomap, *shipped*).
- **D4** — scalable clique guardrail (`moneca_localclique`, planned).
- **D5** — Poisson NMF (`moneca_nmf`, planned).
- **D7** — out-of-core
  ([`rr_from_duckdb()`](https://gmontaletti.github.io/monecascale/reference/rr_from_duckdb.md),
  *shipped*).

## Citation

    Montaletti, G. (2026). monecascale: Scalable Backends for MONECA.
    R package version 0.5.0. https://github.com/gmontaletti/monecascale

    Montaletti, G. (2026). moneca: Mobility Network Clustering Analysis.
    R package version 1.8.0. https://github.com/gmontaletti/MONECA

    Touboel, J., & Larsen, A. G. (2017). Mapping the Social Class Structure:
    From Occupational Mobility to Social Class Categories Using Network Analysis.
    Sociology, 51(6), 1257-1276.

## License

GPL-3
