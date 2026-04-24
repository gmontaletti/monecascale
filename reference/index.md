# Package index

## Scalable Clustering Backends

Alternative backends producing moneca-class output

- [`moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)
  : Hierarchical DC-SBM Backend for MONECA Segmentation

- [`moneca_sbm_install_graphtool()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm_install_graphtool.md)
  :

  Install the graph-tool Python Backend for
  [`moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)

## Flow Clustering

Infomap / Map Equation on mobility counts. Flat Infomap via igraph, with
a recursive wrapper that synthesises a 2-3 level hierarchy. Output is a
plain moneca-class object with a \$flow_diagnostics slot.

- [`moneca_flow()`](https://gmontaletti.github.io/monecascale/reference/moneca_flow.md)
  : Flow-based (Infomap / Map Equation) Backend for MONECA Segmentation

## Bipartite Segmentation

Rectangular input (person x employer, worker x occupation). Returns two
linked moneca-class objects - one row-view, one col-view - each
consumable by moneca’s segment.*() and plot_moneca\_*() stack.

- [`moneca_bipartite()`](https://gmontaletti.github.io/monecascale/reference/moneca_bipartite.md)
  : Bipartite (Rectangular) DC-LBM Backend for MONECA Segmentation
- [`format(`*`<moneca_bipartite>`*`)`](https://gmontaletti.github.io/monecascale/reference/format.moneca_bipartite.md)
  : Format a \`moneca_bipartite\` object
- [`print(`*`<moneca_bipartite>`*`)`](https://gmontaletti.github.io/monecascale/reference/print.moneca_bipartite.md)
  : Print a \`moneca_bipartite\` object
- [`summary(`*`<moneca_bipartite>`*`)`](https://gmontaletti.github.io/monecascale/reference/summary.moneca_bipartite.md)
  : Summary of a \`moneca_bipartite\` object
- [`print(`*`<summary.moneca_bipartite>`*`)`](https://gmontaletti.github.io/monecascale/reference/print.summary.moneca_bipartite.md)
  : Print a \`summary.moneca_bipartite\` object

## Auto-level Selection

Post-hoc selection of a hierarchy level from a moneca_sbm or
moneca_bipartite fit via MDL elbow or level-to-level MI plateau.
Composes with the segment.levels = “auto” wrapper on both backends.

- [`auto_segment_levels()`](https://gmontaletti.github.io/monecascale/reference/auto_segment_levels.md)
  : Auto-select a Hierarchy Level for a \`moneca\` Fit
- [`print(`*`<auto_segment_levels>`*`)`](https://gmontaletti.github.io/monecascale/reference/print.auto_segment_levels.md)
  : Print an \`auto_segment_levels\` Result
- [`format(`*`<auto_segment_levels>`*`)`](https://gmontaletti.github.io/monecascale/reference/format.auto_segment_levels.md)
  : Format an \`auto_segment_levels\` Result

## Out-of-core RR

Build a sparse relative-risk matrix from an out-of-core edge source
(DuckDB / Parquet / CSV) without densifying to an R matrix. Output
composes with the sparse-aware clustering backends.

- [`rr_from_duckdb()`](https://gmontaletti.github.io/monecascale/reference/rr_from_duckdb.md)
  : Build a Sparse Relative-Risk Matrix Out-of-Core via DuckDB
