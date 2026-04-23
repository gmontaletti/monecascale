# Package index

## Scalable Clustering Backends

Alternative backends producing moneca-class output

- [`moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)
  : Hierarchical DC-SBM Backend for MONECA Segmentation

- [`moneca_sbm_install_graphtool()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm_install_graphtool.md)
  :

  Install the graph-tool Python Backend for
  [`moneca_sbm()`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)

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
