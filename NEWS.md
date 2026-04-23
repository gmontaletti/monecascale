# monecascale 0.1.0

## Initial release

Sibling package for [moneca](https://github.com/gmontaletti/MONECA)
providing scalable clustering backends. Extracted from `moneca` branch
`new_directions` on 2026-04-23.

### Features

* `moneca_sbm()` — hierarchical degree-corrected stochastic block model
  clustering as an alternative to `moneca::moneca_fast()`'s clique
  enumeration. Two backends:
    * `"greed"` (default, R-native via `greed::DcSbm`): hierarchical
      dendrogram, Poisson likelihood, sparse `dgCMatrix` end-to-end,
      portable parallel via the `future` framework. Scales to
      ~10^4–10^5 nodes.
    * `"graphtool"` (opt-in, Python via `reticulate`): Peixoto's
      `minimize_nested_blockmodel_dl()` with automatic MDL-based level
      selection. Scales to 10^6+ nodes. Install via
      `moneca_sbm_install_graphtool()`.

* The output of `moneca_sbm()` is a standard `moneca`-class S3 object
  consumed unchanged by `moneca::segment.membership()`,
  `moneca::segment.quality()`, `moneca::plot_moneca_*()` etc. An extra
  `$sbm_diagnostics` slot records backend, MDL/ICL per level, and the
  number of blocks per level.

### Theoretical rationale

`moneca`'s weight matrix `RR = O / E` is structurally the expected-value
residual under a DC-SBM with identity block interaction. Consequence:
`moneca_fast()` is approximately DC-SBM clustering done by clique
enumeration on the thresholded residual, and modern DC-SBM inference
(Peixoto's `graph-tool`, `greed`) is a principled drop-in.

### Empirical validation

On a 20-class synthetic fixture (`generate_mobility_data(n_classes = 20,
immobility_strength = 0.5, class_clustering = 0.75, seed = 2026)`), the
best-matching level pair between `moneca_sbm(backend = "greed")` and
`moneca::moneca_fast()` reaches **ARI = 0.671, NMI = 0.905** — above the
0.55 aspirational target.
