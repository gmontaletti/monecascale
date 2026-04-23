# monecascale 0.2.0

## New features

* `moneca_bipartite()` — D1 of the scaling roadmap. Accepts rectangular
  person × employer (or worker × occupation, etc.) mobility matrices and
  returns an S3 object of class `"moneca_bipartite"` wrapping two linked
  `moneca`-class objects (`$rows`, `$cols`) plus a `$bipartite_diagnostics`
  slot carrying the joint ICL/MDL trace, per-level block counts per side,
  rectangular RR, and per-level `Krow × Kcol` block-interaction matrices.
* Single backend: `greed::DcLbm` (hierarchical Poisson latent block model).
  Scales to ~10^4–10^5 rows / cols per side.
* Output contract: each of `$rows` and `$cols` is consumable unchanged by
  `moneca::segment.membership()`, `moneca::segment.quality()`,
  `moneca::plot_moneca_ggraph()`, `moneca::plot_moneca_hierarchical()`.

## Theoretical note

Rectangular RR = O / ((row · col) / total) is the DcLbm residual under
identity block interaction, extending the RR ≡ DC-SBM bridge from D2.
`moneca_bipartite()` fits DcLbm on scaled-rounded RR and projects each
side via one-mode RR aggregation with `D_other_margin^{-1}` normalisation.

## Empirical validation

On a 60 × 40 synthetic bipartite fixture with 3 latent row blocks and
2 latent col blocks (seed 2026):
- rows: ARI = 0.563, NMI = 0.761
- cols: ARI = 1.000, NMI = 1.000

See `tests/testthat/test-bipartite-parity.R`.

## Dependencies

* Added `mclust` to `Suggests` for the parity test's ARI computation.
* No new `Imports`. `greed` remains `Suggests`.

## Open gaps

* Isolates summary on the bipartite side (zero-degree rows / cols) is
  deferred to 0.2.1.
* Graphtool bipartite backend (via `pclabel`) reserved for scale ≥ 10^5
  per side.

---

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
