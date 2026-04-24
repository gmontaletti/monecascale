# monecascale 0.4.0

## New features

* `moneca_flow()` — D3 of the scaling roadmap. Flow-based hierarchical
  clustering via Infomap and the Map Equation. Single backend:
  `igraph::cluster_infomap()` (R-native, scales to ~10^5 nodes). A
  recursive-flat wrapper synthesises a 2–3 level hierarchy from the flat
  solution. Returns a plain `moneca`-class object with `$flow_diagnostics`
  (codelength per level, node flow assignments).
* `auto_segment_levels()` now recognises flow-backend objects. The MDL
  criterion reads `codelength_per_level` and reports a `codelength` score
  column per level.

## Empirical validation

On a 4-block 60-node synthetic fixture (seed 2026):
- flow: ARI = 1.000, NMI = 1.000

## Documentation

* New vignette `monecascale-flow` covering the Map Equation, the
  flat-vs-hierarchical trade-off, and a comparison of flow recovery
  against `moneca_sbm()` on the same fixture.

## Dependencies

* `igraph (>= 1.3.0)` moved from transitive (via moneca) to explicit
  `Imports`.

---

# monecascale 0.3.0

## New features

* `auto_segment_levels()` — D6 of the scaling roadmap. Post-hoc
  auto-selection of a hierarchy level from a fit, using either:
    * `method = "mdl"` — MDL elbow via the kneedle algorithm
      (Satopää et al. 2011) or `max_second_diff` fallback.
    * `method = "mi_plateau"` — level-to-level mutual information
      plateau, picking the level where marginal gain drops below
      `plateau_tol` of the max MI.
* `segment.levels = "auto"` wrapper on `moneca_sbm()` and
  `moneca_bipartite()` — fits the full hierarchy, calls
  `auto_segment_levels()` post-hoc, trims `segment.list` and
  `mat.list` to the chosen level, and attaches `$auto_level` with
  diagnostics.
* New `auto_method` argument on both backends, defaulting to `"mdl"`.
* Returns an S3 `auto_segment_levels` object with `$level`,
  `$method`, `$diagnostics` (per-level `mdl`, `mi_to_next`,
  `n_blocks`, `score`), `$backend`, and `$call`. Includes
  `print` and `format` methods.

## Internals

* `R/bipartite_utils.R::.build_side_moneca()` now stamps each bipartite
  side with `$bipartite_origin` and a slim `$bipartite_diagnostics_side`,
  so detached sides can feed `auto_segment_levels()` without the parent
  object.

## Coverage

* SBM and bipartite backends fully supported. `moneca_fast()` output
  raises an explicit error pointing at the 0.3.1 roadmap — ex-post MDL
  reconstruction or an upstream bookkeeping PR to moneca is the open
  choice.

## Dependencies

* No new Imports or Suggests.

---

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
