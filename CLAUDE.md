# CLAUDE.md — monecascale

Guidance for Claude Code (claude.ai/code) when working in this repository.

## About monecascale

`monecascale` is a sibling package to [`moneca`](https://github.com/gmontaletti/MONECA)
that hosts **scalable clustering backends** for mobility network analysis. It extends moneca beyond classification-scale matrices (dozens to a few hundred categories) toward settings such as person × employer mobility (10^6 × 10^5).

Every exported entry point returns a `moneca`-class S3 object, so moneca's analysis and plotting stack consumes the output unchanged.

## Relationship to moneca

- **moneca**: reference algorithm (`moneca()`, `moneca_fast()`), weight matrix (`weight.matrix()`), segment metadata (`moneca_segments()`), the full `segment.*()` and `plot_moneca_*()` consumer API.
- **monecascale**: `Imports: moneca`. Adds alternative backends whose output is consumed by moneca unchanged.
- Local path of the sister package: `/Users/giampaolomontaletti/Documents/funzioni/MONECA/`.

## Shipped entry points

- `moneca_sbm()` — hierarchical degree-corrected stochastic block model. Backends:
  - `"greed"` (default, R-native, sparse `dgCMatrix`, scales to ~10^4–10^5 nodes)
  - `"graphtool"` (Python via reticulate; conda-forge install; scales to 10^6+ nodes)
- `moneca_sbm_install_graphtool()` — conda-forge install helper for the graphtool backend.
- `moneca_bipartite()` — hierarchical bipartite co-clustering via `greed::DcLbm`
  on rectangular RR. Returns `moneca_bipartite` wrapping two moneca-class views
  (`$rows`, `$cols`). Shipped in 0.2.0.
- `auto_segment_levels()` — post-hoc hierarchy-level selection via MDL elbow
  or mutual information plateau; also available as `segment.levels = "auto"`
  on `moneca_sbm()` and `moneca_bipartite()`. Shipped in 0.3.0.
- `moneca_flow()` — flow-based hierarchical clustering via Infomap (Map Equation).
  Backend via `igraph::cluster_infomap()` (R-native, scales to ~10^5 nodes). Recursive-flat
  wrapper for 2–3 level hierarchy. Shipped in 0.4.0.

## Planned directions (scaling roadmap)

Landed here unless otherwise noted. See `~/.claude/plans/monecascale-scaling-roadmap.md` and `~/.claude/projects/.../memory/scaling_directions.md`.

- ~~**D1** `moneca_bipartite()`~~ — **SHIPPED 0.2.0**.
- ~~**D6** `auto_segment_levels()`~~ — **SHIPPED 0.3.0** for SBM and bipartite. Fast-backend support added in 0.4.1.
- ~~**D3** `moneca_flow()`~~ — **SHIPPED 0.4.0** Infomap / Map Equation.
- **D4** `moneca_localclique()` — quasi-cliques, k-clique percolation, local Personalized PageRank expansion.
- **D5** `moneca_nmf()` — Poisson NMF.
- **D7** `moneca_duckdb()` / `moneca_arrow()` — out-of-core / streaming backend.

## Package Development Commands

```r
# Interactive
devtools::load_all()

# Document / check / test
devtools::document()
devtools::test()
devtools::check()

# Build / install
devtools::build()
devtools::install()

# Vignettes
devtools::build_vignettes()

# Coverage
covr::package_coverage()
```

## File organization

- `R/monecascale-package.R` — package-level roxygen + `"_PACKAGE"`.
- `R/moneca_sbm.R` — public `moneca_sbm()` entry + backend dispatcher.
- `R/sbm_backend_greed.R` — R-native greed backend.
- `R/sbm_backend_graphtool.R` — Python/graph-tool backend + install helper.
- `R/sbm_utils.R` — shared preprocessing (margins, symmetry, diagonal zeroing) and the `.sbm_fit_to_moneca()` adapter.
- `R/sparse_internal.R` — local copies of two internal moneca helpers (`segment_matrix_sparse`, `sparse_pmin_symmetric`). Avoid calling `moneca:::` from this package.
- `R/moneca_bipartite.R` — public `moneca_bipartite()` entry + backend dispatcher.
- `R/bipartite_backend_greed.R` — greed DcLbm backend.
- `R/bipartite_utils.R` — rectangular RR, one-mode projection, adapter.
- `R/bipartite_methods.R` — `print` / `summary` / `format` for `moneca_bipartite`.
- `R/auto_segment_levels.R` — public `auto_segment_levels()` entry + dispatcher + trim helper.
- `R/auto_level_criteria.R` — MDL (kneedle / max_second_diff) and MI-plateau criteria.
- `R/auto_level_methods.R` — `print` / `format` for `auto_segment_levels`.
- `R/moneca_flow.R` — public `moneca_flow()` entry + backend dispatcher.
- `R/flow_backend_igraph.R` — igraph `cluster_infomap()` backend.
- `R/flow_utils.R` — recursive-flat hierarchy synthesis and adapter.
- `tests/testthat/helper-test-data.R` — shared generators wrapping `moneca::generate_mobility_data()`.
- `vignettes/monecascale-sbm.Rmd` — RR ≡ DC-SBM bridge walkthrough and parity experiment.
- `vignettes/monecascale-bipartite.Rmd` — bipartite walkthrough + recovery experiment.
- `vignettes/monecascale-auto-level.Rmd` — auto-level walkthrough.
- `vignettes/monecascale-flow.Rmd` — Map Equation, flat-vs-hierarchical trade-off, comparison vs SBM.

## Coding conventions

- **2-space indentation**.
- Section comments in R scripts: `# 1. name of section -----`, no `####` dividers.
- **Explicit `moneca::` qualification** whenever calling moneca APIs — makes dependency crossings visible.
- **No `moneca:::internal`** references. Copy internal helpers locally with a comment pointing at the upstream file.
- Variable naming: mixed Italian / English allowed to match the moneca domain (`retribuzione`, `arco`, `chunk_id`, `overlap_flag`).
- Italian prose must use properly accented characters (`è`, `é`, `à`, `ò`, `ù`, `ì`, `perché`, `più`). Primary docs are in English; Italian reserved for translations.
- Neutral technical tone. Avoid emphatic phrasing.

## Output contract (protected behavior)

Every entry point must return a list of class `"moneca"` with at minimum:

- `segment.list` — list of lists. `segment.list[[1]]` atomic (`as.list(1:n_core)`). Levels ≥ 2 are lists of integer-vector cliques with singletons pruned (length ≥ 2).
- `mat.list` — list of matrices. `mat.list[[1]]` is the input with margins; subsequent levels are aggregated via the matching membership vector.
- `small.cell.reduction`, `margins_added`, `density_reduction`, `symmetric_method` — scalars / flags mirroring `moneca::moneca_fast()`.
- `segment_metadata` — the result of `moneca::moneca_segments(out)` called on the otherwise-complete object.

Direction-specific diagnostics (e.g. `$sbm_diagnostics`) go in an extra slot, never inside the core slots above.

## Testing

- `testthat (>= 3.0.0)`. Edition 3.
- Optional-backend tests must `skip_if_not_installed("greed")` or `skip_if_not(reticulate::py_module_available("graph_tool"))`.
- New directions should add: a core-behavior test file, a parity / sanity test file, and a skippable backend test file when an optional dependency is involved.
- Keep test fixtures small (n ≤ 20 categories for unit tests; n = 20–50 for parity fixtures).

## Adding a new direction (D1 / D3 / D4 / D5 / D7)

1. Create `R/moneca_<direction>.R` with the public entry and the backend dispatcher (mirror `moneca_sbm.R`).
2. If needed, add `R/<direction>_backend_*.R` per alternative backend.
3. Reuse `R/sbm_utils.R` helpers where applicable (especially `.ensure_margins_sbm`, `.preprocess_for_sbm`, `.sbm_fit_to_moneca`). If the adapter needs to differ, extract a shared core helper rather than duplicating.
4. Add `tests/testthat/test-<direction>.R` — core behavior + parity vs `moneca::moneca_fast()` where sensible.
5. Add `vignettes/monecascale-<direction>.Rmd`.
6. Register new exports in `_pkgdown.yml` Core Clustering Backends section; add vignette to the articles section.
7. `NEWS.md` entry with a minor version bump (`0.1.0 → 0.2.0` etc.).

## Dependency policy

- `Imports: moneca, Matrix, stats` — hard deps.
- `Suggests`: backend-specific deps (`greed`, `reticulate`), plus test/doc tooling.
- New directions should continue to gate optional backends via `Suggests` + `requireNamespace`. Do not move backend libraries to `Imports` without a strong case.

## Dependency Management

`renv` is initialized for this project. Use:
```r
renv::status()     # check state
renv::snapshot()   # after adding a dep
renv::restore()    # to rehydrate
```

## Repository cleanup

Temporary files, drafts, development artifacts → `../reference/monecascale/` (sibling to the package dir, not inside it). `.gitignore` already excludes `renv/library/`, `renv/staging/`, `vignettes/*.R`, `vignettes/*.html`, etc.

## Agent / skill phase assignment

| Phase          | Agents                                                                 | Skills                                 |
|----------------|------------------------------------------------------------------------|----------------------------------------|
| Design         | `r-project-orchestrator`                                               | `r-package-setup`                      |
| Implementation | `r-function-developer`, `r-documentation-expert`                       | `data-table-skill`, `cvd-accessible-visualization` |
| Verification   | `r-unit-tester`, `r-debug-specialist`, `r-performance-optimizer`       | `r-profiling`                          |
| Release        | `italian-doc-translator`                                               | `commit`, `r-package-maintenance`, `pkgdown-skill-package` |
| Maintenance    | `r-project-maintainer`, `r-project-orchestrator`                       | `r-profiling`                          |

## Notes

- Author / maintainer: Giampaolo Montaletti, `giampaolo.montaletti@gmail.com`, ORCID https://orcid.org/0009-0002-5327-1122.
- When bumping the package version, also update the citation block in `README.md`.
- The `.claude` directory at the repo root (if any) should not be moved during cleanups.
