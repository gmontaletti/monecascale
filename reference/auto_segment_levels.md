# Auto-select a Hierarchy Level for a \`moneca\` Fit

Picks a preferred level in a hierarchical moneca clustering based on one
of two criteria: elbow on the Minimum Description Length (MDL) trace, or
plateau detection on level-to-level mutual information. Dispatches on
\[moneca_sbm()\] output directly, and recurses into the \`\$rows\` and
\`\$cols\` children of a \[moneca_bipartite()\] object.

## Usage

``` r
auto_segment_levels(
  obj,
  method = c("mdl", "mi_plateau"),
  return = c("full", "level"),
  plateau_tol = 0.05,
  mdl_elbow = c("kneedle", "max_second_diff"),
  verbose = FALSE,
  ...
)
```

## Arguments

- obj:

  A \`moneca\` object produced by \[moneca_sbm()\], a
  \`moneca_bipartite\` object produced by \[moneca_bipartite()\], or a
  \`\$rows\` / \`\$cols\` child of such an object.

- method:

  One of \`"mdl"\` (default) or \`"mi_plateau"\`.

- return:

  One of \`"full"\` (default, returns an \`auto_segment_levels\` object)
  or \`"level"\` (returns a bare integer scalar).

- plateau_tol:

  Positive numeric in \`(0, 1)\`. Relative gain threshold under which
  the MI curve is considered flat. Only used when \`method =
  "mi_plateau"\`.

- mdl_elbow:

  One of \`"kneedle"\` (default) or \`"max_second_diff"\`. Only used
  when \`method = "mdl"\`.

- verbose:

  Logical. Emit per-level diagnostic messages.

- ...:

  Reserved for forward compatibility.

## Value

When \`obj\` inherits from \`"moneca_bipartite"\`, a list with named
components \`rows\` and \`cols\`, each an \`auto_segment_levels\` object
(or integer scalar if \`return = "level"\`). Otherwise, when \`return =
"full"\`, an object of class \`auto_segment_levels\` with components: \*
\`level\` - integer scalar, the selected \`segment.list\` index. \*
\`method\` - the criterion used. \* \`diagnostics\` - a data.frame with
columns \`level\`, \`n_blocks\`, \`mdl\`, \`mi_to_next\`, \`score\`. \*
\`backend\` - one of \`"sbm"\`, \`"bipartite_rows"\`,
\`"bipartite_cols"\`. \* \`call\` - the matched call.

When \`return = "level"\`, a single integer scalar.

## Details

The \`"mdl"\` criterion fetches the per-level MDL series from
\`obj\$sbm_diagnostics\$mdl_per_level\` (for SBM backends) or from the
slim per-side slice
\`obj\$bipartite_diagnostics_side\$joint_mdl_per_level\` (for a
\`\$rows\` or \`\$cols\` child of a \`moneca_bipartite\` object). The
series is reordered to coarse-to-fine (ascending block count) and an
elbow is selected by either the Kneedle heuristic or the discrete
second-difference peak.

The \`"mi_plateau"\` criterion computes level-to-level sample mutual
information from \`obj\$segment.list\` directly (no reliance on backend
diagnostics). It picks the smallest level at which the gain \\I_l -
I\_{l-1}\\ falls below \`plateau_tol \* max(I)\`; if no plateau is
detected, it returns the finest level.

For \`moneca::moneca_fast()\` output the function errors out: neither an
MDL series nor a principled MI plateau is available without a
per-iteration trace, which upstream \`moneca_fast()\` does not emit.
Support is planned for monecascale 0.3.1.

## See also

\[moneca_sbm()\], \[moneca_bipartite()\].

## Examples

``` r
if (requireNamespace("greed", quietly = TRUE)) {
  mob <- moneca::generate_mobility_data(n_classes = 20, seed = 1L)
  fit <- moneca_sbm(mob, backend = "greed", seed = 1L)
  auto_segment_levels(fit, method = "mdl")
}
#> <auto_segment_levels> method=mdl backend=sbm picked level=11 of 16
#>  level n_blocks       mdl mi_to_next     score
#>      2       17  923.7724         NA  923.7724
#>      3       16  946.9383         NA  946.9383
#>      4       15  987.3845         NA  987.3845
#>      5       14 1038.0860         NA 1038.0860
#>      6       13 1103.8208         NA 1103.8208
#>      7       12 1166.7615         NA 1166.7615
#>      8       11 1242.3147         NA 1242.3147
#>      9       10 1307.2152         NA 1307.2152
#>     10        9 1446.5529         NA 1446.5529
#>     11        8 1591.6417         NA 1591.6417
#>     12        7 1765.0089         NA 1765.0089
#>     13        6 1954.2520         NA 1954.2520
#>     14        5 2161.0624         NA 2161.0624
#>     15        4 2494.8719         NA 2494.8719
#>     16        3 2835.2535         NA 2835.2535
#>     17        2 3424.7448         NA 3424.7448
```
