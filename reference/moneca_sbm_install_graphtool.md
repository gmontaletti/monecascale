# Install the graph-tool Python Backend for `moneca_sbm()`

Convenience wrapper around
[`reticulate::conda_install()`](https://rstudio.github.io/reticulate/reference/conda-tools.html)
that creates (or updates) a conda environment with `graph-tool`
installed from conda-forge. `graph-tool` is not distributed on PyPI, so
conda is required.

## Usage

``` r
moneca_sbm_install_graphtool(envname = "moneca-sbm", method = "conda")
```

## Arguments

- envname:

  Name of the conda environment to create or update. Defaults to
  `"moneca-sbm"`.

- method:

  Installation method passed to `reticulate`. Defaults to `"conda"`;
  alternatives are rarely useful here.

## Value

Invisibly returns the environment name.

## Details

After installation, point reticulate at the environment before calling
`moneca_sbm(backend = "graphtool")`:

    reticulate::use_condaenv("moneca-sbm", required = TRUE)

## See also

[`moneca_sbm`](https://gmontaletti.github.io/monecascale/reference/moneca_sbm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
moneca_sbm_install_graphtool()
reticulate::use_condaenv("moneca-sbm", required = TRUE)
} # }
```
