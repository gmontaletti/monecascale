# graph-tool backend for moneca_sbm()
# ===================================
#
# Wraps Peixoto's graph-tool (Python) via reticulate. Runs hierarchical
# DC-SBM inference by minimum description length. Suitable for graphs with
# 10^5 to 10^7 nodes; requires a conda environment with graph-tool installed
# (graph-tool is not available on PyPI).

# 1. Installation helper ------------------------------------------------------

#' Install the graph-tool Python Backend for \code{moneca_sbm()}
#'
#' Convenience wrapper around \code{reticulate::conda_install()} that creates
#' (or updates) a conda environment with \code{graph-tool} installed from
#' conda-forge. \code{graph-tool} is not distributed on PyPI, so conda is
#' required.
#'
#' After installation, point reticulate at the environment before calling
#' \code{moneca_sbm(backend = "graphtool")}:
#' \preformatted{
#' reticulate::use_condaenv("moneca-sbm", required = TRUE)
#' }
#'
#' @param envname Name of the conda environment to create or update.
#'   Defaults to \code{"moneca-sbm"}.
#' @param method Installation method passed to \code{reticulate}. Defaults
#'   to \code{"conda"}; alternatives are rarely useful here.
#'
#' @return Invisibly returns the environment name.
#'
#' @seealso \code{\link{moneca_sbm}}
#'
#' @examples
#' \dontrun{
#' moneca_sbm_install_graphtool()
#' reticulate::use_condaenv("moneca-sbm", required = TRUE)
#' }
#'
#' @export
moneca_sbm_install_graphtool <- function(
  envname = "moneca-sbm",
  method = "conda"
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Install reticulate first: install.packages('reticulate').",
      call. = FALSE
    )
  }
  reticulate::conda_install(
    envname = envname,
    packages = "graph-tool",
    channel = "conda-forge",
    method = method
  )
  message(
    "graph-tool installed into conda env '",
    envname,
    "'. ",
    "Activate in R with: reticulate::use_condaenv('",
    envname,
    "', required = TRUE)"
  )
  invisible(envname)
}

# 2. Fit function -------------------------------------------------------------

.sbm_fit_graphtool <- function(
  core_mx,
  deg_corr = TRUE,
  edge_model = "poisson",
  seed = NULL,
  n_init = 5L,
  node_names = NULL,
  verbose = FALSE,
  ...
) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "moneca_sbm(backend = 'graphtool') requires reticulate. ",
      "Install with install.packages('reticulate').",
      call. = FALSE
    )
  }

  gt <- reticulate::import("graph_tool.all", convert = FALSE)

  if (!is.null(seed)) {
    gt$seed_rng(as.integer(seed))
    gt$openmp_set_num_threads(1L)
  }

  # 2.1 Build edge list from the (possibly sparse) core matrix -------------
  if (inherits(core_mx, "sparseMatrix")) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop(
        "Sparse input requires the 'Matrix' package.",
        call. = FALSE
      )
    }
    coo <- as(core_mx, "TsparseMatrix")
    ei <- coo@i
    ej <- coo@j
    ew <- coo@x
  } else {
    nz <- which(core_mx != 0, arr.ind = TRUE)
    ei <- nz[, 1] - 1L
    ej <- nz[, 2] - 1L
    ew <- core_mx[nz]
  }

  # Keep only strictly positive weights
  keep <- is.finite(ew) & ew > 0
  ei <- ei[keep]
  ej <- ej[keep]
  ew <- ew[keep]

  # Poisson SBM needs integer counts
  if (edge_model == "poisson") {
    ew <- as.integer(round(ew))
  }

  # 2.2 Construct the graph-tool Graph --------------------------------------
  is_directed <- !isSymmetric(as(core_mx, "matrix"))
  g <- gt$Graph(directed = is_directed)
  g$add_vertex(as.integer(nrow(core_mx)))

  edges_py <- reticulate::r_to_py(cbind(as.integer(ei), as.integer(ej)))
  g$add_edge_list(edges_py)

  # Attach edge weights as an edge property map
  w_prop_type <- if (edge_model == "poisson") "int" else "double"
  w <- g$new_ep(w_prop_type)
  weights_py <- reticulate::r_to_py(ew)
  w$a <- weights_py
  g$ep[["weight"]] <- w

  # 2.3 Inference -----------------------------------------------------------
  rec_type <- switch(
    edge_model,
    "poisson" = "discrete-poisson",
    "bernoulli" = NULL
  )
  state_args <- if (!is.null(rec_type)) {
    reticulate::dict(
      recs = list(w),
      rec_types = list(rec_type),
      deg_corr = deg_corr
    )
  } else {
    reticulate::dict(deg_corr = deg_corr)
  }

  fits <- vector("list", n_init)
  entropies <- numeric(n_init)
  for (i in seq_len(n_init)) {
    st <- gt$minimize_nested_blockmodel_dl(
      g,
      state_args = state_args
    )
    fits[[i]] <- st
    entropies[i] <- as.numeric(reticulate::py_to_r(st$entropy()))
  }
  best <- fits[[which.min(entropies)]]

  # Stabilize label ordering
  tryCatch(
    gt$inference$order_nested_partition_labels(best),
    error = function(e) invisible(NULL)
  )

  # 2.4 Extract per-level assignments --------------------------------------
  bs_py <- best$get_bs()
  bs <- reticulate::py_to_r(bs_py)
  memberships <- lapply(bs, function(arr) as.integer(as.numeric(arr)) + 1L)
  n_blocks <- vapply(
    memberships,
    function(v) length(unique(v)),
    integer(1)
  )

  # graph-tool levels: 0 = finest, going coarser. We keep only non-trivial
  # levels with n_blocks >= 2 and n_blocks < n_nodes.
  n_nodes <- nrow(core_mx)
  keep_lvl <- which(n_blocks >= 2L & n_blocks <= n_nodes)
  memberships <- memberships[keep_lvl]
  n_blocks <- n_blocks[keep_lvl]

  # MDL: report the best model's total entropy, repeated per kept level
  mdl_total <- as.numeric(reticulate::py_to_r(best$entropy()))
  mdl_vec <- rep(mdl_total, length(memberships))

  list(
    memberships = memberships,
    n_blocks = n_blocks,
    mdl = mdl_vec,
    backend = "graphtool",
    node_names = node_names,
    n_init = n_init
  )
}
