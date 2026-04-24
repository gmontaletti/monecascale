# igraph backend for moneca_flow()
# ================================
#
# Wraps igraph::cluster_infomap() with a recursive wrapper that
# synthesises a 2- or 3-level hierarchy from the flat partition.
#
# The sub-level Infomap runs only serve to produce membership vectors.
# Codelengths for levels >= 2 are computed directly on the full graph
# via .codelength_of_partition() so the per-level trace is comparable.

# 1. Main fit function --------------------------------------------------------

.flow_fit_igraph <- function(
  core_mx,
  depth = 2L,
  nb_trials = 10L,
  directed = TRUE,
  zero_diagonal = TRUE,
  seed = NULL,
  verbose = FALSE
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "moneca_flow(backend = 'igraph') requires the 'igraph' package. ",
      "Install with install.packages('igraph').",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_core <- nrow(core_mx)

  # 2. Build the top-level igraph object ---------------------------------
  g <- igraph::graph_from_adjacency_matrix(
    core_mx,
    mode = if (directed) "directed" else "undirected",
    weighted = TRUE,
    diag = !zero_diagonal
  )

  # Guard: isolated nodes (zero strength) cannot host random-walk flow
  # but igraph tolerates them. We pass the full graph and record them
  # as singleton modules via cluster_infomap()'s own handling.

  # 3. Level 1: flat Infomap on the full graph ---------------------------
  sol1 <- .run_infomap(g, nb_trials = nb_trials, verbose = verbose)
  m1 <- as.integer(sol1$membership)
  codelen1 <- as.numeric(sol1$codelength)

  memberships <- vector("list", depth)
  codelens <- numeric(depth)
  n_blocks <- integer(depth)

  memberships[[1]] <- m1
  codelens[1] <- codelen1
  n_blocks[1] <- max(m1)

  if (verbose) {
    message(sprintf(
      "[.flow_fit_igraph] level 1: K = %d, codelength = %.4f",
      n_blocks[1],
      codelens[1]
    ))
  }

  # 4. Level 2: recurse within each level-1 module -----------------------
  if (depth >= 2L) {
    m2 <- .recursive_infomap(
      g = g,
      parent_membership = memberships[[1]],
      nb_trials = nb_trials,
      verbose = verbose,
      level_label = 2L
    )
    memberships[[2]] <- m2
    codelens[2] <- .codelength_of_partition(
      g,
      membership = m2,
      directed = directed
    )
    n_blocks[2] <- max(m2)

    if (verbose) {
      message(sprintf(
        "[.flow_fit_igraph] level 2: K = %d, codelength = %.4f",
        n_blocks[2],
        codelens[2]
      ))
    }
  }

  # 5. Level 3: recurse within each level-2 module -----------------------
  if (depth >= 3L) {
    m3 <- .recursive_infomap(
      g = g,
      parent_membership = memberships[[2]],
      nb_trials = nb_trials,
      verbose = verbose,
      level_label = 3L
    )
    memberships[[3]] <- m3
    codelens[3] <- .codelength_of_partition(
      g,
      membership = m3,
      directed = directed
    )
    n_blocks[3] <- max(m3)

    if (verbose) {
      message(sprintf(
        "[.flow_fit_igraph] level 3: K = %d, codelength = %.4f",
        n_blocks[3],
        codelens[3]
      ))
    }
  }

  list(
    memberships = memberships,
    n_blocks_per_level = as.integer(n_blocks),
    codelength_per_level = as.numeric(codelens),
    backend = "igraph",
    depth = depth,
    n_init = nb_trials,
    directed = directed,
    node_names = rownames(core_mx)
  )
}

# 6. Thin wrapper around igraph::cluster_infomap() ---------------------------

.run_infomap <- function(g, nb_trials = 10L, verbose = FALSE) {
  # Some igraph versions warn on zero-weight edges or on graphs with
  # isolated nodes; swallow the noise and let the clustering speak.
  ew <- igraph::E(g)$weight
  suppressWarnings(suppressMessages(
    igraph::cluster_infomap(
      g,
      e.weights = ew,
      nb.trials = nb_trials,
      modularity = FALSE
    )
  ))
}
