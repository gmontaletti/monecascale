# Shared utilities for moneca_flow()
# ==================================
#
# Recursive Infomap wrapper, direct map-equation codelength computation,
# and the moneca-class adapter. All helpers in this file are internal.

# 1. Recursive flat Infomap ---------------------------------------------------
#
# For a given `parent_membership` vector on `g`, induce each parent
# module's subgraph, run flat Infomap on it, and stitch the sub-module
# labels into a single level-wide membership vector with a running
# integer counter. Parent modules of size < 2 keep a singleton label.

.recursive_infomap <- function(
  g,
  parent_membership,
  nb_trials = 10L,
  verbose = FALSE,
  level_label = 2L
) {
  n <- igraph::vcount(g)
  parent_membership <- as.integer(parent_membership)
  stopifnot(length(parent_membership) == n)

  child_membership <- integer(n)
  next_id <- 0L

  parent_ids <- sort(unique(parent_membership))
  for (k in parent_ids) {
    idx <- which(parent_membership == k)
    if (length(idx) == 0L) {
      next
    }

    if (length(idx) < 2L) {
      # Singleton parent: assign its node a single new id
      next_id <- next_id + 1L
      child_membership[idx] <- next_id
      next
    }

    g_k <- igraph::induced_subgraph(g, idx)
    if (igraph::ecount(g_k) == 0L) {
      # No internal edges: every node becomes its own sub-module
      sub_m <- seq_len(length(idx))
    } else {
      sol <- .run_infomap(g_k, nb_trials = nb_trials, verbose = FALSE)
      sub_m <- as.integer(sol$membership)
    }

    sub_ids <- sort(unique(sub_m))
    id_map <- integer(max(sub_ids))
    for (s in sub_ids) {
      next_id <- next_id + 1L
      id_map[s] <- next_id
    }
    child_membership[idx] <- id_map[sub_m]
  }

  if (verbose) {
    message(sprintf(
      "[.recursive_infomap] level %d: built %d sub-modules from %d parents",
      level_label,
      next_id,
      length(parent_ids)
    ))
  }

  child_membership
}

# 2. Map-equation codelength for an arbitrary partition ----------------------
#
# Direct computation of the two-level map equation (Rosvall & Bergstrom
# 2008):
#
#   L(M) = q_exit * H(Q) + sum_k p_k_total * H(P_k)
#
# where
#   * p_i is the stationary distribution on nodes (PageRank for
#     robustness against dangling nodes in directed graphs)
#   * q_k is the per-module exit probability (total out-weight of
#     edges leaving module k, normalised by total flow)
#   * q_exit = sum_k q_k is the total module-switching rate
#   * H(Q) is the entropy of the module-index codebook with
#     probabilities proportional to {q_k}
#   * P_k is the within-module codebook for module k: node
#     probabilities proportional to {p_i : i in k} together with the
#     per-module exit probability q_k; p_k_total is its total use rate.
#
# Returns a single non-negative scalar (nats-based, log base 2).

.codelength_of_partition <- function(g, membership, directed = TRUE) {
  membership <- as.integer(membership)
  n <- igraph::vcount(g)
  stopifnot(length(membership) == n)

  if (n == 0L) {
    return(0)
  }
  if (max(membership) == 1L) {
    # Single module: no module-switching. Codelength collapses to the
    # entropy of the stationary distribution.
    p <- .stationary_distribution(g, directed = directed)
    return(.entropy2(p))
  }

  # 2.1 Stationary distribution on nodes -----------------------------------
  p <- .stationary_distribution(g, directed = directed)

  # 2.2 Per-module quantities via edge list --------------------------------
  w <- igraph::E(g)$weight
  if (is.null(w)) {
    w <- rep(1, igraph::ecount(g))
  }
  total_w <- sum(w)

  if (total_w <= 0) {
    # Degenerate graph with no flow: codelength undefined, return 0.
    return(0)
  }

  el <- igraph::as_edgelist(g, names = FALSE)
  src_mod <- membership[el[, 1L]]
  dst_mod <- membership[el[, 2L]]

  K <- max(membership)

  # q_k: probability that the walker leaves module k (out-exits / total flow)
  # Undirected graphs also use a directional edgelist-based count; igraph
  # stores each undirected edge once so this is consistent.
  exit_w <- numeric(K)
  cross <- src_mod != dst_mod
  if (any(cross)) {
    from_mods <- src_mod[cross]
    exit_w_vec <- tapply(w[cross], from_mods, sum)
    exit_w[as.integer(names(exit_w_vec))] <- as.numeric(exit_w_vec)
  }
  q_k <- exit_w / total_w

  # p_k_sum: sum of node visit rates in module k
  p_k_sum <- tapply(p, membership, sum)
  p_k_full <- numeric(K)
  p_k_full[as.integer(names(p_k_sum))] <- as.numeric(p_k_sum)

  # 2.3 Module-switching codebook Q ----------------------------------------
  q_exit <- sum(q_k)
  H_Q <- if (q_exit > 0) .entropy2(q_k / q_exit) else 0

  # 2.4 Within-module codebooks P_k ----------------------------------------
  # Each P_k uses probabilities {q_k} (for exiting) plus {p_i : i in k}.
  # Its total use rate is q_k + p_k_full[k].
  L_within <- 0
  for (k in seq_len(K)) {
    members <- which(membership == k)
    if (length(members) == 0L) {
      next
    }
    p_members <- p[members]
    use_total <- q_k[k] + p_k_full[k]
    if (use_total <= 0) {
      next
    }
    probs <- c(q_k[k], p_members) / use_total
    H_k <- .entropy2(probs)
    L_within <- L_within + use_total * H_k
  }

  L <- q_exit * H_Q + L_within
  as.numeric(L)
}

# 3. PageRank-based stationary distribution ----------------------------------

.stationary_distribution <- function(g, directed = TRUE) {
  if (igraph::vcount(g) == 0L) {
    return(numeric(0))
  }
  w <- igraph::E(g)$weight
  pr <- tryCatch(
    igraph::page_rank(
      g,
      damping = 0.85,
      weights = w,
      directed = directed
    )$vector,
    error = function(e) NULL
  )
  if (is.null(pr) || !all(is.finite(pr)) || sum(pr) <= 0) {
    # Fallback: uniform distribution
    n <- igraph::vcount(g)
    return(rep(1 / n, n))
  }
  pr / sum(pr)
}

# 4. Shannon entropy in bits, tolerating zeros -------------------------------

.entropy2 <- function(p) {
  p <- as.numeric(p)
  p <- p[p > 0]
  if (length(p) == 0L) {
    return(0)
  }
  -sum(p * log2(p))
}

# 5. Main adapter: flow_fit + mx_full => moneca S3 object --------------------

.flow_fit_to_moneca <- function(
  fit,
  mx_full,
  margins_added,
  small.cell.reduction,
  symmetric_method,
  directed
) {
  n_core <- nrow(mx_full) - 1L

  memberships <- fit$memberships
  n_levels <- length(memberships)

  # 5.1 segment.list: atomic level + one clique list per hierarchy level ---
  segment_list <- vector("list", n_levels + 1L)
  segment_list[[1]] <- as.list(seq_len(n_core))
  for (lvl in seq_len(n_levels)) {
    segment_list[[lvl + 1L]] <- .memberships_to_cliques(memberships[[lvl]])
  }

  # 5.2 mat.list: full matrix + aggregated per level -----------------------
  mat_list <- vector("list", n_levels + 1L)
  mat_list[[1]] <- mx_full
  for (lvl in seq_len(n_levels)) {
    mat_list[[lvl + 1L]] <- .aggregate_matrix_by_membership(
      mat_list[[1]],
      as.integer(memberships[[lvl]])
    )
  }

  # 5.3 Diagnostics --------------------------------------------------------
  flow_diagnostics <- list(
    backend = fit$backend %||% "igraph",
    model = "infomap",
    depth = fit$depth %||% n_levels,
    n_blocks_per_level = fit$n_blocks_per_level,
    codelength_per_level = fit$codelength_per_level,
    # Alias so any tooling keyed on `$mdl_per_level` works unchanged.
    mdl_per_level = fit$codelength_per_level,
    n_init = fit$n_init %||% 1L,
    directed = directed,
    node_names = fit$node_names %||% rownames(mx_full)[1:n_core]
  )

  # 5.4 Assemble moneca-class object --------------------------------------
  out <- list(
    segment.list = segment_list,
    mat.list = mat_list,
    small.cell.reduction = small.cell.reduction,
    margins_added = margins_added,
    density_reduction = NULL,
    symmetric_method = symmetric_method,
    flow_diagnostics = flow_diagnostics
  )

  class(out) <- "moneca"
  out$segment_metadata <- moneca::moneca_segments(out)
  out
}
