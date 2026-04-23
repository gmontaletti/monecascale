# greed backend for moneca_sbm()
# ==============================
#
# Wraps greed::greed() with DcSbm()/Sbm() model and exposes its hierarchical
# dendrogram as per-level partitions consumable by the moneca adapter.
#
# greed fits a flat DC-SBM at the ICL-optimal K, then builds a full nested
# partition tree from K down to 1 via hierarchical merges. We expose levels
# from K_final down to 2 (the "2 segments" level is informative; 1 is trivial).

# 1. Main fit function --------------------------------------------------------

.sbm_fit_greed <- function(
  core_mx,
  deg_corr = TRUE,
  edge_model = "poisson",
  seed = NULL,
  n_init = 1L,
  max_K = NULL,
  node_names = NULL,
  verbose = FALSE,
  ...
) {
  if (!requireNamespace("greed", quietly = TRUE)) {
    stop(
      "moneca_sbm(backend = 'greed') requires the 'greed' package. ",
      "Install it with install.packages('greed').",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  model <- if (deg_corr) greed::DcSbm() else greed::Sbm()

  n_core <- nrow(core_mx)
  k_default <- max(5L, as.integer(n_core / 3L))
  k_start <- max_K %||% k_default
  # greed requires K < n_nodes; clamp and cap at a practical ceiling
  k_start <- min(k_start, n_core - 1L, 50L)
  k_start <- max(k_start, 2L)

  # Coerce to integer counts for Poisson-like DC-SBM
  if (edge_model == "poisson") {
    if (inherits(core_mx, "sparseMatrix")) {
      core_mx@x <- as.numeric(round(core_mx@x))
    } else {
      core_mx[] <- round(core_mx)
    }
  }

  # Best-of-n_init: run greed multiple times, keep the highest ICL
  fits <- vector("list", n_init)
  icls <- numeric(n_init)
  for (i in seq_len(n_init)) {
    sol_i <- suppressWarnings(suppressMessages(greed::greed(
      core_mx,
      model = model,
      K = k_start,
      alg = greed::Hybrid(),
      verbose = verbose
    )))
    fits[[i]] <- sol_i
    icls[i] <- tryCatch(
      suppressWarnings(greed::ICL(sol_i)),
      error = function(e) NA_real_
    )
  }
  best_idx <- which.max(icls)
  sol <- fits[[best_idx]]

  K_final <- greed::K(sol)
  if (K_final < 2L) {
    # Everything merged; nothing to expose except the trivial atomic level
    return(list(
      memberships = list(),
      n_blocks = integer(0),
      mdl = numeric(0),
      backend = "greed",
      node_names = node_names,
      n_init = n_init
    ))
  }

  # Collect per-level partitions: K_final, K_final-1, ..., 2
  # Skip any level where K == n_core (all singletons; redundant with the
  # atomic level 1 produced by the adapter).
  k_seq <- seq.int(from = K_final, to = 2L, by = -1L)
  k_seq <- k_seq[k_seq < nrow(core_mx)]
  if (length(k_seq) == 0L) {
    return(list(
      memberships = list(),
      n_blocks = integer(0),
      mdl = numeric(0),
      backend = "greed",
      node_names = node_names,
      n_init = n_init
    ))
  }
  memberships <- vector("list", length(k_seq))
  mdl_vec <- numeric(length(k_seq))
  for (j in seq_along(k_seq)) {
    k <- k_seq[j]
    cut_sol <- tryCatch(
      suppressWarnings(greed::cut(sol, k)),
      error = function(e) NULL
    )
    if (is.null(cut_sol)) {
      memberships[[j]] <- rep(1L, nrow(core_mx))
      mdl_vec[j] <- NA_real_
      next
    }
    memberships[[j]] <- as.integer(greed::clustering(cut_sol))
    mdl_vec[j] <- tryCatch(
      suppressWarnings(-greed::ICL(cut_sol)),
      error = function(e) NA_real_
    )
  }

  # Drop any level that produced no non-singleton clique (e.g. degenerate
  # memberships where every block has size 1 even though nominal K > 1).
  keep <- vapply(
    memberships,
    function(m) any(table(m) >= 2L),
    logical(1)
  )
  memberships <- memberships[keep]
  k_seq <- k_seq[keep]
  mdl_vec <- mdl_vec[keep]

  list(
    memberships = memberships,
    n_blocks = as.integer(k_seq),
    mdl = mdl_vec,
    backend = "greed",
    node_names = node_names,
    n_init = n_init
  )
}
