# Shared test-data generators for monecascale
# ===========================================
#
# Thin wrappers around moneca::generate_mobility_data() with standard
# sizes, keeping tests self-contained and fixtures reproducible.

get_test_data <- function(size = "small", seed = 42) {
  sizes <- list(small = 5, medium = 10, large = 20)
  if (!size %in% names(sizes)) {
    stop("size must be 'small', 'medium', or 'large'")
  }
  moneca::generate_mobility_data(n_classes = sizes[[size]], seed = seed)
}

get_custom_test_data <- function(
  n_classes = 5,
  n_total = 1000,
  immobility_strength = 0.7,
  class_clustering = 0.7,
  noise_level = 0.1,
  seed = 123
) {
  moneca::generate_mobility_data(
    n_classes = n_classes,
    n_total = n_total,
    immobility_strength = immobility_strength,
    class_clustering = class_clustering,
    noise_level = noise_level,
    seed = seed
  )
}

# 1. Bipartite fixture for moneca_bipartite() tests --------------------------
#
# Rectangular Poisson draws from a block-constant intensity matrix with a
# checkerboard high/low pattern over latent (row_block, col_block) pairs.
# Returns a plain dense matrix with attributes `row_block` and `col_block`
# carrying the ground-truth assignments for parity tests.

get_bipartite_test_data <- function(
  n_rows = 30,
  n_cols = 20,
  k_row = 3,
  k_col = 2,
  seed = 2026
) {
  set.seed(seed)
  row_block <- as.integer(
    ((seq_len(n_rows) - 1L) %% k_row) + 1L
  )
  col_block <- as.integer(
    ((seq_len(n_cols) - 1L) %% k_col) + 1L
  )

  Omega <- matrix(0.3, nrow = k_row, ncol = k_col)
  for (r in seq_len(k_row)) {
    for (c in seq_len(k_col)) {
      if ((r + c) %% 2L == 0L) {
        Omega[r, c] <- 3.0
      }
    }
  }

  lambda <- Omega[row_block, col_block]
  mx <- matrix(
    stats::rpois(n_rows * n_cols, lambda = lambda),
    nrow = n_rows,
    ncol = n_cols
  )
  rownames(mx) <- paste0("R", seq_len(n_rows))
  colnames(mx) <- paste0("C", seq_len(n_cols))
  attr(mx, "row_block") <- row_block
  attr(mx, "col_block") <- col_block
  mx
}

# 2. Flow fixture for moneca_flow() tests ------------------------------------
#
# Square Poisson draws from a planted-community model. Nodes are assigned
# modules round-robin; within-module pairs draw Poisson(n * p_in),
# between-module pairs draw Poisson(n * p_out). Diagonal zeroed so the
# fixture carries no self-loops. Returns a dense matrix with the
# ground-truth module membership in `attr(, "module")`.

get_flow_test_data <- function(
  n = 40,
  n_modules = 4,
  p_in = 0.3,
  p_out = 0.02,
  seed = 2026
) {
  set.seed(seed)
  mod <- rep(seq_len(n_modules), length.out = n)
  lambda_in <- n * p_in
  lambda_out <- n * p_out

  mx <- matrix(0L, nrow = n, ncol = n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mx[i, j] <- stats::rpois(
        1,
        lambda = if (mod[i] == mod[j]) lambda_in else lambda_out
      )
    }
  }
  diag(mx) <- 0L

  rownames(mx) <- paste0("N", seq_len(n))
  colnames(mx) <- paste0("N", seq_len(n))
  attr(mx, "module") <- as.integer(mod)
  mx
}
