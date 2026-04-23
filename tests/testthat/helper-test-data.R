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
