# Criterion implementations for auto_segment_levels()
# ===================================================
#
# Two pure-R criteria: MDL elbow detection and level-to-level mutual
# information plateau. Plus shared primitives (kneedle / second diff /
# membership reconstruction / sample MI). All internal helpers.

# 1. MDL elbow criterion ------------------------------------------------------

#' @keywords internal
#' @noRd
.criterion_mdl <- function(obj, backend, elbow, verbose) {
  src <- switch(
    backend,
    flow = list(
      mdl = obj$flow_diagnostics$codelength_per_level,
      n_blocks = obj$flow_diagnostics$n_blocks_per_level,
      score_col = "codelength"
    ),
    sbm = list(
      mdl = obj$sbm_diagnostics$mdl_per_level,
      n_blocks = obj$sbm_diagnostics$n_blocks_per_level,
      score_col = "mdl"
    ),
    bipartite_rows = list(
      mdl = obj$bipartite_diagnostics_side$joint_mdl_per_level,
      n_blocks = obj$bipartite_diagnostics_side$n_blocks_per_level,
      score_col = "mdl"
    ),
    bipartite_cols = list(
      mdl = obj$bipartite_diagnostics_side$joint_mdl_per_level,
      n_blocks = obj$bipartite_diagnostics_side$n_blocks_per_level,
      score_col = "mdl"
    ),
    fast = {
      scored <- .fast_compute_modularity(obj)
      list(
        mdl = -scored$modularity,
        n_blocks = scored$n_blocks,
        score_col = "modularity",
        raw_score = scored$modularity
      )
    },
    stop(
      "Unsupported backend for .criterion_mdl(): ",
      backend,
      call. = FALSE
    )
  )

  mdl <- as.numeric(src$mdl)
  n_blocks <- as.integer(src$n_blocks)
  score_col <- src$score_col

  # Trivial fit: no backend levels recovered. The only available level is
  # the atomic level 1 in segment.list. Return that rather than erroring.
  if (length(mdl) == 0L || length(n_blocks) == 0L) {
    diag_df_empty <- data.frame(
      level = integer(0),
      n_blocks = integer(0),
      mi_to_next = numeric(0),
      score = numeric(0),
      stringsAsFactors = FALSE
    )
    # Insert a typed empty score column under the backend-specific name.
    diag_df_empty[[score_col]] <- numeric(0)
    diag_df_empty <- diag_df_empty[, c(
      "level",
      "n_blocks",
      score_col,
      "mi_to_next",
      "score"
    )]
    return(list(
      level = 1L,
      method = "mdl",
      diagnostics = diag_df_empty
    ))
  }
  if (length(mdl) != length(n_blocks)) {
    stop(
      "mdl and n_blocks lengths disagree (",
      length(mdl),
      " vs ",
      length(n_blocks),
      ").",
      call. = FALSE
    )
  }

  # 2. Standardise ordering to coarse -> fine (ascending n_blocks) ------------
  ord <- order(n_blocks)
  mdl_ord <- mdl[ord]

  L <- length(mdl_ord)

  # 3. Elbow pick in ordered space --------------------------------------------
  idx_ord <- switch(
    elbow,
    kneedle = .kneedle(mdl_ord),
    max_second_diff = .max_second_diff(mdl_ord),
    stop("Unknown elbow mode: ", elbow, call. = FALSE)
  )

  # 4. Map back to the backend index and then to segment.list index ----------
  # segment.list layout: [[1]] atomic + [[2..L+1]] backend levels aligned to
  # the original (unsorted) per-level arrays.
  backend_idx <- ord[idx_ord]
  segment_level <- backend_idx + 1L

  if (verbose) {
    message(sprintf(
      "[.criterion_mdl] backend=%s L=%d elbow=%s picked backend idx=%d -> segment level=%d",
      backend,
      L,
      elbow,
      backend_idx,
      segment_level
    ))
  }

  diag_df <- data.frame(
    level = seq_len(L) + 1L,
    n_blocks = n_blocks,
    mi_to_next = NA_real_,
    score = mdl,
    stringsAsFactors = FALSE
  )
  # Insert the backend-specific score column ("mdl" or "codelength")
  # between n_blocks and mi_to_next so the layout is stable.
  diag_df[[score_col]] <- if (!is.null(src$raw_score)) src$raw_score else mdl
  diag_df <- diag_df[, c(
    "level",
    "n_blocks",
    score_col,
    "mi_to_next",
    "score"
  )]

  list(
    level = as.integer(segment_level),
    method = "mdl",
    diagnostics = diag_df
  )
}

# 5. Mutual-information plateau criterion ------------------------------------

#' @keywords internal
#' @noRd
.criterion_mi_plateau <- function(obj, plateau_tol, verbose) {
  L <- length(obj$segment.list)
  stopifnot(L >= 1L)

  n <- length(obj$segment.list[[1]])

  # 6. Membership vector per level ------------------------------------------
  memb <- vector("list", L)
  memb[[1]] <- seq_len(n)
  if (L >= 2L) {
    for (lvl in 2:L) {
      memb[[lvl]] <- .cliques_to_membership(obj$segment.list[[lvl]], n)
    }
  }

  # 7. Level-to-next MI: I_l = MI(m[l], m[l + 1]) for l in 1..L-1 -----------
  mi <- rep(NA_real_, L)
  if (L >= 2L) {
    for (lvl in seq_len(L - 1L)) {
      mi[lvl] <- .mutual_information(memb[[lvl]], memb[[lvl + 1L]])
    }
  }

  # 8. Plateau pick ----------------------------------------------------------
  picked <- L
  if (L >= 3L && any(is.finite(mi))) {
    I_max <- max(mi, na.rm = TRUE)
    if (is.finite(I_max) && I_max > 0) {
      threshold <- plateau_tol * I_max
      for (lvl in seq(2L, L - 1L)) {
        gain <- mi[lvl] - mi[lvl - 1L]
        if (is.finite(gain) && gain < threshold) {
          picked <- lvl
          break
        }
      }
    }
  }

  if (verbose) {
    message(sprintf(
      "[.criterion_mi_plateau] L=%d tol=%g picked=%d",
      L,
      plateau_tol,
      picked
    ))
  }

  n_blocks <- vapply(
    memb,
    function(m) length(unique(m)),
    integer(1L)
  )

  diag_df <- data.frame(
    level = seq_len(L),
    n_blocks = n_blocks,
    mdl = NA_real_,
    mi_to_next = mi,
    score = mi,
    stringsAsFactors = FALSE
  )

  list(
    level = as.integer(picked),
    method = "mi_plateau",
    diagnostics = diag_df
  )
}

# 9. Kneedle elbow detection -------------------------------------------------

#' @keywords internal
#' @noRd
.kneedle <- function(y) {
  y <- as.numeric(y)
  L <- length(y)
  if (L <= 2L) {
    return(L)
  }
  if (all(y == y[1])) {
    return(1L)
  }

  x_norm <- (seq_len(L) - 1L) / (L - 1L)
  y_min <- min(y)
  y_max <- max(y)
  y_norm <- (y - y_min) / (y_max - y_min)

  decreasing <- y_norm[L] < y_norm[1]
  d <- if (decreasing) {
    y_norm - (1 - x_norm)
  } else {
    y_norm - x_norm
  }

  which.max(abs(d))
}

# 10. Discrete second-difference elbow ---------------------------------------

#' @keywords internal
#' @noRd
.max_second_diff <- function(y) {
  y <- as.numeric(y)
  L <- length(y)
  if (L <= 2L) {
    return(L)
  }
  d2 <- diff(diff(y))
  which.max(abs(d2)) + 1L
}

# 11. Cliques -> membership vector ------------------------------------------

#' @keywords internal
#' @noRd
.cliques_to_membership <- function(cliques, n) {
  m <- seq_len(n)
  for (i in seq_along(cliques)) {
    idx <- cliques[[i]]
    if (length(idx) >= 2L) {
      m[idx] <- n + i
    }
  }
  as.integer(factor(m))
}

# 12. Sample mutual information (bits) --------------------------------------

#' @keywords internal
#' @noRd
.mutual_information <- function(a, b) {
  tab <- table(a, b)
  n_total <- sum(tab)
  if (n_total == 0L) {
    return(0)
  }
  p_ab <- tab / n_total
  p_a <- rowSums(p_ab)
  p_b <- colSums(p_ab)
  nz <- p_ab > 0
  outer_ab <- outer(p_a, p_b)
  sum(p_ab[nz] * log2(p_ab[nz] / outer_ab[nz]))
}

# 13. Modularity trace for moneca_fast() output -----------------------------

#' @keywords internal
#' @noRd
.fast_compute_modularity <- function(obj) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(
      "auto_segment_levels(method = 'mdl') on moneca_fast() output ",
      "requires the 'igraph' package.",
      call. = FALSE
    )
  }

  seg_list <- obj$segment.list
  L <- length(seg_list)
  n_core <- length(seg_list[[1]])

  mx_full <- obj$mat.list[[1]]
  # strip margins if they were appended (mat is (n_core + 1) square)
  core_mx <- if (isTRUE(obj$margins_added) && nrow(mx_full) == n_core + 1L) {
    mx_full[seq_len(n_core), seq_len(n_core), drop = FALSE]
  } else {
    mx_full[seq_len(n_core), seq_len(n_core), drop = FALSE]
  }

  # igraph likes a dense matrix here; at classification scale this is fine
  core_dense <- as.matrix(core_mx)

  g <- igraph::graph_from_adjacency_matrix(
    core_dense,
    mode = "directed",
    weighted = TRUE,
    diag = FALSE
  )
  w <- igraph::E(g)$weight

  # Level 1 is atomic (each node its own block) — modularity is 0 and
  # contributes no elbow information. Skip it; diagnostics align to
  # segment.list levels 2..L.
  if (L < 2L) {
    return(list(modularity = numeric(0), n_blocks = integer(0)))
  }

  mods <- numeric(L - 1L)
  nbs <- integer(L - 1L)
  for (l in seq(2L, L)) {
    memb <- .cliques_to_membership(seg_list[[l]], n_core)
    mods[l - 1L] <- tryCatch(
      igraph::modularity(g, membership = memb, weights = w),
      error = function(e) NA_real_
    )
    nbs[l - 1L] <- length(unique(memb))
  }

  list(modularity = mods, n_blocks = nbs)
}
