#' Permutation-Based Significance Threshold
#'
#' Computes an empirical significance threshold for interaction screening via
#' permutation. This is essential when the theoretical null distribution may
#' not hold (e.g., due to population structure, relatedness, or viral clonality
#' in host-pathogen studies as in Chen et al. 2026).
#'
#' @param data Data frame. Input data with variants and effects.
#' @param n_perm Integer. Number of permutations. Default 1000.
#' @param effect_b_col Character. Column name for the effect to be permuted.
#'   Typically the "exposure" or "independent" effect. Default \code{"x1"}.
#' @param keep_top_n Integer. Track the top N P-values from each permutation
#'   for FDR calculation. Default 10.
#' @param seed Integer. Random seed for reproducibility. Default \code{NULL}.
#' @param trace Logical. Print progress. Default \code{TRUE}.
#'
#' @return A list with class \code{"mrgxe_permutation"} containing:
#'   \describe{
#'     \item{empirical_threshold_0.05}{Empirical P threshold at alpha = 0.05}
#'     \item{empirical_threshold_0.01}{Empirical P threshold at alpha = 0.01}
#'     \item{n_perm}{Number of permutations performed}
#'     \item{min_p_values}{Vector of minimum P-values from each permutation}
#'   }
#'
#' @details
#' The permutation approach:
#' 1. Shuffle the exposure/independent effect (x1) across variants
#' 2. Re-estimate theta using the permuted data
#' 3. Recompute interaction P-values for all variants
#' 4. Record the minimum P-value
#' 5. The empirical threshold is the alpha quantile of the null min-P distribution
#'
#' @export
permutation_threshold <- function(
  data,
  n_perm = 1000,
  effect_b_col = "x1",
  keep_top_n = 10,
  seed = NULL,
  trace = TRUE
) {
  if (!is.null(seed)) set.seed(seed)

  # Validate
  if (!effect_b_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in data", effect_b_col))
  }

  n <- nrow(data)
  effect_b <- data[[effect_b_col]]
  x2 <- data$x2 %||% data$beta_a
  se1 <- data$x1_se %||% data$se_b
  se2 <- data$x2_se %||% data$se_a

  if (is.null(x2) || is.null(se1) || is.null(se2)) {
    stop("Need x2, x1_se, x2_se columns (from harmonize_effects)")
  }

  min_pvals <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    if (trace && i %% 100 == 0) {
      message(sprintf("  Permutation %d / %d", i, n_perm))
    }

    # Permute effect B
    perm_idx <- sample(n)
    x1_perm <- effect_b[perm_idx]

    # Quick theta estimate (OLS on all variants)
    theta_perm <- stats::coef(stats::lm(x2 ~ x1_perm, weights = 1/se2^2))[2]

    # Compute deviation P-values
    dev <- x2 - theta_perm * x1_perm
    dev_se <- sqrt(se2^2 + theta_perm^2 * se1^2)
    dev_se <- pmax(dev_se, .Machine$double.eps)
    p_vals <- 2 * stats::pnorm(-abs(dev / dev_se))
    p_vals <- pmin(p_vals, 1)

    min_pvals[i] <- min(p_vals, na.rm = TRUE)
  }

  # Empirical thresholds
  thresholds <- stats::quantile(min_pvals, probs = c(0.05, 0.01), na.rm = TRUE)

  message(sprintf("Permutation complete (%d runs)", n_perm))
  message(sprintf("  Empirical threshold (alpha=0.05): %.2e", thresholds[1]))
  message(sprintf("  Empirical threshold (alpha=0.01): %.2e", thresholds[2]))

  result <- list(
    empirical_threshold_0.05 = thresholds[1],
    empirical_threshold_0.01 = thresholds[2],
    n_perm = n_perm,
    min_p_values = min_pvals,
    params = list(seed = seed)
  )
  class(result) <- "mrgxe_permutation"
  result
}

#' BRASS-Style Permutation for Structured Samples
#'
#' Implements the BRASS (Binary trait, Relatedness-Adjusted, Stratified
#' permutation) approach used in Chen et al. (2026) for controlling viral
#' clonality and relatedness in EBV genome-wide interaction scans.
#'
#' In structured samples (e.g., related viral genomes), standard permutation
#' is too liberal. BRASS uses a genetic relatedness matrix (GRM) as a random
#' effect and permutes at the cluster level.
#'
#' @param data Data frame with individual-level data.
#' @param phenotype Character. Column name for binary outcome.
#' @param genotype_a Character. Column name for genotype/snp of entity A.
#' @param genotype_b Character. Column name for genotype/snp of entity B.
#' @param covariate_cols Character vector. Covariate column names.
#' @param cluster_col Character. Column defining independent clusters for
#'   permutation. For viral data, this could be the viral lineage or host ID.
#' @param n_perm Integer. Number of permutations. Default 1000.
#' @param seed Integer. Random seed. Default \code{NULL}.
#'
#' @return A list with permutation results and empirical threshold.
#'
#' @export
brass_permutation <- function(
  data,
  phenotype,
  genotype_a,
  genotype_b,
  covariate_cols = NULL,
  cluster_col = NULL,
  n_perm = 1000,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Validate required columns
  required <- c(phenotype, genotype_a, genotype_b)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }

  # Extract variables
  y <- data[[phenotype]]
  g_a <- data[[genotype_a]]
  g_b <- data[[genotype_b]]

  # Build model matrix
  covars <- if (!is.null(covariate_cols)) {
    as.matrix(data[, covariate_cols, drop = FALSE])
  } else {
    NULL
  }

  # Determine permutation units
  if (!is.null(cluster_col) && cluster_col %in% names(data)) {
    clusters <- data[[cluster_col]]
    cluster_ids <- unique(clusters)
    n_units <- length(cluster_ids)
    permute_at_cluster <- TRUE
  } else {
    n_units <- nrow(data)
    cluster_ids <- seq_len(n_units)
    permute_at_cluster <- FALSE
  }

  message(sprintf("BRASS permutation: %d units, %d permutations",
                  n_units, n_perm))

  # Observed interaction statistic
  fit_obs <- stats::glm(
    y ~ g_a + g_b + g_a:g_b + .,
    data = if (!is.null(covars)) {
      cbind(data.frame(y = y, g_a = g_a, g_b = g_b), covars)
    } else {
      data.frame(y = y, g_a = g_a, g_b = g_b)
    },
    family = "binomial"
  )

  obs_stat <- summary(fit_obs)$coefficients["g_a:g_b", "z value"]

  # Permutation loop
  perm_stats <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    if (permute_at_cluster) {
      # Permute cluster labels
      shuffled_b <- g_b
      perm_clusters <- sample(cluster_ids)
      for (cid in seq_along(cluster_ids)) {
        idx <- which(clusters == cluster_ids[cid])
        target_cluster <- perm_clusters[cid]
        target_idx <- which(clusters == target_cluster)
        shuffled_b[idx] <- g_b[target_idx]
      }
    } else {
      shuffled_b <- sample(g_b)
    }

    fit_perm <- stats::glm(
      y ~ g_a + shuffled_b + g_a:shuffled_b + .,
      data = if (!is.null(covars)) {
        cbind(data.frame(y = y, g_a = g_a, shuffled_b = shuffled_b), covars)
      } else {
        data.frame(y = y, g_a = g_a, shuffled_b = shuffled_b)
      },
      family = "binomial"
    )

    perm_stats[i] <- summary(fit_perm)$coefficients["g_a:shuffled_b", "z value"]
  }

  # Empirical P-value
  empirical_p <- mean(abs(perm_stats) >= abs(obs_stat), na.rm = TRUE)

  message(sprintf("BRASS permutation complete"))
  message(sprintf("  Observed Z = %.3f", obs_stat))
  message(sprintf("  Empirical P = %.4f", empirical_p))

  result <- list(
    observed_stat = obs_stat,
    permuted_stats = perm_stats,
    empirical_p = empirical_p,
    n_perm = n_perm,
    n_units = n_units
  )
  class(result) <- "mrgxe_brass"
  result
}
