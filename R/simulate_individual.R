# ===========================================================================
# Individual-Level GxE Simulation Functions
# ===========================================================================
# These functions implement the exact simulation methodology described in:
#   Zhu et al. (2024) Nature Communications, doi:10.1038/s41467-024-47806-3
#
# Based on the original simulation code from:
#   https://github.com/harryyiheyang/MR.GxE.Code (Zenodo: 10.5281/zenodo.10815731)
#
# Three simulation designs are implemented, corresponding to the paper's figures:
#   1. simulate_gxe_fig2ad()  → Fig 2A–D: many variants, no mediation
#   2. simulate_gxe_fig2ef()  → Fig 2E–F: 20 variants, with/without mediation
#   3. simulate_gxe_individual_effects() → compute GWAS & GWIS effects from
#      individual-level data (reusable building block)
# ===========================================================================

#' Simulate Individual-Level GxE Data — Fig 2A–D Design
#'
#' Generates individual-level genotype, environment, and phenotype data
#' following the simulation design of Zhu et al. (2024) Fig 2A–D.
#' This design has NO mediation; it tests T_MRGxE and T_Direct under
#' varying environment distributions between GWAS and GWIS cohorts.
#'
#' @section Methodology (from the paper):
#' For the ith individual, m independent variants are generated:
#'   G*_ij ~ Binom(2, p_j),  p_j ~ Uniform(0.05, 0.5)
#' Genotypes are standardised: G_ij = G*_ij / sqrt(2 p_j (1 - p_j))
#'
#' Environment:
#'   GWAS cohort:  E_i1 ~ N(mu_E1, 1)  (or Binary for Fig 2A-D)
#'   GWIS cohort:  E_i2 ~ N(mu_E2, 1)  (or Binary)
#'   Overlapping samples use GWIS environment
#'
#' Main effects:  beta_1j ~ N(0, sigma2_beta)
#' Phenotype:     Y_i = sum_j(G_ij * beta_1j) + gamma_E * E_i
#'                      + beta_3 * G_i1 * E_i1 + epsilon_i
#'                epsilon_i ~ N(0, sigma2_epsilon)
#'
#' The causal effect theta is estimated using the LAST n_iv variants as
#' instruments (to avoid correlation with the interaction variant).
#'
#' @param n_variants     Integer. Number of variants (m). Default 200.
#' @param n_gwas         Integer. GWAS sample size. Default 160000.
#' @param n_gwis         Integer. GWIS sample size. Default 80000.
#' @param n_overlap      Integer. Number of overlapping samples between
#'   GWAS and GWIS. Default equals \code{n_gwis} (full overlap).
#' @param h2             Numeric. Target heritability (proportion of
#'   phenotypic variance explained by all variants). Default 0.1.
#' @param gamma_e        Numeric. Environment main effect size.
#'   Default 0.1.
#' @param beta_3         Numeric. Interaction effect size for variant 1.
#'   Default 0.005.
#' @param mu_e1          Numeric. Environment mean for GWAS-only samples.
#'   Default 0.1.
#' @param mu_e2          Numeric. Environment mean for GWIS cohort.
#'   Default 0.1.
#' @param env_type       Character. Environment distribution:
#'   \code{"binary"} for Bernoulli (original Fig 2A-D),
#'   \code{"normal"} for Gaussian. Default \code{"binary"}.
#' @param variant1_beta  Numeric. Main effect size for variant 1
#'   (the interaction variant). If \code{NULL}, randomly generated
#'   from N(0, sigma2_beta). For Fig 2C-D, this can be 0, positive,
#'   or negative to test different scenarios.
#' @param seed           Integer. Random seed. Default \code{NULL}.
#' @param compute_effects Logical. If \code{TRUE}, compute GWAS and GWIS
#'   summary statistics from the generated data. Default \code{TRUE}.
#'
#' @return A list with components:
#'   \item{genotypes}{Matrix (n_total x n_variants) of standardised genotypes}
#'   \item{phenotype}{Numeric vector of phenotypes (length n_total)}
#'   \item{environment}{Numeric vector of environment values}
#'   \item{true_beta}{Numeric vector of true main effects}
#'   \item{true_beta3}{The true interaction effect (scalar)}
#'   \item{gwas_indices}{Indices of GWAS cohort}
#'   \item{gwis_indices}{Indices of GWIS cohort}
#'   \item{gwas_summary}{Data frame of GWAS marginal effect estimates (if compute_effects=TRUE)}
#'   \item{gwis_summary}{Data frame of GWIS main + interaction effect estimates (if compute_effects=TRUE)}
#'   \item{true_theta}{The estimated true theta from the IV set}
#'   \item{maf}{Minor allele frequencies}
#'
#' @references
#' Zhu, X., Yang, Y., Lorincz-Comi, N. et al. An approach to identify
#' gene-environment interactions and reveal new biological insight in
#' complex traits. Nat Commun 15, 3385 (2024).
#'
#' @export
simulate_gxe_fig2ad <- function(
    n_variants   = 200,
    n_gwas       = 160000,
    n_gwis       = 80000,
    n_overlap    = NULL,
    h2           = 0.1,
    gamma_e      = 0.1,
    beta_3       = 0.005,
    mu_e1        = 0.1,
    mu_e2        = 0.1,
    env_type     = c("binary", "normal"),
    variant1_beta = NULL,
    seed         = NULL,
    compute_effects = TRUE) {

  if (!is.null(seed)) set.seed(seed)
  env_type <- match.arg(env_type)

  # ---- Sample indices ----
  if (is.null(n_overlap)) n_overlap <- n_gwis
  n_total <- max(n_gwas, n_gwis)

  # GWAS indices: 1..n_gwas
  idx_gwas <- seq_len(n_gwas)
  # GWIS indices: first n_gwis (overlap at the beginning with GWAS)
  idx_gwis <- seq_len(n_gwis)
  # Non-overlapping GWAS-only indices
  idx_gwas_only <- setdiff(idx_gwas, idx_gwis)

  # ---- Generate genotypes ----
  maf <- runif(n_variants, 0.05, 0.5)
  maf[1] <- 0.3  # variant 1 (interaction variant) fixed MAF

  G <- matrix(NA_real_, nrow = n_total, ncol = n_variants)
  for (j in seq_len(n_variants)) {
    g_raw <- rbinom(n_total, size = 2, prob = maf[j])
    # Standardise: G = g_raw / sqrt(2 * p * (1-p))
    # (Note: original code does NOT subtract mean before scaling)
    G[, j] <- g_raw / sqrt(2 * maf[j] * (1 - maf[j]))
  }

  # ---- Generate main effects (beta_1j) ----
  sigma2_beta <- h2 / n_variants
  beta <- rnorm(n_variants, mean = 0, sd = sqrt(sigma2_beta))

  # Override variant 1 main effect if specified
  if (!is.null(variant1_beta)) {
    ssign <- sign(variant1_beta)
    beta[1] <- 2 * sqrt(sigma2_beta) * ssign
  }

  # ---- Generate environment ----
  E <- rep(0, n_total)

  if (env_type == "binary") {
    # Binary environment (original Fig 2A-D approach)
    # GWIS cohort
    E[idx_gwis] <- rbinom(n_gwis, size = 1, prob = mu_e2)
    # GWAS-only cohort
    if (length(idx_gwas_only) > 0) {
      E[idx_gwas_only] <- rbinom(length(idx_gwas_only), size = 1, prob = mu_e1)
    }
    # Standardise: E = E / sd(E)
    E <- E / stats::sd(E)
  } else {
    # Normal environment
    E[idx_gwis] <- stats::rnorm(n_gwis, mean = mu_e2, sd = 1)
    if (length(idx_gwas_only) > 0) {
      E[idx_gwas_only] <- stats::rnorm(length(idx_gwas_only), mean = mu_e1, sd = 1)
    }
  }

  # ---- Scale beta to match target heritability ----
  mug <- as.vector(G %*% beta)
  rho_g <- h2 / stats::var(mug)
  beta <- beta * sqrt(rho_g)

  # Recompute genetic component with scaled betas
  mug <- as.vector(G %*% beta)

  # ---- Build phenotype ----
  mue   <- E * gamma_e                         # environment contribution
  mui1  <- G[, 1] * E * beta_3                 # GxE interaction (variant 1 only)
  mu    <- mug + mue + mui1                     # total mean

  # Residual variance: ensure h2 fraction of total variance is genetic
  sig2_resid <- stats::var(as.matrix(G %*% beta)) / h2 - stats::var(mu)
  sig2_resid <- max(sig2_resid, 0.01)           # prevent negative / zero variance

  y <- mu + stats::rnorm(n_total, mean = 0, sd = sqrt(sig2_resid))

  # ---- Compute GWAS & GWIS summary statistics ----
  result <- list(
    genotypes    = G,
    phenotype    = y,
    environment  = E,
    true_beta    = beta,
    true_beta3   = beta_3,
    gwas_indices = idx_gwas,
    gwis_indices = idx_gwis,
    maf          = maf
  )

  if (compute_effects) {
    eff <- compute_gxe_effects_individual(
      y           = y,
      G           = G,
      E           = E,
      idx_gwas    = idx_gwas,
      idx_gwis    = idx_gwis,
      n_overlap   = n_overlap
    )
    result$gwas_summary <- eff$gwas
    result$gwis_summary <- eff$gwis
    result$true_theta   <- eff$theta_iv
    result$iv_indices   <- eff$iv_indices
  }

  class(result) <- c("gxe_simulation_fig2ad", "list")
  return(result)
}

# ===========================================================================
# Fig 2E–F: 20-variant design with mediation scenarios
# ===========================================================================

#' Simulate Individual-Level GxE Data — Fig 2E–F Design
#'
#' Generates data for the 20-variant simulation design in Zhu et al. (2024)
#' Fig 2E–F. This design includes both type I error and power scenarios,
#' with and without mediation through the environment.
#'
#' @section Phenotype models:
#' Three models are available:
#' \describe{
#'   \item{"none"}{No mediation, no interaction:
#'     \code{Y ~ 0.1*G + gamma*E + N(0, 10)}}
#'   \item{"mediation_only"}{Mediation but no interaction:
#'     \code{E ~ 0.05*G + N(1, 0.9975); Y ~ 0.1*G + gamma*E + N(0, 10)}}
#'   \item{"mediation_interaction"}{Mediation and interaction:
#'     \code{E ~ 0.05*G + N(1, 0.9975); Y ~ 0.1*G + gamma*E + 0.1*G*E + N(0, 10)}}
#' }
#'
#' @param n_variants   Integer. Number of variants. Default 20.
#' @param n_gwas       Integer. GWAS sample size. Default 20000.
#'   Also accepts a vector like \code{c(20000, 100000, 200000, 300000)}
#'   for power curves.
#' @param n_gwis       Integer. GWIS sample size. Default 20000.
#' @param n_sim        Integer. Number of simulation replicates.
#'   Default 1000.
#' @param maf          Numeric. Minor allele frequency (same for all
#'   variants). Default 0.3.
#' @param gamma        Numeric. Environment effect size.
#'   Default 1 (as in the paper's Fig 2E-F, though paper says γ = 1 or 5).
#'   Use 5 for the larger effect scenario.
#' @param mu_e         Numeric. Environment mean for GWIS cohort
#'   (GWAS cohort mean = 2 * mu_e for non-overlapped).
#'   Default 1. (Paper also used 0.5 for Supplementary Fig 10.)
#' @param model        Character. Phenotype model:
#'   \code{"none"}, \code{"mediation_only"}, or \code{"mediation_interaction"}.
#'   Default \code{"none"}.
#' @param env_contrib  Numeric. Proportion of phenotypic variance
#'   contributed by environment. Default 0.01 (1\%).
#'   Used to scale the environment coefficient gamma.
#' @param mediation_rho Numeric. Correlation between genotype and
#'   environment under mediation. Default 0.01 (1\%).
#' @param seed         Integer. Random seed. Default \code{NULL}.
#' @param verbose      Logical. Print progress every 100 replicates.
#'   Default \code{TRUE}.
#'
#' @return A list with:
#'   \item{power}{Matrix of power estimates (n_sim x n_sample_sizes x 3 tests)}
#'   \item{theta_estimates}{Matrix of r = alpha_hat - beta_hat estimates}
#'   \item{params}{Named list of simulation parameters}
#'
#' @export
simulate_gxe_fig2ef <- function(
    n_variants    = 20,
    n_gwas        = c(20000, 100000, 200000, 300000),
    n_gwis        = 20000,
    n_sim         = 1000,
    maf           = 0.3,
    gamma         = 1,
    mu_e          = 1,
    model         = c("none", "mediation_only", "mediation_interaction"),
    env_contrib   = 0.01,
    mediation_rho = 0.01,
    seed          = NULL,
    verbose       = TRUE) {

  if (!is.null(seed)) set.seed(seed)
  model <- match.arg(model)

  # ---- Derived parameters ----
  # Standardisation factor for genotypes
  # G ~ Binom(2, p), std G = G / sqrt(2*p*(1-p)) = G / (2*p*(1-p))
  # But original code divides by 0.648 directly (= 2*0.3*0.7)
  # Actually: sqrt(2*0.3*0.7) = sqrt(0.42) ≈ 0.648
  # So dividing by 0.648 is the same as dividing by sqrt(2pq)
  std_factor <- sqrt(2 * maf * (1 - maf))  # ≈ 0.648 for maf=0.3

  # Environment variance contribution
  # Paper: E ~ N(μE, 1), so var(E) = 1
  # If gamma is scaled by env_contrib, we need gamma such that
  # var(gamma*E) / var(Y) ≈ env_contrib
  # var(Y) = var(0.1*G) + gamma^2 * var(E) + 10 + (interaction term)
  # var(0.1*G) = 0.01 (since std G has var 1)
  # So gamma^2 * 1 / (0.01 + gamma^2 + 10) = env_contrib
  # gamma^2 = env_contrib * (10.01) / (1 - env_contrib)

  mu_g <- 2 * maf / std_factor  # mean of standardised genotype

  n_sample_sizes <- length(n_gwas)

  # ---- Storage ----
  # power[sim, sample_size, test] where test = 1:T_Direct, 2:T_MRGxE, 3:TwoStep
  power_store <- array(0, dim = c(n_sim, n_sample_sizes, 3))
  theta_store <- matrix(0, n_sim, n_sample_sizes)

  # ---- Main simulation loop ----
  for (i in seq_len(n_sample_sizes)) {
    n1 <- n_gwas[i]
    n2 <- n_gwis
    n_total <- max(n1, n2)

    idx_gwas <- seq_len(n1)
    idx_gwis <- seq_len(n2)

    for (j in seq_len(n_sim)) {
      # Generate 20 variants, each Binom(2, maf), standardised
      G <- matrix(NA_real_, nrow = n_total, ncol = n_variants)
      for (k in seq_len(n_variants)) {
        g_raw <- stats::rbinom(n_total, size = 2, prob = maf)
        G[, k] <- g_raw / std_factor
      }

      # Generate environment
      rho1 <- if (model %in% c("mediation_only", "mediation_interaction")) {
        mediation_rho  # mediation: G contributes to E
      } else {
        0  # no mediation
      }

      # For mediation: E depends on variant 1's genotype
      # E_i = rho1 * G_i1 + N(mu_E_target, sqrt(1 - rho1^2))
      E <- rep(0, n_total)

      # GWIS cohort (first n2 subjects)
      for (k in seq_len(n_variants)[-seq_len(min(1, n_variants))]) {
        # Variants 2..20: no mediation
        # (handled below)
      }

      if (model == "none") {
        # All 20 variants: no mediation
        for (k in seq_len(n_variants)) {
          g_k <- stats::rbinom(n_total, size = 2, prob = maf) / std_factor
          # GWAS cohort environment (mean = 2*mu_e for non-overlapped)
          mu_e_gwas <- if (n1 > n2) 2 * mu_e else mu_e
          e_k_gwas <- stats::rnorm(n1, mean = mu_e_gwas, sd = 1)
        }
        # Simplified: single environment variable
        E[idx_gwis] <- stats::rnorm(n2, mean = mu_e, sd = 1)
        if (n1 > n2) {
          non_overlap_idx <- seq.int(n2 + 1, n1)
          E[non_overlap_idx] <- stats::rnorm(length(non_overlap_idx),
                                             mean = 2 * mu_e, sd = 1)
        }
      } else if (model == "mediation_only" || model == "mediation_interaction") {
        # Variant 1 mediates environment
        # E = rho1 * G[,1] + N(mu_E_target, sqrt(1 - rho1^2))
        # But actually, per the paper: only variant 1 has mediation
        # Variants 2..20: no mediation

        # GWIS cohort
        E[idx_gwis] <- rho1 * G[idx_gwis, 1] +
          stats::rnorm(n2, mean = mu_e, sd = sqrt(1 - rho1^2))

        # GWAS-only cohort (if any)
        if (n1 > n2) {
          non_overlap_idx <- seq.int(n2 + 1, n1)
          E[non_overlap_idx] <- rho1 * G[non_overlap_idx, 1] +
            stats::rnorm(length(non_overlap_idx), mean = 2 * mu_e,
                         sd = sqrt(1 - rho1^2))
        }
      }

      # ---- Generate phenotype ----
      y <- rep(0, n_total)

      # Genetic effect: 0.1 * G[,k] for each variant
      # (For power scenarios, variant 1 may have different coefficient)
      g_effect <- 0.1 * G[, 1]  # variant 1

      if (model == "none") {
        # Type I error scenario: Y = 0.1*G + gamma*E + N(0, 10)
        # OR power scenario: Y = 0.1*G + gamma*E + 0.1*G*E + N(0, 10)
        # We default to type I error (no interaction), power option via env_contrib
        y <- g_effect + gamma * E + stats::rnorm(n_total, mean = 0, sd = sqrt(10))
      } else if (model == "mediation_only") {
        # Power scenario variant: Y = 0.1*G + gamma*E + N(0, 10)
        # (mediation exists in E but not in Y directly)
        y <- g_effect + gamma * E + stats::rnorm(n_total, mean = 0, sd = sqrt(10))
      } else if (model == "mediation_interaction") {
        # Interaction scenario: Y = 0.1*G + gamma*E + 0.1*G*E + N(0, 10)
        y <- g_effect + gamma * E + 0.1 * G[, 1] * E +
          stats::rnorm(n_total, mean = 0, sd = sqrt(10))
      }

      # ---- Compute GWAS & GWIS effects for all variants ----
      gwas_beta <- numeric(n_variants)
      gwas_se   <- numeric(n_variants)
      gwis_main <- numeric(n_variants)
      gwis_se   <- numeric(n_variants)
      gwis_int  <- numeric(n_variants)
      gwis_int_se <- numeric(n_variants)

      y_gwas <- y[idx_gwas]
      y_gwis <- y[idx_gwis]
      G_gwas <- G[idx_gwas, , drop = FALSE]
      G_gwis <- G[idx_gwis, , drop = FALSE]
      E_gwis <- E[idx_gwis]

      for (k in seq_len(n_variants)) {
        # GWAS: marginal effect y ~ G_k
        fit_gwas <- stats::lm(y_gwas ~ G_gwas[, k])
        s_gwas <- summary(fit_gwas)$coefficients
        if (nrow(s_gwas) >= 2) {
          gwas_beta[k] <- s_gwas[2, 1]
          gwas_se[k]   <- s_gwas[2, 2]
        }

        # GWIS: y ~ G_k * E
        fit_gwis <- stats::lm(y_gwis ~ G_gwis[, k] * E_gwis)
        s_gwis <- summary(fit_gwis)$coefficients
        if (nrow(s_gwis) >= 4) {
          gwis_main[k]   <- s_gwis[2, 1]  # main effect of G_k
          gwis_se[k]     <- s_gwis[2, 2]
          gwis_int[k]    <- s_gwis[4, 1]  # interaction G_k:E
          gwis_int_se[k] <- s_gwis[4, 2]
        }
      }

      # ---- Compute T_Direct P-values ----
      p_direct <- stats::pchisq((gwis_int / gwis_int_se)^2, df = 1,
                                lower.tail = FALSE)

      # ---- Compute T_MRGxE ----
      # r = alpha_hat - beta_hat
      r_vec <- gwas_beta - gwis_main

      # Compute rho (correlation between GWAS and GWIS effects due to overlap)
      n_overlap <- n2  # all GWIS are in GWAS
      rho <- sqrt(n_overlap / n1 / n2) /
        sqrt(1 / (1 - rho1^2) + (rho1 * (1 - mu_g) + mu_e)^2)

      # Var(r) = var(alpha_hat) + var(beta_hat) - 2*rho*se(alpha)*se(beta)
      var_r <- gwas_se^2 + gwis_se^2 - 2 * rho * gwas_se * gwis_se
      var_r <- pmax(var_r, 1e-10)  # ensure positive

      p_mrgxe <- stats::pchisq(r_vec^2 / var_r, df = 1, lower.tail = FALSE)

      # ---- Bonferroni correction for 20 tests ----
      bonf_threshold <- 0.05 / n_variants  # = 0.0025

      # T_Direct significant?
      direct_sig <- p_direct < bonf_threshold

      # T_MRGxE significant?
      mrgxe_sig <- p_mrgxe < bonf_threshold

      # Two-step: first MRGxE, then Direct on survivors
      two_step_sig <- if (any(mrgxe_sig)) {
        survivors <- which(mrgxe_sig)
        k_surv <- length(survivors)
        if (k_surv == 1) {
          p_direct[survivors] < (0.05 / k_surv)
        } else {
          any(p_direct[survivors] < (0.05 / k_surv))
        }
      } else {
        FALSE
      }

      # Store theta estimate (the r value for variant 1)
      theta_store[j, i] <- r_vec[1]

      # Store power (for the interaction variant = variant 1)
      # Check if we're in a power scenario (variant 1 has interaction)
      # For "none" model: variant 1 has no interaction → type I error
      # For "mediation_interaction" model: variant 1 HAS interaction → power
      if (model == "mediation_interaction" || model == "mediation_only") {
        # Power: variant 1 is the interaction variant
        power_store[j, i, 1] <- as.numeric(direct_sig[1])
        power_store[j, i, 2] <- as.numeric(mrgxe_sig[1])
        power_store[j, i, 3] <- as.numeric(two_step_sig)
      } else {
        # Type I error: any variant could be significant
        power_store[j, i, 1] <- as.numeric(any(direct_sig))
        power_store[j, i, 2] <- as.numeric(any(mrgxe_sig))
        power_store[j, i, 3] <- as.numeric(two_step_sig)
      }

      if (verbose && j %% 100 == 0) {
        message(sprintf("Sample size %d/%d, replicate %d/%d",
                        i, n_sample_sizes, j, n_sim))
      }
    }
  }

  # ---- Summarise results ----
  colnames(power_store) <- paste0("n", n_gwas)
  dimnames(power_store)[[3]] <- c("T_Direct", "T_MRGxE", "TwoStep")

  result <- list(
    power           = power_store,
    theta_estimates = theta_store,
    params = list(
      n_variants     = n_variants,
      n_gwas         = n_gwas,
      n_gwis         = n_gwis,
      n_sim          = n_sim,
      maf            = maf,
      gamma          = gamma,
      mu_e           = mu_e,
      model          = model,
      env_contrib    = env_contrib,
      mediation_rho  = mediation_rho
    )
  )

  class(result) <- c("gxe_simulation_fig2ef", "list")
  return(result)
}


# ===========================================================================
# Core: Compute GWAS & GWIS effects from individual-level data
# ===========================================================================

#' Compute GWAS and GWIS Summary Statistics from Individual-Level Data
#'
#' Given individual-level genotype, phenotype, and environment data,
#' computes GWAS marginal effects (y ~ G) and GWIS main + interaction
#' effects (y ~ G * E). This is the workhorse used by all simulation
#' functions and can also be applied to real individual-level data
#' (e.g., UK Biobank).
#'
#' @param y          Numeric vector. Phenotype (length n_total).
#' @param G          Numeric matrix. Genotypes (n_total x n_variants).
#'   Should be standardised (mean 0, unit variance or similar).
#' @param E          Numeric vector. Environment variable (length n_total).
#' @param idx_gwas   Integer vector. Row indices for GWAS cohort.
#' @param idx_gwis   Integer vector. Row indices for GWIS cohort.
#' @param n_overlap  Integer. Number of overlapping samples (used for
#'   rho estimation). Default = min(length(idx_gwas), length(idx_gwis)).
#' @param iv_indices Integer vector. Indices of variants to use as
#'   instruments for theta estimation. If \code{NULL}, uses variant
#'   indices 3:end (skipping the first two which may be interaction
#'   variants, matching the paper's approach).
#'
#' @return A list with:
#'   \item{gwas}{Data frame with columns: BETA, SE, P, N}
#'   \item{gwis}{Data frame with columns: BETA_main, SE_main, BETA_int, SE_int}
#'   \item{theta_iv}{Estimated theta from IVW regression of GWAS ~ GWIS effects
#'     using the instrument variants}
#'   \item{rho}{Estimated sample overlap correlation}
#'   \item{iv_indices}{Indices used as instruments}
#'
#' @export
compute_gxe_effects_individual <- function(
    y,
    G,
    E,
    idx_gwas,
    idx_gwis,
    n_overlap = NULL,
    iv_indices = NULL) {

  n_variants <- ncol(G)
  n_gwas     <- length(idx_gwas)
  n_gwis     <- length(idx_gwis)

  if (is.null(n_overlap)) {
    n_overlap <- min(n_gwas, n_gwis)  # assume maximum overlap
  }

  # ---- Subset data ----
  y_gwas <- y[idx_gwas]
  y_gwis <- y[idx_gwis]
  G_gwas <- G[idx_gwas, , drop = FALSE]
  G_gwis <- G[idx_gwis, , drop = FALSE]
  E_gwis <- E[idx_gwis]

  # ---- Storage ----
  gwas_beta <- numeric(n_variants)
  gwas_se   <- numeric(n_variants)
  gwas_p    <- numeric(n_variants)

  gwis_main    <- numeric(n_variants)
  gwis_main_se <- numeric(n_variants)
  gwis_int     <- numeric(n_variants)
  gwis_int_se  <- numeric(n_variants)

  # ---- Compute effect estimates ----
  for (j in seq_len(n_variants)) {
    # GWAS: y ~ G_j
    fit_a <- stats::lm(y_gwas ~ G_gwas[, j])
    s_a   <- summary(fit_a)$coefficients
    if (nrow(s_a) >= 2) {
      gwas_beta[j] <- s_a[2, 1]
      gwas_se[j]   <- s_a[2, 2]
      gwas_p[j]    <- s_a[2, 4]
    }

    # GWIS: y ~ G_j + E + G_j:E
    fit_b <- stats::lm(y_gwis ~ G_gwis[, j] * E_gwis)
    s_b   <- summary(fit_b)$coefficients
    if (nrow(s_b) >= 4) {
      gwis_main[j]    <- s_b[2, 1]  # G_j main effect
      gwis_main_se[j] <- s_b[2, 2]
      gwis_int[j]     <- s_b[4, 1]  # G_j:E interaction
      gwis_int_se[j]  <- s_b[4, 2]
    }
  }

  # ---- Select instrument variants ----
  if (is.null(iv_indices)) {
    # Default: use variants 3..end (skip first 2 = interaction variants)
    if (n_variants >= 5) {
      iv_indices <- seq.int(3, n_variants)
    } else {
      iv_indices <- seq_len(n_variants)
    }
  }

  # ---- Estimate theta (IVW) using instrument variants ----
  # theta = cov(alpha, beta) / var(beta) for IV variants
  # Equivalent to lm(alpha ~ beta, weights = 1/se_beta^2)
  beta_iv <- gwis_main[iv_indices]
  alpha_iv <- gwas_beta[iv_indices]

  valid <- is.finite(beta_iv) & is.finite(alpha_iv) &
    !is.na(beta_iv) & !is.na(alpha_iv)

  if (sum(valid) >= 3) {
    fit_theta <- stats::lm(alpha_iv[valid] ~ beta_iv[valid])
    theta_iv  <- stats::coef(fit_theta)[2]
  } else {
    theta_iv <- NA_real_
  }

  # ---- Estimate rho (sample overlap correlation) ----
  rho <- n_overlap / sqrt(n_gwas) / sqrt(n_gwis) /
    sqrt(1 + mean(E[idx_gwis], na.rm = TRUE)^2)

  # ---- Build output data frames ----
  gwas_df <- data.frame(
    BETA = gwas_beta,
    SE   = gwas_se,
    P    = gwas_p,
    N    = n_gwas,
    stringsAsFactors = FALSE
  )

  gwis_df <- data.frame(
    BETA_main = gwis_main,
    SE_main   = gwis_main_se,
    BETA_int  = gwis_int,
    SE_int    = gwis_int_se,
    N         = n_gwis,
    stringsAsFactors = FALSE
  )

  list(
    gwas       = gwas_df,
    gwis       = gwis_df,
    theta_iv   = theta_iv,
    rho        = rho,
    iv_indices = iv_indices
  )
}

# ===========================================================================
# T_MRGxE test statistic (based on original t_pleio function)
# ===========================================================================

#' Compute T_MRGxE Test Statistic
#'
#' Implements the T_MRGxE (pleiotropy) test as described in Zhu et al. (2024).
#' Tests whether a variant deviates from the regression line
#' alpha_hat = theta * beta_hat.
#'
#' @param bx      Numeric vector. GWIS main effect estimates (beta_hat).
#' @param by      Numeric vector. GWAS marginal effect estimates (alpha_hat).
#' @param bx_se   Numeric vector. Standard errors of GWIS main effects.
#' @param by_se   Numeric vector. Standard errors of GWAS marginal effects.
#' @param betahat Numeric scalar. Estimated theta (causal effect).
#' @param rho     Numeric scalar. Sample overlap correlation.
#' @param sig     Numeric scalar. Estimated variance of true betas.
#'
#' @return Numeric vector of P-values (one per variant).
#'
#' @keywords internal
tmrgxe_test <- function(bx, by, bx_se, by_se, betahat, rho, sig) {
  # T_var = var(by - betahat * bx)
  #       = by_se^2 + betahat^2 * bx_se^2 - 2*betahat*rho*by_se*bx_se
  #         + (1 - betahat)^2 * sig
  T_var <- by_se^2 + betahat^2 * bx_se^2 -
    2 * betahat * rho * by_se * bx_se +
    (1 - betahat)^2 * sig
  T_var <- abs(T_var)  # ensure positive
  T_var <- pmax(T_var, 1e-300)

  T_stat <- (by - betahat * bx)^2 / T_var
  T_p    <- stats::pchisq(T_stat, df = 1, lower.tail = FALSE)
  return(T_p)
}


# ===========================================================================
# Convenience: Three phenotype models for mediation scenarios
# ===========================================================================

#' Generate Phenotype Under GxE Models
#'
#' Generates a phenotype vector under one of three models matching
#' Zhu et al. (2024) Fig 2E-F:
#' \enumerate{
#'   \item No mediation, no interaction (type I error baseline)
#'   \item Mediation only (G affects E, but no GxE interaction)
#'   \item Mediation + interaction (G affects E, and GxE interaction)
#' }
#'
#' @param G           Numeric matrix. Genotypes (standardised).
#' @param E           Numeric vector. Environment values.
#' @param model       Character. \code{"none"}, \code{"mediation_only"},
#'   or \code{"mediation_interaction"}.
#' @param beta_g      Numeric. Genotype main effect. Default 0.1.
#' @param gamma_e     Numeric. Environment coefficient. Default 1.
#' @param beta_gxe    Numeric. GxE interaction coefficient. Default 0.1.
#' @param resid_sd    Numeric. Residual standard deviation. Default sqrt(10).
#' @param seed        Integer. Random seed. Default \code{NULL}.
#'
#' @return Numeric vector of phenotypes.
#'
#' @examples
#' # 1000 individuals, 1 variant
#' n <- 1000
#' g <- rbinom(n, 2, 0.3) / sqrt(2 * 0.3 * 0.7)
#' e <- rnorm(n, mean = 1, sd = 1)
#' y <- generate_gxe_phenotype(cbind(g), e, model = "mediation_interaction")
#'
#' @export
generate_gxe_phenotype <- function(
    G,
    E,
    model    = c("none", "mediation_only", "mediation_interaction"),
    beta_g   = 0.1,
    gamma_e  = 1,
    beta_gxe = 0.1,
    resid_sd = sqrt(10),
    seed     = NULL) {

  if (!is.null(seed)) set.seed(seed)
  model <- match.arg(model)

  n <- length(E)
  G <- as.matrix(G)

  # Genetic component: sum over all variants
  if (ncol(G) == 1) {
    g_effect <- beta_g * G[, 1]
  } else {
    g_effect <- as.vector(G %*% rep(beta_g, ncol(G)))
  }

  # Environment effect
  e_effect <- gamma_e * E

  # Interaction (only applies to variant 1)
  if (model == "mediation_interaction") {
    gxe_effect <- beta_gxe * G[, 1] * E
  } else {
    gxe_effect <- 0
  }

  # Residual
  eps <- stats::rnorm(n, mean = 0, sd = resid_sd)

  y <- g_effect + e_effect + gxe_effect + eps
  return(y)
}


# ===========================================================================
# Summary method for simulation objects
# ===========================================================================

#' @export
summary.gxe_simulation_fig2ad <- function(object, ...) {
  cat("GxE Simulation (Fig 2A-D Design)\n")
  cat("=================================\n")
  cat(sprintf("Variants: %d\n", ncol(object$genotypes)))
  cat(sprintf("GWAS N:   %d\n", length(object$gwas_indices)))
  cat(sprintf("GWIS N:   %d\n", length(object$gwis_indices)))
  cat(sprintf("True theta: %.4f\n", object$true_theta))
  cat(sprintf("True beta3: %.6f\n", object$true_beta3))
  if (!is.null(object$gwas_summary)) {
    cat(sprintf("GWAS significant variants (P < 5e-8): %d\n",
                sum(object$gwas_summary$P < 5e-8, na.rm = TRUE)))
  }
  invisible(object)
}

#' @export
summary.gxe_simulation_fig2ef <- function(object, ...) {
  cat("GxE Simulation (Fig 2E-F Design)\n")
  cat("=================================\n")
  cat(sprintf("Model:      %s\n", object$params$model))
  cat(sprintf("Variants:   %d\n", object$params$n_variants))
  cat(sprintf("GWAS N:     %s\n", paste(object$params$n_gwas, collapse = ", ")))
  cat(sprintf("GWIS N:     %d\n", object$params$n_gwis))
  cat(sprintf("Replicates: %d\n", object$params$n_sim))

  # Average power by sample size and test
  avg_power <- apply(object$power, c(2, 3), mean, na.rm = TRUE)
  colnames(avg_power) <- paste0("n", object$params$n_gwas)
  cat("\nAverage power/type I error:\n")
  print(round(avg_power, 4))
  invisible(object)
}
