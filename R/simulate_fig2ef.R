# ===========================================================================
# Individual-Level GxE Simulation Functions (Fig 2E-F Design)
# ===========================================================================
# Exact reproduction of Zhu et al. (2024) Nature Communications,
# Fig 2E-F and Supplementary Fig 10 simulation methodology.
#
# Based on original code: Simulation_Fig_2E_and_2F.R
# from https://github.com/harryyiheyang/MR.GxE.Code
# ===========================================================================

#' Simulate GxE with Mediation — Fig 2E-F Design (Exact Reproduction)
#'
#' Reproduces the exact simulation design from Zhu et al. (2024) Fig 2E-F
#' and Supplementary Fig 10. Tests 20 independent variants under six
#' scenarios combining type I error and power, with and without mediation.
#'
#' @section Simulation scenarios:
#' \describe{
#'   \item{Scenario 1 (type I)}{No mediation, E contributes 1\%.
#'     E independent of G. Y ~ 0.1*G + γ*E + N(0,10)}
#'   \item{Scenario 2 (type I)}{Mediation, E contributes 1\%.
#'     Variant 1: E ~ 0.01*G + N(2*μE, ~1); Y ~ 0.2*G + E + N(0,10).
#'     Variants 2-20: E independent of G; Y ~ 0.2*G + E + N(0,10)}
#'   \item{Scenario 3 (type I)}{Mediation, E contributes 5\%.
#'     Variant 1: E ~ 0.05*G + N(2*μE, ~1); Y ~ 0.2*G + √5*E + N(0,10).
#'     Variants 2-20: E independent of G; Y ~ 0.2*G + √5*E + N(0,10)}
#'   \item{Scenario 4 (power)}{No mediation, E contributes 1\%.
#'     Variant 1: Y ~ 0.1*G + γ*E + 0.1*G*E + N(0,10).
#'     Variants 2-20: Y ~ 0.1*G + γ*E + N(0,10)}
#'   \item{Scenario 5 (power)}{Mediation, E contributes 1\%.
#'     Variant 1: E ~ 0.01*G + N(2*μE, ~1); Y ~ 0.2*G + E + 0.2*G*E + N(0,10).
#'     Variants 2-20: E independent of G; Y ~ 0.2*G + E + N(0,10)}
#'   \item{Scenario 6 (power)}{Mediation, E contributes 5\%.
#'     Variant 1: E ~ 0.05*G + N(2*μE, ~1); Y ~ 0.1*G + √5*E + 0.2*G*E + N(0,10).
#'     Variants 2-20: E independent of G; Y ~ 0.2*G + √5*E + N(0,10)}
#' }
#'
#' @param n_variants  Integer. Number of variants tested. Default 20.
#' @param n_gwas      Integer vector. GWAS sample sizes.
#'   Default \code{c(20000, 100000, 200000, 300000)}.
#' @param n_gwis      Integer. GWIS sample size. Default 20000.
#' @param n_sim       Integer. Number of simulation replicates per
#'   sample size. Default 1000.
#' @param maf         Numeric. Minor allele frequency. Default 0.3.
#' @param gamma       Numeric. Environment effect (γ). Default 1.
#'   The paper tested γ = 1 and γ = 5.
#' @param mu_e        Numeric. Environment mean for GWIS cohort.
#'   Default 1. (Use 0.5 for Supplementary Fig 10.)
#' @param scenarios   Integer vector. Which scenarios to run (1-6).
#'   Default \code{1:6}. 1-3 are type I error, 4-6 are power.
#' @param seed        Integer. Random seed. Default \code{NULL}.
#' @param n_cores     Integer. Number of parallel cores. Default 1
#'   (sequential). Set to \code{detectCores()/2} for parallel.
#' @param verbose     Logical. Print progress. Default \code{TRUE}.
#'
#' @return A list with components:
#'   \item{results}{Array (n_scenarios x n_sample_sizes x 3), average
#'     type I error / power for T_Direct, T_MRGxE, and TwoStep tests}
#'   \item{scenario_labels}{Character vector of scenario descriptions}
#'   \item{params}{Named list of simulation parameters}
#'
#' @references
#' Zhu, X. et al. (2024) Nature Communications, 15:3385.
#'
#' @export
simulate_gxe_fig2ef_exact <- function(
    n_variants = 20,
    n_gwas     = c(20000, 100000, 200000, 300000),
    n_gwis     = 20000,
    n_sim      = 1000,
    maf        = 0.3,
    gamma      = 1,
    mu_e       = 1,
    scenarios  = 1:6,
    seed       = NULL,
    n_cores    = 1,
    verbose    = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  # ---- Constants from the paper ----
  # Genotype standardisation factor: sqrt(2*p*(1-p))
  std_factor <- sqrt(2 * maf * (1 - maf))  # ≈ 0.648 for p=0.3

  # Mean of standardised genotype: E[G] = 2p / sqrt(2*p*(1-p))
  mu_g <- 2 * maf / std_factor  # ≈ 0.9259 for p=0.3

  n_sample_sizes <- length(n_gwas)
  n_scenarios    <- length(scenarios)

  # ---- Scenario definitions ----
  # Each scenario = list(type, rho1, beta_g, gamma_e, beta_gxe)
  # type: "type1" or "power"
  # rho1: G→E correlation (0 for no mediation)
  # beta_g: genotype main effect coefficient
  # gamma_e: environment coefficient (unscaled; E can be sqrt(5)*E later)
  # beta_gxe: GxE interaction coefficient (= 0 for type I error)

  scenario_defs <- list(
    # Type I error scenarios
    # beta_g: variant 1 genotype coefficient
    # beta_g_null: variants 2-20 genotype coefficient (may differ for scenario 6)
    "1" = list(type = "type1", rho1 = 0,    beta_g = 0.1, beta_g_null = 0.1, env_scale = 1,   beta_gxe = 0,   label = "No mediation, E 1%"),
    "2" = list(type = "type1", rho1 = 0.01, beta_g = 0.2, beta_g_null = 0.2, env_scale = 1,   beta_gxe = 0,   label = "Mediation 1%, E 1%"),
    "3" = list(type = "type1", rho1 = 0.05, beta_g = 0.2, beta_g_null = 0.2, env_scale = sqrt(5), beta_gxe = 0, label = "Mediation 5%, E 5%"),
    # Power scenarios
    "4" = list(type = "power", rho1 = 0,    beta_g = 0.1, beta_g_null = 0.1, env_scale = 1,   beta_gxe = 0.1, label = "No mediation, E 1%"),
    "5" = list(type = "power", rho1 = 0.01, beta_g = 0.2, beta_g_null = 0.2, env_scale = 1,   beta_gxe = 0.2, label = "Mediation 1%, E 1%"),
    "6" = list(type = "power", rho1 = 0.05, beta_g = 0.1, beta_g_null = 0.2, env_scale = sqrt(5), beta_gxe = 0.2, label = "Mediation 5%, E 5%")
  )

  # ---- Storage: n_scenarios x n_sample_sizes x 3 tests ----
  results <- array(0, dim = c(n_scenarios, n_sample_sizes, 3))
  dimnames(results) <- list(
    scenario = sapply(scenario_defs[as.character(scenarios)], `[[`, "label"),
    n_gwas   = paste0("n=", n_gwas),
    test     = c("T_Direct", "T_MRGxE", "TwoStep")
  )

  # ---- Bonferroni threshold for 20 variants ----
  bonf_20 <- 0.05 / n_variants  # = 0.0025

  # ---- Main simulation ----
  for (s_idx in seq_along(scenarios)) {
    scen_num <- scenarios[s_idx]
    scen     <- scenario_defs[[as.character(scen_num)]]

    if (verbose) message(sprintf("\n=== Scenario %d: %s ===", scen_num, scen$label))

    for (i in seq_len(n_sample_sizes)) {
      n1 <- n_gwas[i]  # GWAS sample size
      n2 <- n_gwis      # GWIS sample size

      # Storage for this sample size: n_sim x 3 tests
      pv_store <- matrix(0, n_sim, 3)  # T_Direct, T_MRGxE, TwoStep
      theta_store <- numeric(n_sim)

      for (j in seq_len(n_sim)) {
        # Storage for 20 variants: 20 x 2 P-values (direct, mrgxe)
        pm <- matrix(0, n_variants, 2)

        for (k in seq_len(n_variants)) {
          # ---- Generate genotype for this variant ----
          g_raw <- stats::rbinom(n1, size = 2, prob = maf)
          g <- g_raw / std_factor  # standardised

          # ---- Generate environment ----
          # First variant may have mediation; variants 2-20 may not
          rho1 <- if (k == 1) scen$rho1 else 0

          e <- numeric(n1)

          # GWIS cohort (first n2 subjects): mean = mu_e
          e[seq_len(n2)] <- rho1 * g[seq_len(n2)] +
            stats::rnorm(n2, mean = mu_e, sd = sqrt(1 - rho1^2))

          # GWAS-only subjects (if any): mean = 2 * mu_e
          if (n1 > n2) {
            non_overlap_idx <- seq.int(n2 + 1, n1)
            e[non_overlap_idx] <- rho1 * g[non_overlap_idx] +
              stats::rnorm(length(non_overlap_idx), mean = 2 * mu_e,
                           sd = sqrt(1 - rho1^2))
          }

          mu_e1_obs <- mean(e[seq_len(n2)], na.rm = TRUE)

          # ---- Generate phenotype ----
          # Y = beta_g * G + gamma_e_unscaled * E + beta_gxe * G * E + N(0, 10)
          # For variant 1: use scen$beta_g; for variants 2-20: use scen$beta_g_null
          # (Scenario 6 has different beta_g for variant 1 vs null variants, per original code)
          bg <- if (k == 1) scen$beta_g else scen$beta_g_null
          y <- bg * g + scen$env_scale * e +
            scen$beta_gxe * g * e + stats::rnorm(n1, mean = 0, sd = sqrt(10))

          # ---- GWAS: y ~ g (marginal effect) ----
          fit1 <- stats::lm(y ~ g)
          s1   <- summary(fit1)$coefficients
          alpha_hat <- s1[2, 1]  # marginal effect
          alpha_se  <- s1[2, 2]

          # ---- GWIS: y2 ~ g2 * e2 (main + interaction effects) ----
          y2 <- y[seq_len(n2)]
          g2 <- g[seq_len(n2)]
          e2 <- e[seq_len(n2)]

          fit2 <- stats::lm(y2 ~ g2 * e2)
          s2   <- summary(fit2)$coefficients
          beta_hat <- s2[2, 1]   # main effect of G in GWIS
          beta_se  <- s2[2, 2]
          int_hat  <- s2[4, 1]   # GxE interaction
          int_se   <- s2[4, 2]

          # ---- T_Direct P-value ----
          pv_direct <- stats::pchisq((int_hat / int_se)^2, df = 1,
                                     lower.tail = FALSE)

          # ---- T_MRGxE ----
          # r = alpha_hat - beta_hat
          r_val <- alpha_hat - beta_hat

          # Rho (sample overlap correlation)
          # Formula from paper: rho = sqrt(n2/n1) / sqrt(1/(1-rho1^2) + (rho1*(1-mu_g)+mu_e)^2)
          rho <- sqrt(n2 / n1) / sqrt(1 / (1 - rho1^2) +
                                        (rho1 * (1 - mu_g) + mu_e1_obs)^2)

          # var(r) = var(alpha) + var(beta) - 2*rho*se(alpha)*se(beta)
          var_r <- alpha_se^2 + beta_se^2 - 2 * rho * alpha_se * beta_se
          var_r <- max(var_r, 1e-100)

          pv_mrgxe <- stats::pchisq(r_val^2 / var_r, df = 1,
                                    lower.tail = FALSE)

          # Store
          pm[k, 1] <- pv_direct
          pm[k, 2] <- pv_mrgxe

          # Store theta estimate for variant 1
          if (k == 1) theta_store[j] <- r_val
        }

        # ---- Power / Type I error for this replicate ----
        # T_Direct: min p-value across 20 variants < bonf_20
        # T_MRGxE:  min p-value across 20 variants < bonf_20
        # TwoStep:  first screen with T_MRGxE (bonf_20), then test
        #           survivors with T_Direct (Bonferroni corrected)

        pv_store[j, 1] <- as.numeric(min(pm[, 1]) < bonf_20)
        pv_store[j, 2] <- as.numeric(min(pm[, 2]) < bonf_20)

        # Two-step procedure
        survivors <- which(pm[, 2] < bonf_20)
        k_surv <- length(survivors)
        if (k_surv > 0) {
          if (k_surv == 1) {
            pv_store[j, 3] <- as.numeric(min(pm[survivors, 1]) < (0.05 / k_surv))
          } else {
            pv_store[j, 3] <- as.numeric(min(pm[survivors, 1]) < (0.05 / k_surv))
          }
        } else {
          pv_store[j, 3] <- 0
        }

        if (verbose && j %% 100 == 0) {
          message(sprintf("  Scenario %d, n=%d, rep %d/%d, power so far: %.3f, %.3f, %.3f",
                          scen_num, n1, j, n_sim,
                          mean(pv_store[1:j, 1]), mean(pv_store[1:j, 2]),
                          mean(pv_store[1:j, 3])))
        }
      }

      # Store average power/type I error
      results[s_idx, i, ] <- colMeans(pv_store)
    }
  }

  list(
    results        = results,
    scenario_labels = sapply(scenario_defs[as.character(scenarios)], `[[`, "label"),
    params = list(
      n_variants = n_variants,
      n_gwas     = n_gwas,
      n_gwis     = n_gwis,
      n_sim      = n_sim,
      maf        = maf,
      gamma      = gamma,
      mu_e       = mu_e,
      scenarios  = scenarios
    )
  )
}


# ===========================================================================
# Simplified wrapper: run a single mediation scenario
# ===========================================================================

#' Run a Single Mediation Scenario
#'
#' Convenience wrapper around \code{\link{simulate_gxe_fig2ef_exact}} that
#' runs one specific scenario, suitable for quick testing or vignettes.
#'
#' @param scenario  Integer 1-6 (see \code{\link{simulate_gxe_fig2ef_exact}}).
#' @param n_gwas    Integer vector. GWAS sample sizes.
#' @param n_gwis    Integer. GWIS sample size.
#' @param n_sim     Integer. Number of replicates.
#' @param ...       Additional arguments passed to
#'   \code{\link{simulate_gxe_fig2ef_exact}}.
#'
#' @return A data frame with columns: n_gwas, T_Direct, T_MRGxE, TwoStep.
#'
#' @export
run_single_scenario <- function(scenario, n_gwas = c(20000, 100000, 200000, 300000),
                                 n_gwis = 20000, n_sim = 1000, ...) {
  res <- simulate_gxe_fig2ef_exact(
    scenarios = scenario,
    n_gwas    = n_gwas,
    n_gwis    = n_gwis,
    n_sim     = n_sim,
    ...
  )

  data.frame(
    n_gwas    = n_gwas,
    T_Direct  = res$results[1, , 1],
    T_MRGxE   = res$results[1, , 2],
    TwoStep   = res$results[1, , 3],
    scenario  = res$scenario_labels[1],
    stringsAsFactors = FALSE
  )
}


# ===========================================================================
# Power curve plot function
# ===========================================================================

#' Plot Power Curves from Fig 2E-F Simulation
#'
#' Creates a publication-quality power/type-I error comparison plot
#' matching the style of Zhu et al. (2024) Fig 2E-F.
#'
#' @param sim_result  Result from \code{\link{simulate_gxe_fig2ef_exact}}.
#' @param type        \code{"type1"} for type I error or \code{"power"}
#'   for power plot.
#'
#' @return A ggplot object.
#'
#' @export
plot_fig2ef_power <- function(sim_result, type = c("type1", "power")) {
  type <- match.arg(type)
  y_label <- if (type == "type1") "Type I Error" else "Power"

  # Convert array to data frame
  n_scen <- dim(sim_result$results)[1]
  n_n    <- dim(sim_result$results)[2]

  df <- data.frame(
    scenario = rep(sim_result$scenario_labels, each = n_n * 3),
    n_gwas   = rep(rep(sim_result$params$n_gwas, each = 3), n_scen),
    test     = rep(c("T_Direct", "T_MRGxE", "TwoStep"), n_scen * n_n),
    value    = as.vector(sim_result$results),
    stringsAsFactors = FALSE
  )

  # Standard error for type I error (binomial proportion)
  se_05 <- 2 * sqrt(0.05 * 0.95 / sim_result$params$n_sim)

  # Filter: first 3 scenarios are type I, last 3 are power
  if (type == "type1") {
    df <- df[df$scenario %in% sim_result$scenario_labels[1:3], ]
    hline_val <- 0.05
    hline_label <- "0.05"
  } else {
    df <- df[df$scenario %in% sim_result$scenario_labels[4:6], ]
    hline_val <- NULL
    hline_label <- NULL
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting. Install with install.packages('ggplot2')")
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = factor(n_gwas), y = value,
                                         fill = test)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge",
                      color = "black", width = 0.7) +
    ggplot2::facet_grid(~ scenario) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid      = ggplot2::element_blank(),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0.02, 0.98),
      legend.title    = ggplot2::element_blank(),
      panel.border    = ggplot2::element_rect(size = 1.5),
      legend.text    = ggplot2::element_text(size = 10),
      axis.text      = ggplot2::element_text(size = 10),
      axis.title     = ggplot2::element_text(size = 12),
      strip.text     = ggplot2::element_text(size = 11)
    ) +
    ggplot2::geom_hline(yintercept = 0.05, linetype = "dashed",
                        color = "red", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = c("#45d9fd", "#e74c3c", "#2ecc71")) +
    ggplot2::labs(
      x = "GWAS sample size",
      y = y_label,
      fill = "Method"
    ) +
    ggplot2::ylim(0, max(1, max(df$value, na.rm = TRUE) * 1.1))

  p
}
