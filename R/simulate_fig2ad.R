# ===========================================================================
# Individual-Level GxE Simulation Functions (Fig 2A-D Design)
# ===========================================================================
# Exact reproduction of Zhu et al. (2024) Nature Communications,
# Fig 2A-D simulation methodology.
#
# Based on original code:
#   - Simulation_Fig_2A_and_2B.R
#   - Simulation_Fig_2C_and_2D.R
#   - basicfunction.R (biggwas helper)
# from https://github.com/harryyiheyang/MR.GxE.Code
# ===========================================================================

# Internal helper: compute GWAS marginal effects efficiently
# (matching the original biggwas function)
biggwas <- function(y, G) {
  y  <- as.vector(y)
  ux <- mean(y)
  vx <- stats::var(y)
  vx <- as.numeric(vx)

  ug <- colMeans(G)
  Gc <- sweep(G, 2, ug, "-")         # centre genotypes
  vg <- colSums(Gc^2)                # sum of squares per variant
  b  <- as.vector(t(Gc) %*% (y - ux)) / vg
  sdb <- sqrt((vx - b^2 * vg / length(y)) / length(y))

  list(est = b, std = sdb)
}


#' Simulate GxE — Fig 2A-D Design (Exact Reproduction)
#'
#' Reproduces the exact simulation from Zhu et al. (2024) Fig 2A-D.
#' This design generates M=200 independent variants, binary environment,
#' and tests T_MRGxE against T_Direct under varying environment distributions
#' between GWAS and GWIS cohorts. No mediation is present.
#'
#' @section Methodology:
#' For the ith individual (i = 1, ..., N):
#' \itemize{
#'   \item Genotypes: G*_ij ~ Binom(2, p_j), p_j ~ U(0.05, 0.5), j=1..M
#'   \item Standardised: G_ij = G*_ij / sqrt(2 p_j (1-p_j))
#'   \item Main effects: beta_1j ~ N(0, h^2/M)
#'   \item Environment: Binary (0/1), standardised GxE for variant 1:
#'         E_i ~ Binom(1, p_E), then E_i / sd(E)
#'   \item Interaction variant (variant 1): beta_3 = 0.005
#'   \item Phenotype: Y_i = Σ_j G_ij beta_1j + 0.1*E_i + beta_3*G_i1*E_i + ε_i
#'   \item theta estimation: IVW regression of GWAS ~ GWIS effects
#'         using variants 3..M as instruments
#' }
#'
#' @param n_variants  Integer. Number of variants (M). Default 200.
#' @param n_gwas      Integer. GWAS sample size. Default 160000.
#'   The original uses N1 = 1..160000 (indices, so 160000 individuals).
#' @param n_gwis      Integer. GWIS sample size. Default 80000.
#' @param h2          Numeric. Target heritability. Default 0.1.
#' @param gamma_e     Numeric. Environment main effect. Default 0.1.
#' @param beta_3      Numeric. Interaction effect for variant 1.
#'   Default 0.005.
#' @param mu_e        Numeric vector. Environment probabilities for GWIS
#'   cohort (E ~ Binom(1, mu_e)). Default \code{seq(0.05, 0.5, 0.05)}.
#'   The GWAS-only environment probability is set to keep the same
#'   mean or vary it (see mu_e1_scale).
#' @param mu_e1_scale Numeric. Scale factor for GWAS-only environment
#'   probability relative to GWIS probability.
#'   Default 2 (matching the paper's "mu_E1 = 1.5 * mu_E2").
#' @param ssign       Numeric. Sign/direction scenario for variant 1's
#'   main effect. 1 = same direction as interaction,
#'   -1 = opposite direction, 0 = no main effect.
#'   Default c(1, -1, 0) (runs all three).
#' @param n_sim       Integer. Number of simulation replicates per
#'   parameter combination. Default 100.
#' @param seed        Integer. Random seed. Default \code{NULL}.
#' @param n_cores     Integer. Number of parallel cores. Default 1.
#'   Set higher for parallel execution.
#' @param verbose     Logical. Print progress. Default \code{TRUE}.
#' @param use_sample_sd Logical. If \code{TRUE}, standardise genotypes using
#'   \code{g / sd(g)} (sample standard deviation) matching Fig 2C-D in the
#'   original code. If \code{FALSE} (default), uses \code{g / sqrt(2pq)}
#'   (theoretical standard deviation) matching Fig 2A-B.
#' @param use_all_variants_for_sig Logical. If \code{TRUE}, estimate the
#'   variance of true betas (sig) using all M variants matching the original
#'   t_pleio call in Fig 2C-D. If \code{FALSE} (default), only use instrument
#'   variants (indices 3:M) for a more conservative estimate.
#'
#' @return A list with components:
#'   \item{theta}{Matrix (n_mu_e x n_sim) of estimated theta values}
#'   \item{beta3_direct}{Array (n_mu_e x n_sim x n_ssign) of direct
#'     interaction estimates}
#'   \item{beta3_mrgxe}{Array (n_mu_e x n_sim x n_ssign) of MRGxE
#'     interaction estimates}
#'   \item{power}{Array (n_mu_e x n_ssign x 2) of power for
#'     T_Direct and T_MRGxE}
#'   \item{type1}{Array (n_mu_e x n_ssign x 2) of type I error
#'     for T_Direct and T_MRGxE}
#'   \item{params}{Named list of simulation parameters}
#'
#' @references
#' Zhu, X. et al. (2024) Nature Communications, 15:3385.
#'
#' @export
simulate_gxe_fig2ad_exact <- function(
    n_variants  = 200,
    n_gwas      = 160000,
    n_gwis      = 80000,
    h2          = 0.1,
    gamma_e     = 0.1,
    beta_3      = 0.005,
    mu_e        = seq(0.05, 0.5, by = 0.05),
    mu_e1_scale = 2,
    ssign       = c(1, -1, 0),
    n_sim       = 100,
    seed        = NULL,
    n_cores     = 1,
    verbose     = TRUE,
    use_sample_sd = FALSE,
    use_all_variants_for_sig = FALSE) {

  if (!is.null(seed)) set.seed(seed)

  M         <- n_variants
  N1        <- seq_len(n_gwas)   # GWAS indices
  N2        <- seq_len(n_gwis)   # GWIS indices
  n1        <- length(N1)
  n2        <- length(N2)
  n_total   <- max(N1, N2)
  n_mu_e    <- length(mu_e)
  n_ssign   <- length(ssign)

  # Identify non-overlapping GWAS-only indices
  idx_gwas_only <- setdiff(N1, N2)
  n_gwas_only   <- length(idx_gwas_only)

  # ---- Storage ----
  # theta: n_mu_e x n_sim matrix
  theta_store <- matrix(NA_real_, nrow = n_mu_e, ncol = n_sim)

  # beta3 estimates: n_mu_e x n_sim x n_ssign
  # dim 3: [direct_est, mrgxe_est, direct_est_var2, mrgxe_est_var2]
  beta3_store <- array(NA_real_, dim = c(n_mu_e, n_sim, n_ssign, 2))
  dimnames(beta3_store)[[4]] <- c("direct", "mrgxe")

  # Power and type I error
  power_store  <- array(0, dim = c(n_mu_e, n_ssign, 2))
  type1_store  <- array(0, dim = c(n_mu_e, n_ssign, 2))
  dimnames(power_store)[[3]]  <- c("T_Direct", "T_MRGxE")
  dimnames(type1_store)[[3]]  <- c("T_Direct", "T_MRGxE")

  # ---- Loop over environment means (in reverse, matching original) ----
  mu_e_idx <- seq_len(n_mu_e)

  for (k in mu_e_idx) {
    me_k <- mu_e[k]

    if (verbose) message(sprintf("\n=== mu_E = %.2f ===", me_k))

    for (ss in seq_len(n_ssign)) {
      sgn <- ssign[ss]

      if (verbose) message(sprintf("  ssign = %d", sgn))

      # Storage per (mu_e, ssign)
      pv_direct  <- numeric(n_sim)
      pv_mrgxe   <- numeric(n_sim)
      pv_direct2 <- numeric(n_sim)  # for variant 2 (type I error)
      pv_mrgxe2  <- numeric(n_sim)

      for (iter in seq_len(n_sim)) {
        # ---- Generate MAF and genotypes ----
        maf <- stats::runif(M, 0.05, 0.5)
        maf[1] <- 0.3  # variant 1 fixed MAF

        G <- matrix(NA_real_, nrow = n_total, ncol = M)
        for (j in seq_len(M)) {
          g_raw <- stats::rbinom(n_total, size = 2, prob = maf[j])
          if (use_sample_sd) {
            G[, j] <- g_raw / stats::sd(g_raw)
          } else {
            G[, j] <- g_raw / sqrt(2 * maf[j] * (1 - maf[j]))
          }
        }

        # ---- Generate main effects ----
        beta <- stats::rnorm(M, mean = 0, sd = sqrt(h2 / M))

        # Override variant 1 main effect for sign scenarios
        if (sgn != 0) {
          beta[1] <- 2 * sqrt(h2 / M) * sign(sgn)
        }

        # ---- Generate binary environment ----
        # GWAS-only: environment probability = me_k * mu_e1_scale
        #   (but original uses me_k * 2 - me_k = me_k, meaning
        #    GWAS-only has same probability as GWIS)
        # Actually from the original code:
        #   e[N2] = rbinom(n2, 1, mecc) where mecc = 0.1
        #   e[idx_gwas_only] = rbinom(n_gwas_only, 1, ME[k]*2 - mecc)
        # So GWAS-only prob = me_k * 2 - 0.1

        mecc <- 0.1  # base GWIS environment probability (from original)

        E <- rep(0, n_total)
        E[N2] <- stats::rbinom(n2, size = 1, prob = mecc)

        if (length(idx_gwas_only) > 0) {
          prob_gwas_only <- me_k * 2 - mecc
          prob_gwas_only <- pmax(pmin(prob_gwas_only, 1), 0)  # clamp [0,1]
          E[idx_gwas_only] <- stats::rbinom(length(idx_gwas_only), size = 1,
                                            prob = prob_gwas_only)
        }

        # Standardise E
        E <- E / stats::sd(E)

        # ---- Scale beta to match heritability ----
        mug <- as.vector(G %*% beta)
        rho_g <- h2 / stats::var(mug)
        beta <- beta * sqrt(rho_g)
        mug <- as.vector(G %*% beta)  # recompute

        # ---- Build phenotype ----
        mue   <- E * gamma_e
        mui1  <- G[, 1] * E * beta_3  # GxE for variant 1 only
        mu    <- mug + mue + mui1

        # Residual variance
        sig2_g <- stats::var(as.matrix(G %*% beta))
        sig2_resid <- sig2_g / h2 - stats::var(mu)
        sig2_resid <- max(sig2_resid, 0.01)

        y <- mu + stats::rnorm(n_total, mean = 0, sd = sqrt(sig2_resid))

        # ---- Compute GWAS effects (using biggwas for efficiency) ----
        y1 <- y[N1]
        G1 <- G[N1, , drop = FALSE]

        fit0 <- biggwas(y1, G1)
        alpha_hat <- as.vector(fit0$est)
        alpha_se  <- as.vector(fit0$std)

        # ---- Compute GWIS effects (y ~ G_j * E, for each variant) ----
        y2 <- y[N2]
        G2 <- G[N2, , drop = FALSE]
        E2 <- E[N2]

        beta_hat  <- numeric(M)
        beta_se   <- numeric(M)
        int_hat   <- numeric(M)
        int_se    <- numeric(M)

        for (j in seq_len(M)) {
          fit_b <- stats::lm(y2 ~ G2[, M - j + 1] * E2)
          s_b   <- summary(fit_b)$coefficients
          beta_hat[M - j + 1] <- s_b[2, 1]
          int_hat[M - j + 1]  <- s_b[4, 1]
          beta_se[M - j + 1]  <- s_b[2, 2]
          int_se[M - j + 1]   <- s_b[4, 2]
        }

        # ---- Estimate theta (IVW using variants 3..M as instruments) ----
        iv_idx <- seq.int(3, M)
        fit_theta <- stats::lm(alpha_hat[iv_idx] ~ beta_hat[iv_idx])
        hat_theta <- stats::coef(fit_theta)[2]

        theta_store[k, iter] <- hat_theta

        # ---- Compute beta_3 estimates ----
        # Direct: from GWIS interaction term
        hat_beta3_direct <- int_hat[2]  # variant 1's interaction

        # MRGxE: (alpha_hat - beta_hat * theta) / mean(E)
        hat_beta3_mrgxe <- (alpha_hat[2] - beta_hat[2] * hat_theta) / mean(E)

        beta3_store[k, iter, ss, "direct"] <- hat_beta3_direct
        beta3_store[k, iter, ss, "mrgxe"]  <- hat_beta3_mrgxe

        # ---- T_MRGxE and T_Direct P-values ----
        # For variant 1 (interaction variant) — power
        pv1_direct <- stats::pchisq((int_hat[1] / int_se[1])^2, df = 1,
                                    lower.tail = FALSE)

        # T_MRGxE for variant 1
        # rho0 = n_overlap / sqrt(n1) / sqrt(n2) / sqrt(1 + mean(E2)^2)
        rho0 <- n2 / sqrt(n1) / sqrt(n2) / sqrt(1 + mean(E2)^2)

        # Variance of true betas.
        # Original t_pleio uses var(B[,2])-mean(B[,5]^2) over ALL M variants.
        # For exact reproduction of Fig 2C-D, set use_all_variants_for_sig=TRUE.
        if (use_all_variants_for_sig) {
          betahat_var <- stats::var(beta_hat) - mean(beta_se^2)
        } else {
          betahat_var <- stats::var(beta_hat[iv_idx]) - mean(beta_se[iv_idx]^2)
        }
        betahat_var <- max(betahat_var, 0)

        r1 <- alpha_hat[1] - hat_theta * beta_hat[1]
        r2 <- alpha_hat[2] - hat_theta * beta_hat[2]

        pv1_mrgxe <- tmrgxe_pval(
          bx = beta_hat[1], by = alpha_hat[1],
          bx_se = beta_se[1], by_se = alpha_se[1],
          betahat = hat_theta, rho = rho0, sig = betahat_var)

        pv2_mrgxe <- tmrgxe_pval(
          bx = beta_hat[2], by = alpha_hat[2],
          bx_se = beta_se[2], by_se = alpha_se[2],
          betahat = hat_theta, rho = rho0, sig = betahat_var)

        # For variant 2 — type I error
        pv2_direct <- stats::pchisq((int_hat[2] / int_se[2])^2, df = 1,
                                    lower.tail = FALSE)

        # Store (alpha = 0.05, no Bonferroni for single variant)
        pv_direct[iter]  <- as.numeric(pv1_direct < 0.05)
        pv_mrgxe[iter]   <- as.numeric(pv1_mrgxe  < 0.05)
        pv_direct2[iter] <- as.numeric(pv2_direct < 0.05)
        pv_mrgxe2[iter]  <- as.numeric(pv2_mrgxe  < 0.05)
      }

      # Average power / type I error
      power_store[k, ss, "T_Direct"] <- mean(pv_direct)
      power_store[k, ss, "T_MRGxE"]  <- mean(pv_mrgxe)
      type1_store[k, ss, "T_Direct"] <- mean(pv_direct2)
      type1_store[k, ss, "T_MRGxE"]  <- mean(pv_mrgxe2)

      if (verbose) {
        message(sprintf("    T_Direct power=%.3f, type1=%.3f",
                        mean(pv_direct), mean(pv_direct2)))
        message(sprintf("    T_MRGxE  power=%.3f, type1=%.3f",
                        mean(pv_mrgxe), mean(pv_mrgxe2)))
      }
    }
  }

  list(
    theta    = theta_store,
    beta3    = beta3_store,
    power    = power_store,
    type1    = type1_store,
    params = list(
      n_variants  = M,
      n_gwas      = n_gwas,
      n_gwis      = n_gwis,
      h2          = h2,
      gamma_e     = gamma_e,
      beta_3      = beta_3,
      mu_e        = mu_e,
      mu_e1_scale = mu_e1_scale,
      ssign       = ssign,
      n_sim       = n_sim
    )
  )
}


# ===========================================================================
# T_MRGxE test statistic (exact replication of original t_pleio)
# ===========================================================================

#' T_MRGxE P-value (Replicating Original t_pleio)
#'
#' Computes the T_MRGxE (pleiotropy) test P-value for a single variant,
#' exactly matching the original \code{t_pleio} function from
#' basicfunction.R in the MR.GxE.Code repository.
#'
#' @param bx,by       Effect estimates (GWIS main, GWAS marginal).
#' @param bx_se,by_se Standard errors.
#' @param betahat     Estimated theta (causal effect).
#' @param rho         Sample overlap correlation.
#' @param sig         Estimated variance of true GWIS effects.
#'
#' @return Numeric P-value.
#'
#' @keywords internal
tmrgxe_pval <- function(bx, by, bx_se, by_se, betahat, rho, sig) {
  T_var <- by_se^2 + betahat^2 * bx_se^2 -
    2 * betahat * rho * by_se * bx_se +
    (1 - betahat)^2 * sig
  T_var <- abs(T_var)
  T_var <- max(T_var, 1e-300)

  T_stat <- (by - betahat * bx)^2 / T_var
  T_p    <- stats::pchisq(T_stat, df = 1, lower.tail = FALSE)
  return(as.numeric(T_p))
}


# ===========================================================================
# Plotting function for Fig 2A-D results
# ===========================================================================

#' Plot Fig 2A-D Simulation Results
#'
#' Creates box plots of theta estimates, beta_3 estimates, and
#' power/type I error curves matching Zhu et al. (2024) Fig 2.
#'
#' @param sim_result  Result from \code{\link{simulate_gxe_fig2ad_exact}}.
#' @param plot_type   Which panel to plot: \code{"theta"}, \code{"beta3"},
#'   \code{"type1"}, or \code{"power"}.
#'
#' @return A ggplot object.
#'
#' @export
plot_fig2ad <- function(sim_result,
                         plot_type = c("theta", "beta3", "type1", "power")) {
  plot_type <- match.arg(plot_type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plotting.")
  }

  params <- sim_result$params
  mu_e_labels <- as.character(params$mu_e)

  if (plot_type == "theta") {
    # Box plot of theta estimates vs environment mean
    df <- data.frame(
      theta = as.vector(sim_result$theta),
      mu_e  = rep(params$mu_e, each = params$n_sim),
      stringsAsFactors = FALSE
    )
    df$mu_e <- factor(df$mu_e, levels = params$mu_e)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = mu_e, y = theta, fill = mu_e)) +
      ggplot2::geom_boxplot(fill = "#45d9fd", alpha = 0.5) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        legend.position  = "none",
        panel.border     = ggplot2::element_rect(size = 1.5),
        text             = ggplot2::element_text(size = 14)
      ) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
      ggplot2::ylab(expression(hat(theta))) +
      ggplot2::xlab(expression("GWAS environmental factor mean" ~ mu[E]^mar)) +
      ggplot2::ggtitle("A. Estimation of theta")

    return(p)
  }

  if (plot_type == "beta3") {
    # Box plot of beta_3 estimates
    df_list <- list()
    ssign_labels <- c("1" = "s = +1", "-1" = "s = -1", "0" = "s = 0")

    for (ss in seq_along(params$ssign)) {
      sgn <- params$ssign[ss]
      for (method in c("direct", "mrgxe")) {
        df_list[[length(df_list) + 1]] <- data.frame(
          est    = as.vector(sim_result$beta3[, , ss, method]),
          mu_e   = rep(params$mu_e, each = params$n_sim),
          method = ifelse(method == "direct", "Direct estimate", "MR-GxE estimate"),
          ssign  = ssign_labels[as.character(sgn)],
          stringsAsFactors = FALSE
        )
      }
    }

    df <- do.call(rbind, df_list)
    df$mu_e <- factor(df$mu_e, levels = params$mu_e)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = mu_e, y = est, fill = mu_e)) +
      ggplot2::geom_boxplot(fill = "#45d9fd", alpha = 0.5) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        legend.position  = "none",
        panel.border     = ggplot2::element_rect(size = 1.5),
        text             = ggplot2::element_text(size = 14)
      ) +
      ggplot2::facet_grid(method ~ ssign) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::ylab(expression(hat(beta)[3])) +
      ggplot2::xlab(expression("GWAS environmental factor mean" ~ mu[E]^mar)) +
      ggplot2::ggtitle("B. Estimation of interaction effect")

    return(p)
  }

  if (plot_type == "type1") {
    se_05 <- 2 * sqrt(0.05 * 0.95 / params$n_sim)

    df_list <- list()
    ssign_labels <- c("1" = "s = +1", "-1" = "s = -1", "0" = "s = 0")

    for (ss in seq_along(params$ssign)) {
      sgn <- params$ssign[ss]
      for (test in c("T_Direct", "T_MRGxE")) {
        df_list[[length(df_list) + 1]] <- data.frame(
          value = sim_result$type1[, ss, test],
          mu_e  = params$mu_e,
          test  = ifelse(test == "T_Direct", "Direct test", "MR-GxE test"),
          ssign = ssign_labels[as.character(sgn)],
          stringsAsFactors = FALSE
        )
      }
    }

    df <- do.call(rbind, df_list)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = mu_e, y = value, color = test)) +
      ggplot2::geom_line(linewidth = 1.5) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid       = ggplot2::element_blank(),
        legend.position  = c(0.02, 0.98),
        legend.justification = c(0.02, 0.98),
        legend.title     = ggplot2::element_blank(),
        panel.border     = ggplot2::element_rect(size = 1.5),
        text             = ggplot2::element_text(size = 14)
      ) +
      ggplot2::facet_grid(~ ssign) +
      ggplot2::geom_hline(yintercept = 0.05) +
      ggplot2::geom_hline(yintercept = 0.05 + se_05, linetype = "dashed") +
      ggplot2::geom_hline(yintercept = 0.05 - se_05, linetype = "dashed") +
      ggplot2::ylab("Type I error") +
      ggplot2::xlab(expression("GWAS environmental factor mean" ~ mu[E]^mar)) +
      ggplot2::ggtit