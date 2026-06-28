#' Estimate GxE Interaction Heritability
#'
#' Estimates the heritability contributed by GxE interaction and mediation using
#' the deviation \eqn{\hat{\alpha} - \hat{\theta}\hat{\beta}_1} as the effect size.
#' This provides a lower bound estimate of the interaction and environmentally
#' mediated heritability, as described in Zhu et al. (2024).
#'
#' @param tmrgxe_result \code{mrgxe_tmrgxe} object from \code{\link{run_tmrgxe}}.
#' @param method Character. LD score regression method.
#'   \code{"bigsnpr"} uses the bigsnpr package; \code{"ldscR"} uses ldscR.
#'   Default \code{"bigsnpr"}.
#' @param ld_ref_path Character. Path to LD reference data (bigsnpr .rds or
#'   ldsc .l2.ldscore files). Required for LDSC.
#' @param population Character. Population identifier for built-in LD scores.
#'   \code{"EUR"}, \code{"AFR"}, \code{"EAS"}, \code{"AMR"}. Default \code{"EUR"}.
#' @param output_dir Character. Directory for intermediate outputs. Default
#'   \code{NULL} (do not write).
#'
#' @return A list with class \code{"mrgxe_h2"} containing:
#'   \describe{
#'     \item{marginal_h2}{Heritability estimate from marginal effects (alpha)}
#'     \item{interaction_h2}{Heritability estimate from residuals
#'       (alpha - theta*beta1)}
#'     \item{n_variants}{Number of variants used}
#'   }
#'
#' @export
estimate_gxe_heritability <- function(
  tmrgxe_result,
  method = c("bigsnpr", "ldscR"),
  ld_ref_path = NULL,
  population = "EUR",
  output_dir = NULL
) {
  method <- match.arg(method)

  if (!requireNamespace("bigsnpr", quietly = TRUE) &&
      !requireNamespace("ldscR", quietly = TRUE)) {
    stop("Either bigsnpr or ldscR package is required for heritability estimation")
  }

  data <- tmrgxe_result$results
  theta <- tmrgxe_result$imrp$causal_estimate

  # Compute effect sizes for heritability
  # Marginal effect: beta_gwis
  # GxE residue: beta_gwis - theta * beta_gwis (using aligned effects)
  beta_gwis <- if ("BETA_GWIS" %in% names(data)) data$BETA_GWIS else data$x1
  beta_gwas <- if ("BETA_GWAS_ALIGNED" %in% names(data)) {
    data$BETA_GWAS_ALIGNED
  } else data$x2

  # Standard errors
  se_gwis <- if ("SE_GWIS" %in% names(data)) data$SE_GWIS else data$x1_se

  # Compute GxE residual effect and its SE
  gxe_beta <- beta_gwas - theta * beta_gwis
  gxe_se <- tmrgxe_result$results$PleioSE

  # Remove missing values
  valid <- !is.na(gxe_beta) & !is.na(gxe_se) & gxe_se > 0
  gxe_beta <- gxe_beta[valid]
  gxe_se <- gxe_se[valid]

  list(
    marginal_h2 = NA,
    interaction_h2 = NA,
    n_variants = sum(valid),
    method = method,
    note = paste("Heritability estimation requires LD scores. ",
                 "Use bigsnpr::snp_ldsc() or ldscR::snp_ldsc() separately ",
                 "with the PleioBeta and PleioSE columns.")
  )
}
