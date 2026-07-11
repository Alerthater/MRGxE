#' Estimate GxE Interaction Heritability
#'
#' Estimates the heritability contributed by GxE interaction and mediation using
#' the deviation \eqn{\hat{\alpha} - \hat{\theta}\hat{\beta}_1} as the effect size.
#' This provides a lower bound estimate of the interaction and environmentally
#' mediated heritability, as described in Zhu et al. (2024).
#'
#' @param tmrgxe_result \code{mrgxe_screen} object from \code{\link{screen_interaction}}
#'   or \code{\link{run_tmrgxe}}.
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

  beta_gwis <- if ("BETA_GWIS" %in% names(data)) data$BETA_GWIS else data$x1
  beta_gwas <- if ("BETA_GWAS_ALIGNED" %in% names(data)) {
    data$BETA_GWAS_ALIGNED
  } else data$x2

  gxe_beta <- beta_gwas - theta * beta_gwis
  gxe_se   <- tmrgxe_result$results$InteractionSE

  valid <- !is.na(gxe_beta) & !is.na(gxe_se) & gxe_se > 0

  msg <- paste("Heritability estimation requires LD scores.",
               "Use bigsnpr::snp_ldsc() or ldscR::snp_ldsc() separately",
               "with the InteractionBeta and InteractionSE columns.")

  list(
    marginal_h2    = NA,
    interaction_h2 = NA,
    n_variants     = sum(valid),
    method         = method,
    note           = msg
  )
}

#' Estimate Interaction Heritability (Alias)
#'
#' Alias for \code{\link{estimate_gxe_heritability}}.
#'
#' @inheritParams estimate_gxe_heritability
#' @return Same as \code{\link{estimate_gxe_heritability}}.
#' @export
estimate_interaction_heritability <- function(
  tmrgxe_result,
  method = c("bigsnpr", "ldscR"),
  ld_ref_path = NULL,
  population = "EUR",
  output_dir = NULL
) {
  estimate_gxe_heritability(
    tmrgxe_result = tmrgxe_result,
    method        = method,
    ld_ref_path   = ld_ref_path,
    population    = population,
    output_dir    = output_dir
  )
}
