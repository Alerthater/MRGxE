#' Run IMRP to Estimate Theta
#'
#' Estimates the causal/calibration parameter theta using IMRP, which serves as
#' the regression slope between effect B (exposure) and effect A (outcome).
#'
#' @param data Data frame (from \code{harmonize_effects}) containing columns
#'   \code{x1} (effect B) and \code{x2} (effect A).
#' @param instruments Data frame of selected instrumental variables.
#' @param rho Numeric. Correlation coefficient (from \code{\link{estimate_rho}}).
#'   May be 0 for designs where sample overlap is not relevant.
#' @param method Character. MR method: \code{"IVW"} (default) or \code{"MR-Egger"}.
#' @param signif_threshold Numeric. Threshold for pleiotropic outlier detection. Default 0.05.
#' @param design Optional \code{mrgxe_design} object.
#' @param ... Additional arguments.
#'
#' @return An object of class \code{"mrgxe_imrp"} containing theta estimate and
#'   global pleiotropy statistics.
#'
#' @export
run_imrp <- function(
  data,
  instruments,
  rho,
  method = "IVW",
  signif_threshold = 0.05,
  design = NULL,
  ...
) {
  if (!requireNamespace("IMRP", quietly = TRUE)) {
    stop("The IMRP package is required. Install with:\n",
         "  remotes::install_github('XiaofengZhuCase/IMRP')")
  }

  # Determine column names
  x1_col <- if ("x1" %in% names(instruments)) "x1" else
    if ("beta_b" %in% names(instruments)) {
      # Convert to standardized if needed
      "beta_b"
    } else stop("Cannot find exposure column in instruments")

  x2_col <- if ("x2" %in% names(instruments)) "x2" else
    if ("beta_a" %in% names(instruments)) "beta_a"
    else stop("Cannot find outcome column in instruments")

  se1_col <- if ("x1_se" %in% names(instruments)) "x1_se" else
    if ("se_b" %in% names(instruments)) "se_b"
    else stop("Cannot find exposure SE column")

  se2_col <- if ("x2_se" %in% names(instruments)) "x2_se" else
    if ("se_a" %in% names(instruments)) "se_a"
    else stop("Cannot find outcome SE column")

  iv_df <- as.data.frame(instruments)
  message(sprintf("Running IMRP with %d instruments, rho = %.4f",
                  nrow(iv_df), rho))

  mr <- IMRP::MR_pleio(
    BetaOutcome = x2_col,
    BetaExposure = x1_col,
    SdOutcome = se2_col,
    SdExposure = se1_col,
    data = iv_df,
    SignifThreshold = signif_threshold,
    rho = rho,
    method = method
  )

  causal_beta <- as.numeric(mr$CausalEstimate)
  causal_se <- as.numeric(mr$SdCausalEstimate)
  causal_p <- as.numeric(mr$Causal_p)

  message(sprintf("Theta = %.4f (SE = %.4f, P = %.2e)",
                  causal_beta, causal_se, causal_p))

  result <- list(
    mr_result = mr,
    causal_estimate = causal_beta,
    causal_se = causal_se,
    causal_p = causal_p,
    global_p_pre = as.numeric(mr$GlobalPvalue_pre),
    global_p_aft = as.numeric(mr$GlobalPvalue_aft),
    instruments_used = iv_df,
    rho = rho,
    params = list(method = method, signif_threshold = signif_threshold)
  )
  class(result) <- "mrgxe_imrp"
  result
}

#' Search Pleiotropy with Effect Estimates
#'
#' Computes genome-wide interaction deviation statistics. Generalizes
#' \code{search_pleio_with_effect} to any pair of effect estimates.
#'
#' @param data Data frame with \code{x1}, \code{x2}, \code{x1_se}, \code{x2_se},
#'   or generic columns.
#' @param rho Numeric. Correlation between the two effect estimates.
#' @param causal_estimate Numeric. Theta from \code{\link{run_imrp}}.
#' @param causal_se Numeric. Standard error of theta.
#' @param x1_col Character. Column name for exposure effect.
#' @param x2_col Character. Column name for outcome effect.
#' @param se1_col Character. Column name for exposure SE.
#' @param se2_col Character. Column name for outcome SE.
#'
#' @return A list with \code{pleio_p}, \code{beta} (deviation), \code{se}, \code{z}.
#'
#' @export
search_pleio_with_effect <- function(
  data,
  rho,
  causal_estimate,
  causal_se,
  x1_col = "x1",
  x2_col = "x2",
  se1_col = "x1_se",
  se2_col = "x2_se"
) {
  required <- c(x1_col, x2_col, se1_col, se2_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns in data: ", paste(missing, collapse = ", "))
  }

  beta_out <- data[[x2_col]]
  beta_exp <- data[[x1_col]]
  se_out <- data[[se2_col]]
  se_exp <- data[[se1_col]]

  deviation <- beta_out - causal_estimate * beta_exp

  deviation_se <- sqrt(
    se_out^2 +
      causal_estimate^2 * se_exp^2 +
      beta_exp^2 * causal_se^2 -
      2 * causal_estimate * rho * se_exp * se_out
  )

  deviation_se <- pmax(deviation_se, .Machine$double.eps)
  pval <- 2 * stats::pnorm(-abs(deviation / deviation_se))
  pval <- pmin(pval, 1)

  list(
    pleio_p = pval,
    beta = deviation,
    se = deviation_se,
    z = deviation / deviation_se
  )
}

#' Genome-wide Interaction Screening
#'
#' The main screening function: runs IMRP on instruments to estimate theta,
#' then screens all variants for deviation from the regression line.
#' Works for any study design (GxE, genome-to-genome, custom).
#'
#' @param data Data frame of harmonized effect estimates.
#' @param instruments Data frame of selected instruments.
#' @param rho Numeric. Correlation coefficient.
#' @param standardize Logical. Whether effects were standardized.
#'   Default \code{TRUE} for GxE, \code{FALSE} for genome-to-genome.
#' @param method Character. MR method for IMRP.
#' @param signif_threshold Numeric. Significance threshold for interaction.
#' @param design Optional \code{mrgxe_design} object (for labels/reporting).
#' @param ... Additional arguments to \code{\link{run_imrp}}.
#'
#' @return An object of class \code{"mrgxe_screen"} with:
#'   \describe{
#'     \item{results}{Full data with interaction statistics appended}
#'     \item{imrp}{IMRP result}
#'     \item{n_significant}{Count of significant variants}
#'     \item{threshold}{Threshold used}
#'     \item{design}{Study design used}
#'   }
#'
#' @export
screen_interaction <- function(
  data,
  instruments,
  rho,
  standardize = TRUE,
  method = "IVW",
  signif_threshold = 5e-8,
  design = NULL,
  ...
) {
  # Tag with design
  if (is.null(design)) design <- study_design("custom")

  # Prepare standardized analysis columns
  dat <- .prepare_analysis_data(data, standardize = standardize)

  # Standardize instruments similarly
  iv_dat <- .prepare_analysis_data(instruments, standardize = standardize)

  # Run IMRP
  imrp <- run_imrp(
    data = dat,
    instruments = iv_dat,
    rho = rho,
    method = method,
    standardize = FALSE,  # already done
    design = design,
    ...
  )

  # Screen all variants
  pleio <- search_pleio_with_effect(
    data = dat,
    rho = rho,
    causal_estimate = imrp$causal_estimate,
    causal_se = imrp$causal_se
  )

  # Attach results
  dat$InteractionBeta <- as.numeric(pleio$beta)
  dat$InteractionSE <- as.numeric(pleio$se)
  dat$InteractionZ <- as.numeric(pleio$z)
  dat$InteractionP <- as.numeric(pleio$pleio_p)
  dat$is_significant <- dat$InteractionP < signif_threshold

  sig_count <- sum(dat$is_significant, na.rm = TRUE)
  if (sig_count == 0) {
    message("No significant variants at P < ", format(signif_threshold, scientific = TRUE))
  } else {
    message(sprintf("Found %d significant variants at %.0e", sig_count, signif_threshold))
  }

  result <- list(
    results = dat,
    imrp = imrp,
    n_significant = sig_count,
    significant_variants = dat[which(dat$is_significant), ],
    threshold = signif_threshold,
    design = design,
    params = list(
      standardize = standardize,
      method = method,
      signif_threshold = signif_threshold,
      rho = rho
    )
  )
  class(result) <- "mrgxe_screen"
  result
}

# ---------------------------------------------------------------------------
# Legacy TMRGxE wrapper (backward compatible)
# ---------------------------------------------------------------------------

#' @rdname screen_interaction
#' @export
run_tmrgxe <- function(
  data,
  instruments,
  rho,
  standardize = TRUE,
  method = "IVW",
  signif_threshold = 5e-8,
  ...
) {
  design <- study_design("gxe")
  screen_interaction(
    data = data,
    instruments = instruments,
    rho = rho,
    standardize = standardize,
    method = method,
    signif_threshold = signif_threshold,
    design = design,
    ...
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.prepare_analysis_data <- function(data, standardize = TRUE) {
  if (standardize) {
    # Check if already standardized (has x1/x2)
    if (all(c("x1", "x2", "x1_se", "x2_se") %in% names(data))) {
      return(data)
    }

    # Try to standardize from raw columns
    if (all(c("beta_a", "beta_b", "se_a", "se_b") %in% names(data))) {
      n_a <- if ("n_a" %in% names(data)) data$n_a else 1
      n_b <- if ("n_b" %in% names(data)) data$n_b else 1
      data$x1 <- data$beta_b / data$se_b / sqrt(n_b)
      data$x2 <- data$beta_a / data$se_a / sqrt(n_a)
      data$x1_se <- 1 / sqrt(n_b)
      data$x2_se <- 1 / sqrt(n_a)
      return(data)
    }

    # Legacy GxE format
    if (all(c("BETA_GWIS", "BETA_GWAS_ALIGNED", "SE_GWIS", "SE_GWAS") %in% names(data))) {
      n_gwis <- if ("N_GWIS" %in% names(data)) data$N_GWIS else 5000
      n_gwas <- if ("N_GWAS" %in% names(data)) data$N_GWAS else 5000
      data$x1 <- data$BETA_GWIS / data$SE_GWIS / sqrt(n_gwis)
      data$x2 <- data$BETA_GWAS_ALIGNED / data$SE_GWAS / sqrt(n_gwas)
      data$x1_se <- 1 / sqrt(n_gwis)
      data$x2_se <- 1 / sqrt(n_gwas)
      return(data)
    }

    stop("Cannot standardize: required columns not found")
  }

  # No standardization: use raw or existing columns
  if (all(c("x1", "x2", "x1_se", "x2_se") %in% names(data))) return(data)
  if (all(c("beta_b", "beta_a", "se_b", "se_a") %in% names(data))) {
    data$x1 <- data$beta_b
    data$x2 <- data$beta_a
    data$x1_se <- data$se_b
    data$x2_se <- data$se_a
    return(data)
  }
  if (all(c("BETA_GWIS", "BETA_GWAS_ALIGNED", "SE_GWIS", "SE_GWAS") %in% names(data))) {
    data$x1 <- data$BETA_GWIS
    data$x2 <- data$BETA_GWAS_ALIGNED
    data$x1_se <- data$SE_GWIS
    data$x2_se <- data$SE_GWAS
    return(data)
  }

  stop("Cannot prepare analysis data: unrecognized column structure")
}

#' @export
print.mrgxe_imrp <- function(x, ...) {
  cat(sprintf("IMRP: theta = %.4f (SE = %.4f, P = %.2e)\n",
              x$causal_estimate, x$causal_se, x$causal_p))
  cat(sprintf("  Global P (pre/aft): %.2e / %.2e\n",
              x$global_p_pre, x$global_p_aft))
  cat(sprintf("  Instruments: %d\n", nrow(x$instruments_used)))
}

#' @export
print.mrgxe_screen <- function(x, ...) {
  cat("Interaction screening result\n")
  cat(sprintf("  Variants screened: %d\n", nrow(x$results)))
  cat(sprintf("  Significant:       %d\n", x$n_significant))
  if (!is.null(x$design)) cat(sprintf("  Design: %s\n", x$design$design_type))
}
