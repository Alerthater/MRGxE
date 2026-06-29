#' Two-Step Procedure for GxE Validation
#'
#' Implements the two-step validation procedure described in Zhu et al. (2024):
#' Step 1 screens with TMRGxE, Step 2 validates significant variants using the
#' direct interaction test (TDirect) from the GWIS data. This substantially
#' reduces the multiple testing burden.
#'
#' @param tmrgxe_result Object of class \code{mrgxe_tmrgxe} from
#'   \code{\link{run_tmrgxe}}.
#' @param gwis_data Data frame. Original GWIS data containing columns needed for
#'   the direct interaction test (interaction effect and its SE if available).
#'   If \code{NULL}, the function will attempt to extract from the TMRGxE result.
#' @param step1_threshold Numeric. TMRGxE significance threshold for step 1.
#'   Default \code{5e-8}.
#' @param step2_alpha Numeric. Overall alpha level for step 2 Bonferroni
#'   correction. Default 0.05.
#' @param interaction_beta_col Character. Column name for the GxE interaction
#'   effect from GWIS. If available, enables TDirect calculation. Default
#'   \code{NULL} (not required).
#' @param interaction_se_col Character. Column name for the interaction effect SE.
#'   Default \code{NULL}.
#'
#' @return A list with class \code{"mrgxe_twostep"} containing:
#'   \describe{
#'     \item{step1_significant}{Data frame of variants significant in step 1}
#'     \item{step2_results}{Data frame with step 2 TDirect results}
#'     \item{step2_significant}{Data frame of variants passing step 2}
#'     \item{step2_threshold}{Bonferroni-corrected threshold for step 2}
#'     \item{params}{Parameters used}
#'   }
#'
#' @details
#' When direct interaction data is not available, step 2 results will not be
#' computed and the function returns only step 1 findings with a note.
#'
#' @examples
#' \dontrun{
#' validation <- two_step_validation(screen)
#' print(validation$step2_significant)
#' }
#'
#' @export
two_step_validation <- function(
  tmrgxe_result,
  gwis_data = NULL,
  step1_threshold = 5e-8,
  step2_alpha = 0.05,
  interaction_beta_col = NULL,
  interaction_se_col = NULL
) {
  if (!inherits(tmrgxe_result, "mrgxe_screen") &&
      !inherits(tmrgxe_result, "mrgxe_tmrgxe")) {
    stop("tmrgxe_result must be an object of class 'mrgxe_screen' or 'mrgxe_tmrgxe'")
  }

  results <- tmrgxe_result$results

  # Determine which P column to use
  # screen_interaction() uses "InteractionP"; legacy used "PleioP_MR"
  p_col <- if ("InteractionP" %in% names(results)) "InteractionP" else
           if ("PleioP_MR" %in% names(results)) "PleioP_MR" else
           stop("Cannot find interaction P-value column in results")

  # Step 1: Identify significant variants
  step1_sig <- results[results[[p_col]] < step1_threshold &
                         !is.na(results[[p_col]]), ]

  if (nrow(step1_sig) == 0) {
    message("No variants passed step 1 (TMRGxE) screening")
    return(list(
      step1_significant = step1_sig,
      step2_results = NULL,
      step2_significant = NULL,
      step2_threshold = NULL,
      message = "No significant variants at step 1"
    ))
  }

  message(sprintf("Step 1: %d significant variants (P < %.0e)",
                  nrow(step1_sig), step1_threshold))

  # Step 2: Direct test
  # Bonferroni correction for step 1 significant variants
  x <- nrow(step1_sig)
  step2_threshold <- step2_alpha / max(x, 1)

  step2_results <- NULL
  step2_sig <- NULL

  if (!is.null(interaction_beta_col) && !is.null(interaction_se_col)) {
    # TDirect from GWIS interaction terms
    if (interaction_beta_col %in% names(step1_sig) &&
        interaction_se_col %in% names(step1_sig)) {
      tdirect_z <- step1_sig[[interaction_beta_col]] /
        step1_sig[[interaction_se_col]]
      tdirect_p <- 2 * stats::pnorm(-abs(tdirect_z))

      step2_results <- data.frame(
        SNP = step1_sig$SNP,
        CHR = step1_sig$CHR,
        BP = step1_sig$BP,
        TDirect_Z = tdirect_z,
        TDirect_P = tdirect_p,
        TDirect_sig = tdirect_p < step2_threshold,
        stringsAsFactors = FALSE
      )

      step2_sig <- step2_results[step2_results$TDirect_sig, ]
      message(sprintf("Step 2 threshold: %.2e (correcting for %d tests)",
                      step2_threshold, x))
      message(sprintf("Step 2: %d significant (P < %.2e)",
                      nrow(step2_sig), step2_threshold))
    } else {
      message("Interaction columns not found in results; skipping step 2 TDirect calculation")
    }
  } else {
    message("No interaction columns provided; step 2 requires TDirect data.")
    message("Set interaction_beta_col and interaction_se_col pointing to the GWIS interaction effect.")
  }

  result <- list(
    step1_significant = step1_sig,
    step2_results = step2_results,
    step2_significant = step2_sig,
    step2_threshold = step2_threshold,
    params = list(
      step1_threshold = step1_threshold,
      step2_alpha = step2_alpha
    )
  )
  class(result) <- "mrgxe_twostep"
  result
}

#' @export
pr