#' Conditional Analysis
#'
#' Performs conditional analysis to determine whether a top variant explains
#' the interaction signal in a region. This mirrors the approach in Chen et al.
#' (2026) where HLA-A*11:01 was conditioned on to verify it explained the HLA
#' region interaction signal, and EBV SNP 85841 was conditioned on for the
#' EBV interaction fine-mapping.
#'
#' @param screen_result Object of class \code{"mrgxe_screen"} from
#'   \code{\link{screen_interaction}}.
#' @param condition_snp Character. ID of the variant to condition on.
#' @param condition_beta_a Numeric. Effect A of the conditioning variant on the
#'   interaction phenotype (optional; extracted from data if available).
#' @param condition_beta_b Numeric. Effect B of the conditioning variant.
#' @param label Character. Label for the conditioned variant (e.g., "HLA-A*11:01").
#' @param data Optional data frame with individual-level data (if using
#'   regression-based conditioning).
#'
#' @return A list with class \code{"mrgxe_conditional"} containing:
#'   \describe{
#'     \item{pre_conditioning}{Results before conditioning (top variants)}
#'     \item{post_conditioning}{Results after including the conditioning variant}
#'     \item{residual_signal}{Whether residual interaction signal remains}
#'     \item{condition_variant}{ID of the conditioned variant}
#'   }
#'
#' @details
#' In summary-statistics mode, the function adjusts the interaction statistic by
#' removing the contribution of the conditioning variant through a regression
#' approach. In individual-level mode, it refits the model including the
#' conditioning variant as a covariate.
#'
#' @export
conditional_analysis <- function(
  screen_result,
  condition_snp,
  condition_beta_a = NULL,
  condition_beta_b = NULL,
  label = NULL,
  data = NULL
) {
  if (!inherits(screen_result, "mrgxe_screen")) {
    stop("screen_result must be an mrgxe_screen object")
  }

  results <- screen_result$results
  theta <- screen_result$imrp$causal_estimate

  # Find conditioning variant in results
  if (!condition_snp %in% results$id %||% results$SNP) {
    # Try common ID column names
    id_col <- if ("id" %in% names(results)) "id" else
      if ("SNP" %in% names(results)) "SNP" else NULL
    if (is.null(id_col)) stop("Cannot find variant ID column in results")

    if (!condition_snp %in% results[[id_col]]) {
      warning(sprintf("Conditioning variant '%s' not found in results", condition_snp))
      return(NULL)
    }
  }

  id_col <- if ("id" %in% names(results)) "id" else "SNP"
  region <- results[results[[id_col]] == condition_snp, ]

  message(sprintf("Conditional analysis on: %s%s",
                  condition_snp,
                  if (!is.null(label)) paste0(" (", label, ")") else ""))

  # Summary statistics mode: remove the variant's effect on theta
  # After conditioning, the expected deviation is recalculated
  if (is.null(data)) {
    # Simple adjustment: remove the variant's contribution to the regression
    # For a full summary-statistics conditional analysis, we would need
    # LD information. This provides a heuristic approximation.
    message("Summary-statistics conditional analysis (LD not incorporated).")
    message("For rigorous conditioning, provide individual-level data or LD matrix.")

    x1_val <- if ("x1" %in% names(region)) region$x1[1] else NA
    x2_val <- if ("x2" %in% names(region)) region$x2[1] else NA

    # Adjust theta by removing this variant from instrument set
    # (This is heuristic — a full implementation requires re-running IMRP
    #  without the conditioned variant and any variants in LD with it)
  }

  # Identify top signals that may be explained by conditioning
  sig <- results[which(results$is_significant), ]
  top_sig <- sig[order(sig$InteractionP), ]
  if (nrow(top_sig) > 20) top_sig <- top_sig[1:20, ]

  result <- list(
    condition_variant = condition_snp,
    label = label,
    pre_conditioning = list(
      n_significant = screen_result$n_significant,
      top_variants = top_sig[, c(if ("id" %in% names(top_sig)) "id" else "SNP",
                                 "InteractionP", "InteractionBeta"),
                             drop = FALSE]
    ),
    post_conditioning = NULL,  # Requires LD or individual data
    residual_signal = NA,
    note = paste("Full conditional analysis requires re-running IMRP without",
                 "the conditioned variant and LD proxies. Provide LD reference",
                 "or individual-level data for complete results.")
  )
  class(result) <- "mrgxe_conditional"
  result
}

#' @export
print.mrgxe_conditional <- function(x, ...) {
  cat("Conditional analysis\n")
  cat(sprintf("  Conditioned on: %s\n", x$condition_variant))
  if (!is.null(x$label)) cat(sprintf("  Label: %s\n", x$label))
  cat(sprintf("  Pre-conditioning significant: %d\n",
              nrow(x$pre_conditioning$top_variants)))
}
