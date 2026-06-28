#' Generate Interaction Screening Report
#'
#' @param screen_result \code{mrgxe_screen} object.
#' @param harm_result Optional harmonization result.
#' @param iv_result Optional instrument selection result.
#' @param rho_result Optional rho estimate.
#' @param perm_result Optional permutation result.
#' @param output_file Character. Output file path. If NULL, returns lines.
#' @param title Character. Report title.
#'
#' @export
generate_report <- function(
  screen_result,
  harm_result = NULL,
  iv_result = NULL,
  rho_result = NULL,
  perm_result = NULL,
  output_file = NULL,
  title = "Interaction Screening Report"
) {
  design <- screen_result$design %||% list(design_type = "unknown")
  imrp <- screen_result$imrp

  lines <- c(
    sprintf("# %s", title),
    "",
    sprintf("**Design**: %s", design$design_type),
    sprintf("**Generated**: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "---",
    "",
    "## 1. Summary",
    "",
    sprintf("- **Variants screened**: %d", nrow(screen_result$results)),
    sprintf("- **Significant** (P < %.0e): %d",
            screen_result$threshold, screen_result$n_significant),
    sprintf("- **Theta**: %.4f (SE = %.4f, P = %.2e)",
            imrp$causal_estimate, imrp$causal_se, imrp$causal_p),
    sprintf("- **Global P (pre/aft)**: %.2e / %.2e",
            imrp$global_p_pre, imrp$global_p_aft),
    ""
  )

  # Design-specific reporting
  lines <- c(lines,
             sprintf("- **Effect A**: %s", design$label_a %||% "Effect A"),
             sprintf("- **Effect B**: %s", design$label_b %||% "Effect B"),
             sprintf("- **Interaction meaning**: %s",
                     design$interaction_meaning %||%
                     "Deviation from regression line"),
             "")

  # QC
  if (!is.null(harm_result) && !is.null(harm_result$qc_counts)) {
    lines <- c(lines, "## 2. Quality Control", "", "| Step | Count |", "|------|-------|")
    for (i in seq_len(nrow(harm_result$qc_counts))) {
      lines <- c(lines, sprintf("| %s | %d |",
                 harm_result$qc_counts$step[i], harm_result$qc_counts$n[i]))
    }
    lines <- c(lines, "")
  }

  # Instruments
  if (!is.null(iv_result)) {
    lines <- c(lines,
               sprintf("## 3. Instruments"),
               "",
               sprintf("- Count: %d", iv_result$n_instruments),
               sprintf("- P threshold: %.0e", iv_result$p_threshold_used),
               "")
  }

  # Rho
  if (!is.null(rho_result) && !is.null(rho_result$rho)) {
    lines <- c(lines,
               sprintf("## 4. Rho"), "",
               sprintf("- rho = %.4f (N null = %d)", rho_result$rho, rho_result$n_null),
               "")
  }

  # Permutation
  if (!is.null(perm_result)) {
    lines <- c(lines,
               "## 5. Permutation Thresholds", "",
               sprintf("- Permutations: %d", perm_result$n_perm),
               sprintf("- Empirical threshold (0.05): %.2e", perm_result$empirical_threshold_0.05),
               sprintf("- Empirical threshold (0.01): %.2e", perm_result$empirical_threshold_0.01),
               "")
  }

  # Significant variants
  if (screen_result$n_significant > 0) {
    sig <- screen_result$significant_variants
    lines <- c(lines, "## 6. Significant Variants", "",
               "| ID | P-value | Beta | SE |", "|-----|---------|------|-----|")

    id_col <- if ("id" %in% names(sig)) "id" else if ("SNP" %in% names(sig)) "SNP" else NULL
    sig_show <- sig[order(sig$InteractionP), ]
    if (nrow(sig_show) > 50) {
      lines <- c(lines, sprintf("*Showing top 50 of %d*", nrow(sig_show)))
      sig_show <- sig_show[1:50, ]
    }

    for (i in seq_len(nrow(sig_show))) {
      id_val <- if (!is.null(id_col)) sig_show[[id_col]][i] else i
      lines <- c(lines, sprintf("| %s | %.2e | %.4e | %.4e |",
                 id_val, sig_show$InteractionP[i],
                 sig_show$InteractionBeta[i], sig_show$InteractionSE[i]))
    }
    lines <- c(lines, "")
  }

  if (!is.null(output_file)) {
    writeLines(lines, output_file)
    message(sprintf("Report written to %s", output_file))
    invisible(lines)
  } else {
    lines
  }
}
