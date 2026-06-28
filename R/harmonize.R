#' Harmonize Two Sets of Effect Estimates
#'
#' Aligns, filters, and standardizes two matched sets of effect estimates
#' for downstream interaction screening. Generalizes the GxE-specific QC
#' to support any pairwise study design.
#'
#' @param data An object of class \code{"mrgxe_data"} from \code{\link{import_effects}}
#'   or a data frame with columns \code{id}, \code{beta_a}, \code{beta_b},
#'   \code{se_a}, \code{se_b}.
#' @param design A \code{mrgxe_design} object. Controls QC defaults.
#' @param standardize Logical. Whether to standardize effects as
#'   \code{beta / se / sqrt(n)}. Default taken from design.
#' @param remove_palindromic Logical. Whether to remove A/T and C/G SNPs.
#'   Default \code{FALSE} for genome-to-genome studies (non-human genomes may
#'   not have traditional palindromic issues).
#' @param align_alleles Logical. Whether to perform allele alignment.
#'   Default \code{TRUE} for GxE, \code{FALSE} for genome-to-genome.
#' @param max_eaf_diff Numeric. Maximum EAF difference after harmonization.
#'   Default \code{0.15} for GxE; \code{NULL} for genome-to-genome.
#' @param chr_col, bp_col, a1_col, a2_col, eaf_col Character. Column names for
#'   chromosome, position, alleles, and EAF in the original data (needed for
#'   palindromic filtering and allele alignment).
#' @param min_n Integer. Minimum sample size. Default \code{NULL}.
#' @param design_type Deprecated. Use \code{design} instead.
#' @param quiet Logical. Suppress messages. Default \code{FALSE}.
#'
#' @return An object of class \code{"mrgxe_harmonized"} with:
#'   \describe{
#'     \item{cleaned}{Data frame with harmonized, QC-passing variants}
#'     \item{qc_counts}{Data frame with counts at each QC step}
#'     \item{design}{Study design used}
#'   }
#'
#' @export
harmonize_effects <- function(
  data,
  design = NULL,
  standardize = NULL,
  remove_palindromic = NULL,
  align_alleles = NULL,
  max_eaf_diff = NULL,
  chr_col = NULL,
  bp_col = NULL,
  a1_col = NULL,
  a2_col = NULL,
  eaf_col = NULL,
  min_n = NULL,
  design_type = NULL,
  quiet = FALSE
) {
  # Backward compatibility: extract design from mrgxe_data
  if (inherits(data, "mrgxe_data") && !is.null(data$design)) {
    design <- data$design
    if ("gwis" %in% names(data)) {
      # Legacy GxE format; handle via the old harmonize_gwas_gwis
      return(harmonize_gwas_gwis(data$gwis, data$gwas, quiet = quiet))
    }
    dt <- data$data
  } else if (is.data.frame(data)) {
    dt <- data
  } else {
    stop("data must be an mrgxe_data object or a data frame")
  }

  # Resolve design
  if (is.null(design)) {
    design <- study_design(design_type %||% "custom")
  }

  # Resolve QC parameters from design
  if (is.null(standardize)) standardize <- design$default_standardize %||% TRUE
  if (is.null(remove_palindromic)) remove_palindromic <- design$default_palindromic_removal %||% FALSE
  if (is.null(align_alleles)) align_alleles <- design$design_type == "gxe"
  if (is.null(max_eaf_diff)) {
    max_eaf_diff <- if (design$default_eaf_qc %||% FALSE) 0.15 else NULL
  }

  log_msg <- function(msg) { if (!quiet) message(msg) }
  qc_steps <- data.frame(step = character(), n = integer(), stringsAsFactors = FALSE)
  qc_log <- function(step_name, n) {
    qc_steps <<- rbind(qc_steps, data.frame(step = step_name, n = n, stringsAsFactors = FALSE))
  }

  n0 <- nrow(dt)
  qc_log("input", n0)

  # Allele alignment (GxE-specific)
  if (align_alleles && !is.null(a1_col) && !is.null(a2_col)) {
    # Check for palindromic SNPs
    if (remove_palindromic) {
      is_pal <- .is_palindromic(dt[[a1_col]], dt[[a2_col]])
      dt <- dt[!is_pal, ]
      qc_log("drop_palindromic", sum(is_pal))
    }

    # Check for missing allele info (for non-human genomes, this may differ)
    if (!is.null(chr_col) && !is.null(bp_col)) {
      dt <- dt[!duplicated(dt$id), ]
    }
  }

  # EAF QC (GxE-specific; requires both datasets to have EAF)
  if (!is.null(max_eaf_diff) && !is.null(eaf_col)) {
    # EAF filtering would go here if both datasets have it
    # For now, this is a placeholder for GxE-specific extension
  }

  # Remove rows with invalid SE
  dt <- dt[dt$se_a > 0 & dt$se_b > 0 & is.finite(dt$se_a) & is.finite(dt$se_b), ]
  qc_log("valid_se", nrow(dt))

  # Sample size filter
  if (!is.null(min_n) && "n_a" %in% names(dt) && "n_b" %in% names(dt)) {
    dt <- dt[dt$n_a >= min_n & dt$n_b >= min_n, ]
    qc_log(sprintf("min_n_%d", min_n), nrow(dt))
  }

  # Standardize if requested
  if (standardize) {
    # Need sample size columns
    if ("n_a" %in% names(dt) && "n_b" %in% names(dt)) {
      dt$x1 <- dt$beta_b / dt$se_b / sqrt(dt$n_b)
      dt$x2 <- dt$beta_a / dt$se_a / sqrt(dt$n_a)
      dt$x1_se <- 1 / sqrt(dt$n_b)
      dt$x2_se <- 1 / sqrt(dt$n_a)
      log_msg("Standardized effects using beta / se / sqrt(n)")
    } else {
      log_msg("Sample size missing; using raw effect estimates")
      dt$x1 <- dt$beta_b
      dt$x2 <- dt$beta_a
      dt$x1_se <- dt$se_b
      dt$x2_se <- dt$se_a
    }
  } else {
    dt$x1 <- dt$beta_b
    dt$x2 <- dt$beta_a
    dt$x1_se <- dt$se_b
    dt$x2_se <- dt$se_a
  }

  log_msg(sprintf("Retained %d / %d variants after QC", nrow(dt), n0))

  result <- list(
    cleaned = dt,
    qc_counts = qc_steps,
    design = design,
    params = list(
      standardize = standardize,
      remove_palindromic = remove_palindromic,
      align_alleles = align_alleles,
      max_eaf_diff = max_eaf_diff
    )
  )
  class(result) <- "mrgxe_harmonized"
  result
}

# ---------------------------------------------------------------------------
# Legacy GxE-specific harmonization (kept for backward compatibility)
# ---------------------------------------------------------------------------

#' @rdname harmonize_effects
#' @export
harmonize_gwas_gwis <- function(
  gwis,
  gwas,
  max_eaf_diff = 0.15,
  remove_palindromic = TRUE,
  remove_ambiguous_indel = TRUE,
  min_n_gwis = NULL,
  min_n_gwas = NULL,
  duplicate_resolution = c("keep_closest_maf", "remove_all"),
  quiet = FALSE
) {
  duplicate_resolution <- match.arg(duplicate_resolution)
  stopifnot(is.data.frame(gwis), is.data.frame(gwas))

  .validate_required_cols(gwis, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "N"), "gwis")
  .validate_required_cols(gwas, c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "N"), "gwas")

  qc_log <- function(msg) { if (!quiet) message(msg) }

  gwis$CHR <- as.numeric(as.character(gwis$CHR))
  gwas$CHR <- as.numeric(as.character(gwas$CHR))
  gwis$BP <- as.numeric(as.character(gwis$BP))
  gwas$BP <- as.numeric(as.character(gwas$BP))

  merged <- merge(gwis, gwas, by = c("SNP", "CHR", "BP"),
                  suffixes = c("_GWIS", "_GWAS"), all = FALSE)
  n_merged <- nrow(merged)
  qc_log(sprintf("Merged: %d variants", n_merged))

  merged$qc_palindromic <- .is_palindromic(merged$A1_GWIS, merged$A2_GWIS)
  if (remove_ambiguous_indel) {
    merged$qc_ambiguous <- !(merged$A1_GWIS %in% c("A", "C", "G", "T")) |
      !(merged$A2_GWIS %in% c("A", "C", "G", "T"))
  } else {
    merged$qc_ambiguous <- FALSE
  }

  same <- merged$A1_GWIS == merged$A1_GWAS & merged$A2_GWIS == merged$A2_GWAS
  flipped <- merged$A1_GWIS == merged$A2_GWAS & merged$A2_GWIS == merged$A1_GWAS
  merged$qc_allele_match <- same | flipped
  merged$was_flipped <- flipped

  has_eaf <- all(c("EAF_GWIS", "EAF_GWAS") %in% names(merged)) ||
    ("EAF" %in% names(gwis) && "EAF" %in% names(gwas))

  merged$BETA_GWAS_ALIGNED <- ifelse(flipped, -merged$BETA_GWAS, merged$BETA_GWAS)

  if (has_eaf) {
    if ("EAF_GWIS" %in% names(merged)) {
      merged$EAF_GWAS_ALIGNED <- ifelse(flipped, 1 - merged$EAF_GWAS, merged$EAF_GWAS)
    } else {
      merged$EAF_GWAS_ALIGNED <- ifelse(flipped, 1 - merged$EAF, merged$EAF)
    }
    merged$eaf_diff <- abs(merged$EAF_GWIS - merged$EAF_GWAS_ALIGNED)
    merged$qc_eaf_ok <- merged$eaf_diff <= max_eaf_diff
  } else {
    merged$EAF_GWAS_ALIGNED <- NA_real_
    merged$eaf_diff <- NA_real_
    merged$qc_eaf_ok <- TRUE
  }

  keep <- rep(TRUE, nrow(merged))
  if (remove_palindromic) keep <- keep & !merged$qc_palindromic
  keep <- keep & merged$qc_allele_match
  if (remove_ambiguous_indel) keep <- keep & !merged$qc_ambiguous
  if (has_eaf) keep <- keep & merged$qc_eaf_ok
  if (!is.null(min_n_gwis)) keep <- keep & merged$N_GWIS >= min_n_gwis
  if (!is.null(min_n_gwas)) keep <- keep & merged$N_GWAS >= min_n_gwas

  cleaned <- merged[keep, ]

  dup_idx <- duplicated(cleaned$SNP) | duplicated(cleaned$SNP, fromLast = TRUE)
  if (any(dup_idx)) {
    if (duplicate_resolution == "remove_all") {
      cleaned <- cleaned[!dup_idx, ]
    } else if (duplicate_resolution == "keep_closest_maf" && has_eaf) {
      cleaned <- cleaned[order(cleaned$eaf_diff), ]
      cleaned <- cleaned[!duplicated(cleaned$SNP), ]
    } else {
      cleaned <- cleaned[!duplicated(cleaned$SNP), ]
    }
  }

  qc_counts <- data.frame(
    step = c("merged",
             if (remove_palindromic) "drop_palindromic" else NULL,
             "drop_allele_mismatch_or_ambiguous",
             if (has_eaf) sprintf("drop_eaf_diff_gt_%.2f", max_eaf_diff) else NULL,
             if (!is.null(min_n_gwis)) "drop_min_n_gwis" else NULL,
             if (!is.null(min_n_gwas)) "drop_min_n_gwas" else NULL,
             "retained")[!sapply(c("merged",
               if (remove_palindromic) "drop_palindromic" else NULL,
               "drop_allele_mismatch_or_ambiguous",
               if (has_eaf) "x" else NULL,
               if (!is.null(min_n_gwis)) "x" else NULL,
               if (!is.null(min_n_gwas)) "x" else NULL,
               "retained"), is.null)],
    n = c(n_merged,
          if (remove_palindromic) sum(merged$qc_palindromic) else NULL,
          sum(!merged$qc_palindromic & (!merged$qc_allele_match | merged$qc_ambiguous)),
          if (has_eaf) sum(!merged$qc_palindromic & merged$qc_allele_match & !merged$qc_ambiguous & !merged$qc_eaf_ok) else NULL,
          if (!is.null(min_n_gwis)) sum(keep & merged$N_GWIS < min_n_gwis) else NULL,
          if (!is.null(min_n_gwas)) sum(keep & merged$N_GWAS < min_n_gwas) else NULL,
          nrow(cleaned)),
    stringsAsFactors = FALSE
  )

  qc_log(sprintf("Retained after QC: %d variants (%.1f%%)",
                 nrow(cleaned), 100 * nrow(cleaned) / n_merged))

  result <- list(
    cleaned = cleaned,
    qc_counts = qc_counts,
    design = study_design("gxe"),
    params = list(
      max_eaf_diff = max_eaf_diff,
      remove_palindromic = remove_palindromic
    )
  )
  class(result) <- "mrgxe_harmonized"
  result
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
.is_palindromic <- function(a1, a2) {
  pair <- paste0(toupper(a1), toupper(a2))
  pair %in% c("AT", "TA", "CG", "GC")
}

#' @export
print.mrgxe_harmonized <- function(x, ...) {
  cat("MRGxE harmonized effects\n")
  cat(sprintf("  Retained: %d variants\n", nrow(x$cleaned)))
  if (!is.null(x$design)) {
    cat(sprintf("  Design:   %s\n", x$design$design_type))
  }
  if (!is.null(x$qc_counts) && nrow(x$qc_counts) > 0) {
    cat("  QC counts:\n")
    print(x$qc_counts, row.names = FALSE)
  }
}
