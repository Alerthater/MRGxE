#' Import and Harmonize Effect Estimates
#'
#' Reads two sets of genetic effect estimates (from files or data frames),
#' aligns them by variant identifier, and returns a standardized data matrix
#' ready for downstream analysis. This replaces the GxE-specific
#' \code{import_summary_stats} with a more general interface.
#'
#' @param data_a Data frame or character. Effect estimates "A" (y-axis).
#'   If character, treated as file path.
#' @param data_b Data frame or character. Effect estimates "B" (x-axis).
#'   If character, treated as file path.
#' @param design A \code{mrgxe_design} object from \code{\link{study_design}}.
#'   Determines default parameters and QC rules.
#' @param id_col Character. Column name for variant identifier in both datasets.
#'   Default \code{"SNP"}.
#' @param beta_a_col Character. Column name for effect A estimate.
#' @param beta_b_col Character. Column name for effect B estimate.
#' @param se_a_col Character. Column name for standard error of effect A.
#' @param se_b_col Character. Column name for standard error of effect B.
#' @param n_a_col Character. Column name for sample size of effect A (optional).
#' @param n_b_col Character. Column name for sample size of effect B (optional).
#' @param p_a_col Character. Column name for P-value of effect A (optional).
#' @param p_b_col Character. Column name for P-value of effect B (optional).
#' @param sep Character. Field separator when reading files. Default \code{"\t"}.
#' @param header Logical. Whether files have headers. Default \code{TRUE}.
#'
#' @return An object of class \code{"mrgxe_data"} containing:
#'   \describe{
#'     \item{data}{Merged data frame with columns \code{id}, \code{beta_a},
#'       \code{beta_b}, \code{se_a}, \code{se_b}, and optionally \code{n_a},
#'       \code{n_b}, \code{p_a}, \code{p_b}}
#'     \item{design}{The study design object}
#'     \item{n_variants}{Number of matched variants}
#'   }
#'
#' @examples
#' \dontrun{
#' gxe <- study_design("gxe")
#' d <- import_effects("gwas.txt", "gwis.txt", design = gxe,
#'                     beta_a_col = "BETA", beta_b_col = "Effect",
#'                     se_a_col = "SE", se_b_col = "StdErr")
#' }
#'
#' @export
import_effects <- function(
  data_a,
  data_b,
  design = NULL,
  id_col = "SNP",
  beta_a_col = NULL,
  beta_b_col = NULL,
  se_a_col = NULL,
  se_b_col = NULL,
  n_a_col = NULL,
  n_b_col = NULL,
  p_a_col = NULL,
  p_b_col = NULL,
  chr_col = NULL,
  bp_col = NULL,
  sep = "\t",
  header = TRUE
) {
  if (is.null(design)) {
    design <- study_design("custom")
  }

  # Read files if paths are provided
  if (is.character(data_a)) {
    data_a <- utils::read.table(data_a, sep = sep, header = header,
                                 stringsAsFactors = FALSE, comment.char = "")
  }
  if (is.character(data_b)) {
    data_b <- utils::read.table(data_b, sep = sep, header = header,
                                 stringsAsFactors = FALSE, comment.char = "")
  }

  # Ensure both are data frames
  data_a <- as.data.frame(data_a)
  data_b <- as.data.frame(data_b)

  # Validate ID column existence
  if (!id_col %in% names(data_a)) stop(sprintf("id_col '%s' not found in data_a", id_col))
  if (!id_col %in% names(data_b)) stop(sprintf("id_col '%s' not found in data_b", id_col))

  # Auto-detect column names if not provided
  # (Look for likely candidates based on the design type)
  beta_a_col <- beta_a_col %||% .detect_beta_col(data_a, "a")
  beta_b_col <- beta_b_col %||% .detect_beta_col(data_b, "b")
  se_a_col <- se_a_col %||% .detect_se_col(data_a, "a")
  se_b_col <- se_b_col %||% .detect_se_col(data_b, "b")
  chr_col <- chr_col %||% .detect_chr_col(data_a)
  bp_col <- bp_col %||% .detect_bp_col(data_a)

  if (is.null(beta_a_col)) stop("Could not auto-detect beta_a column; provide beta_a_col")
  if (is.null(beta_b_col)) stop("Could not auto-detect beta_b column; provide beta_b_col")
  if (is.null(se_a_col)) stop("Could not auto-detect se_a column; provide se_a_col")
  if (is.null(se_b_col)) stop("Could not auto-detect se_b column; provide se_b_col")

  # Build columns to select from each dataset
  cols_a <- unique(c(id_col, beta_a_col, se_a_col, n_a_col, p_a_col, chr_col, bp_col))
  cols_b <- unique(c(id_col, beta_b_col, se_b_col, n_b_col, p_b_col))

  # Merge by ID
  merged <- merge(
    data_a[, intersect(cols_a, names(data_a)), drop = FALSE],
    data_b[, intersect(cols_b, names(data_b)), drop = FALSE],
    by = id_col, suffixes = c("_A", "_B"), all = FALSE
  )

  # Helper: get merged column name after merge suffix mangling
  .col_after_merge <- function(base_col, suffix) {
    if (base_col == id_col) return(id_col)
    c1 <- paste0(base_col, suffix)   # e.g. "beta_a_A"
    c2 <- base_col                   # e.g. "beta_a" (when no name collision)
    if (c1 %in% names(merged)) return(c1)
    if (c2 %in% names(merged)) return(c2)
    stop(sprintf("Column '%s' (or '%s') not found after merge", c1, c2))
  }

  # Standardize output
  out <- data.frame(
    id = merged[[id_col]],
    beta_a = as.numeric(merged[[.col_after_merge(beta_a_col, "_A")]]),
    beta_b = as.numeric(merged[[.col_after_merge(beta_b_col, "_B")]]),
    se_a = as.numeric(merged[[.col_after_merge(se_a_col, "_A")]]),
    se_b = as.numeric(merged[[.col_after_merge(se_b_col, "_B")]]),
    stringsAsFactors = FALSE
  )

  # Add genomic coordinates from data_a side
  if (!is.null(chr_col)) {
    chr_merged <- .col_after_merge(chr_col, "_A")
    out$CHR <- as.numeric(merged[[chr_merged]])
  }
  if (!is.null(bp_col)) {
    bp_merged <- .col_after_merge(bp_col, "_A")
    out$BP <- as.numeric(merged[[bp_merged]])
  }

  # Add optional columns
  if (!is.null(n_a_col) && n_a_col %in% names(data_a)) {
    out$n_a <- as.numeric(merged[[n_a_col]])
  }
  if (!is.null(n_b_col) && n_b_col %in% names(data_b)) {
    out$n_b <- as.numeric(merged[[n_b_col]])
  }
  if (!is.null(p_a_col) && p_a_col %in% names(data_a)) {
    out$p_a <- as.numeric(merged[[p_a_col]])
  }
  if (!is.null(p_b_col) && p_b_col %in% names(data_b)) {
    out$p_b <- as.numeric(merged[[p_b_col]])
  }

  # Remove missing values
  out <- out[!is.na(out$beta_a) & !is.na(out$beta_b) &
             !is.na(out$se_a) & !is.na(out$se_b) &
             is.finite(out$se_a) & out$se_a > 0 &
             is.finite(out$se_b) & out$se_b > 0, ]

  message(sprintf("Imported %d matched variants", nrow(out)))

  result <- list(
    data = out,
    design = design,
    n_variants = nrow(out)
  )
  class(result) <- "mrgxe_data"
  result
}

# ---------------------------------------------------------------------------
# Legacy wrapper
# ---------------------------------------------------------------------------

#' @rdname import_effects
#' @export
import_summary_stats <- function(
  gwas_path,
  gwis_path,
  gwas_format = c("custom", "GLGC", "UKBB"),
  gwis_format = c("custom", "CHARGE"),
  col_map = NULL,
  sep = "\t",
  header = TRUE,
  env_variable = NULL
) {
  # Legacy wrapper for the GxE-specific import
  # Delegates to import_effects with appropriate column mapping
  gwas_format <- match.arg(gwas_format)
  gwis_format <- match.arg(gwis_format)

  if (gwas_format == "custom" && is.null(col_map)) {
    stop("col_map is required when gwas_format = 'custom'")
  }

  gwas_map <- .get_column_map(gwas_format, col_map, "gwas")
  gwis_map <- .get_column_map(gwis_format, col_map, "gwis")

  # Read files directly for backward compat
  gwas_raw <- utils::read.table(gwas_path, sep = sep, header = header,
                                stringsAsFactors = FALSE, comment.char = "")
  gwis_raw <- utils::read.table(gwis_path, sep = sep, header = header,
                                stringsAsFactors = FALSE, comment.char = "")

  gwis <- .standardize_data(gwis_raw, gwis_map, "GWIS")
  gwas <- .standardize_data(gwas_raw, gwas_map, "GWAS")

  .validate_required_cols(gwis, c("SNP", "BETA", "SE", "N"), "GWIS")
  .validate_required_cols(gwas, c("SNP", "BETA", "SE", "N"), "GWAS")

  na_idx_gwis <- is.na(gwis$BETA) | is.na(gwis$SE) | is.na(gwis$N) | gwis$N <= 0
  na_idx_gwas <- is.na(gwas$BETA) | is.na(gwas$SE) | is.na(gwas$N) | gwas$N <= 0
  gwis <- gwis[!na_idx_gwis, ]
  gwas <- gwas[!na_idx_gwas, ]

  if (nrow(gwis) == 0) stop("No valid rows remaining in GWIS data")
  if (nrow(gwas) == 0) stop("No valid rows remaining in GWAS data")

  message(sprintf("Imported GWIS: %d variants", nrow(gwis)))
  message(sprintf("Imported GWAS: %d variants", nrow(gwas)))

  result <- list(
    gwis = gwis,
    gwas = gwas,
    params = list(
      gwas_path = gwas_path,
      gwis_path = gwis_path,
      env_variable = env_variable,
      gwas_format = gwas_format,
      gwis_format = gwis_format,
      call = match.call()
    )
  )
  class(result) <- "mrgxe_data"
  result
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.detect_beta_col <- function(df, label) {
  candidates <- c("BETA", "Beta", "beta", "Effect", "EFFECT",
                  "BETA_GWAS", "BETA_GWIS", "METAL_Effect",
                  if (label == "a") c("BETA_A", "beta_a") else c("BETA_B", "beta_b"))
  for (c in candidates) {
    if (c %in% names(df)) return(c)
  }
  NULL
}

.detect_se_col <- function(df, label) {
  candidates <- c("SE", "StdErr", "se", "std_err",
                  "SE_GWAS", "SE_GWIS", "METAL_StdErr",
                  if (label == "a") c("SE_A", "se_a") else c("SE_B", "se_b"))
  for (c in candidates) {
    if (c %in% names(df)) return(c)
  }
  NULL
}

.detect_chr_col <- function(df) {
  candidates <- c("CHR", "Chr", "chr", "CHROM", "Chromosome")
  for (c in candidates) {
    if (c %in% names(df)) return(c)
  }
  NULL
}

.detect_bp_col <- function(df) {
  candidates <- c("BP", "Bp", "bp", "POS", "Pos", "position")
  for (c in candidates) {
    if (c %in% names(df)) return(c)
  }
  NULL
}

.get_column_map <- function(format, col_map, data_type) {
  if (format %in% c("GLGC", "UKBB")) {
    base <- list(
      gwas = list(
        snp = "rsID", chr = "CHR", bp = "POS",
        a1 = "Effect_allele", a2 = "Other_allele",
        beta = "Beta", se = "SE", n = "N",
        eaf = "EAF", p = "Pvalue"
      ),
      gwis = list(
        snp = "MarkerName", chr = "CHROM", bp = "POS_b37",
        a1 = "Allele1", a2 = "Allele2",
        beta = "Effect", se = "StdErr", n = "TotalSampleSize",
        eaf = "Freq1", p = "Pvalue"
      )
    )
    return(base[[data_type]])
  }
  col_map
}

.standardize_data <- function(raw, map, suffix) {
  col_map <- list(
    SNP = map$snp,
    CHR = map$chr,
    BP = map$bp,
    A1 = map$a1,
    A2 = map$a2,
    BETA = map$beta,
    SE = map$se,
    N = map$n,
    EAF = map$eaf,
    P = map$p
  )

  out <- data.frame(
    SNP = as.character(raw[[col_map$SNP]]),
    CHR = as.numeric(as.character(raw[[col_map$CHR]])),
    BP = as.numeric(as.character(raw[[col_map$BP]])),
    A1 = toupper(as.character(raw[[col_map$A1]])),
    A2 = toupper(as.character(raw[[col_map$A2]])),
    BETA = as.numeric(raw[[col_map$BETA]]),
    SE = as.numeric(raw[[col_map$SE]]),
    N = as.numeric(raw[[col_map$N]]),
    stringsAsFactors = FALSE
  )

  if (!is.null(col_map$EAF) && col_map$EAF %in% names(raw)) {
    out$EAF <- as.numeric(raw[[col_map$EAF]])
  }
  if (!is.null(col_map$P) && col_map$P %in% names(raw)) {
    out$P <- as.numeric(raw[[col_map$P]])
  }

  attr(out, "data_type") <- suffix
  out
}

.validate_required_cols <- function(df, cols, name) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns in %s data: %s",
                 name, paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

#' @export
print.mrgxe_data <- function(x, ...) {
  cat("MRGxE imported data\n")
  cat(sprintf("  Variants:  %d\n", x$n_variants %||% nrow(x$data %||% x$gwis)))
  if (!is.null(x$design)) {
    cat(sprintf("  Design:    %s\n", x$design$design_type))
    cat(sprintf("  Effect A:  %s\n", x$design$label_a))
    cat(sprintf("  Effect B:  %s\n", x$design$label_b))
  }
}
