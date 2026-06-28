#' Select Independent Instruments
#'
#' @inheritParams MRGxE::select_instruments
#' @export
select_instruments <- function(
  data,
  p_col = NULL,
  p_threshold = 5e-8,
  ld_method = c("distance", "plink", "none"),
  window_kb = 500,
  r2_threshold = 0.1,
  plink_path = "plink",
  ld_ref_path = NULL,
  min_instruments = 10,
  relax_threshold = 5e-6,
  quiet = FALSE
) {
  ld_method <- match.arg(ld_method)
  log_msg <- function(msg) { if (!quiet) message(msg) }

  # Auto-detect P-value column
  if (is.null(p_col)) {
    for (candidate in c("P", "p_a", "p_b", "P_GWAS", "P_GWIS", "pvalue")) {
      if (candidate %in% names(data)) { p_col <- candidate; break }
    }
  }
  if (is.null(p_col)) stop("Could not auto-detect P-value column; provide p_col")

  if (!p_col %in% names(data)) {
    stop("Column '", p_col, "' not found in data. Available: ",
         paste(names(data), collapse = ", "))
  }

  candidates <- data[which(data[[p_col]] < p_threshold), ]
  used_threshold <- p_threshold
  log_msg(sprintf("Variants with %s < %.0e: %d", p_col, p_threshold, nrow(candidates)))

  if (nrow(candidates) < min_instruments && !is.null(relax_threshold)) {
    candidates <- data[which(data[[p_col]] < relax_threshold), ]
    used_threshold <- relax_threshold
    log_msg(sprintf("Relaxing threshold to %.0e (found %d candidates)",
                    relax_threshold, nrow(candidates)))
  }
  if (nrow(candidates) < 2) {
    stop("Too few significant variants (", nrow(candidates), ")")
  }

  instruments <- switch(ld_method,
    distance = .clump_by_distance(candidates, p_col = p_col, window_kb = window_kb),
    plink = .clump_by_plink(candidates, data, p_col, window_kb, r2_threshold,
                            plink_path, ld_ref_path),
    none = candidates
  )

  instruments <- instruments[order(instruments$CHR, instruments$BP), ]
  n_iv <- nrow(instruments)
  if (n_iv < 2) stop("Only ", n_iv, " independent instrument(s) found.")

  log_msg(sprintf("Selected %d independent instruments", n_iv))

  result <- list(
    instruments = instruments,
    n_instruments = n_iv,
    p_threshold_used = used_threshold,
    all_candidates = candidates,
    params = list(p_threshold = p_threshold, ld_method = ld_method,
                  window_kb = window_kb, min_instruments = min_instruments)
  )
  class(result) <- "mrgxe_instruments"
  result
}

.clump_by_distance <- function(candidates, p_col, window_kb = 500) {
  if (nrow(candidates) == 0) return(candidates[FALSE, ])
  candidates <- candidates[order(candidates$CHR, candidates[[p_col]], candidates$BP), ]
  selected <- candidates[FALSE, ]
  window_bp <- window_kb * 1000

  for (chr in sort(unique(candidates$CHR))) {
    chr_dat <- candidates[candidates$CHR == chr, ]
    chr_dat <- chr_dat[order(chr_dat[[p_col]], chr_dat$BP), ]
    kept_bp <- numeric(0)
    for (i in seq_len(nrow(chr_dat))) {
      far <- length(kept_bp) == 0 || all(abs(chr_dat$BP[i] - kept_bp) > window_bp)
      if (far) { selected <- rbind(selected, chr_dat[i, ]); kept_bp <- c(kept_bp, chr_dat$BP[i]) }
    }
  }
  selected[order(selected$CHR, selected$BP), ]
}

.clump_by_plink <- function(candidates, data, p_col, window_kb, r2, plink_path, ld_ref) {
  if (is.null(ld_ref)) stop("ld_ref_path required for PLINK clumping")
  tmp <- tempfile("clump")
  utils::write.table(data.frame(SNP = data$SNP, P = data[[p_col]]),
            file = paste0(tmp, ".input"), row.names = FALSE, quote = FALSE, sep = "\t")
  cmd <- sprintf('"%s" --bfile "%s" --clump "%s.input" --clump-p1 %.0e --clump-p2 1 --clump-r2 %.2f --clump-kb %d --out "%s" 2>&1',
                 plink_path, ld_ref, tmp, min(candidates[[p_col]], 1e-10), r2, window_kb, tmp)
  system(cmd, ignore.stdout = TRUE)
  clumped_file <- paste0(tmp, ".clumped")
  if (!file.exists(clumped_file)) return(.clump_by_distance(candidates, p_col, window_kb))
  clumped <- utils::read.table(clumped_file, header = TRUE, stringsAsFactors = FALSE)
  selected <- candidates[candidates$SNP %in% clumped$SNP, ]
  unlink(paste0(tmp, "*"))
  selected
}

#' Estimate Correlation rho
#'
#' @inheritParams MRGxE::estimate_rho
#' @export
estimate_rho <- function(
  data,
  p_a_col = NULL, p_b_col = NULL,
  null_p_threshold = 0.05,
  use_standardized = TRUE
) {
  # Auto-detect P columns
  if (is.null(p_a_col)) {
    for (c in c("P", "p_a", "P_GWAS", "pvalue")) {
      if (c %in% names(data)) { p_a_col <- c; break }
    }
  }
  if (is.null(p_b_col)) {
    for (c in c("P", "p_b", "P_GWIS", "pvalue")) {
      if (c %in% names(data)) { p_b_col <- c; break }
    }
  }
  if (is.null(p_a_col) || is.null(p_b_col)) stop("Cannot find P-value columns")

  null_idx <- data[[p_a_col]] > null_p_threshold &
    data[[p_b_col]] > null_p_threshold
  null_idx[is.na(null_idx)] <- FALSE
  null_set <- data[null_idx, ]
  if (nrow(null_set) < 10) warning("Few null variants for rho estimation")

  if (use_standardized) {
    x <- null_set$beta_b %||% null_set$BETA_GWIS %||% null_set$x1
    y <- null_set$beta_a %||% null_set$BETA_GWAS_ALIGNED %||% null_set$x2
  } else {
    x <- null_set$beta_b %||% null_set$BETA_GWIS
    y <- null_set$beta_a %||% null_set$BETA_GWAS_ALIGNED
  }

  rho <- stats::cor(x, y, use = "complete.obs")
  if (!is.finite(rho)) { warning("rho not finite; setting to 0"); rho <- 0 }
  message(sprintf("rho = %.4f (n_null = %d)", rho, nrow(null_set)))

  result <- list(rho = rho, n_null = nrow(null_set), null_p_threshold = null_p_threshold)
  class(result) <- "mrgxe_rho"
  result
}

#' @export
print.mrgxe_instruments <- function(x, ...) {
  cat(sprintf("Instruments: %d (P < %.0e)\n", x$n_instruments, x$p_threshold_used))
}
#' @export
print.mrgxe_rho <- function(x, ...) {
  cat(sprintf("rho = %.4f (N = %d)\n", x$rho, x$n_null))
}
