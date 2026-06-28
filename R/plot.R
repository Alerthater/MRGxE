#' Interaction Manhattan Plot
#'
#' @param results Data frame with P-value and chromosome/position columns.
#' @param p_col Character. P-value column name. Default "InteractionP".
#' @param chr_col, bp_col Character. Chromosome/position columns.
#' @param threshold Numeric. Significance threshold.
#' @param title Character. Plot title.
#' @param ... Additional arguments.
#'
#' @export
plot_manhattan <- function(
  results,
  p_col = "InteractionP",
  chr_col = "CHR",
  bp_col = "BP",
  threshold = 5e-8,
  title = "Interaction Manhattan plot",
  ...
) {
  # Find required columns
  if (!p_col %in% names(results)) {
    # Try alternatives
    for (alt in c("PleioP_MR", "p", "P")) {
      if (alt %in% names(results)) { p_col <- alt; break }
    }
  }
  if (!chr_col %in% names(results)) chr_col <- "CHR"
  if (!bp_col %in% names(results)) bp_col <- "BP"

  # Extract valid data
  plot_dat <- results[!is.na(results[[p_col]]) & is.finite(results[[p_col]]) &
                      results[[p_col]] > 0, ]
  if (nrow(plot_dat) == 0) stop("No valid data to plot")

  plot_dat <- plot_dat[order(plot_dat[[chr_col]], plot_dat[[bp_col]]), ]

  # Chromosome offsets
  chr_sizes <- stats::setNames(
    tapply(plot_dat[[bp_col]], plot_dat[[chr_col]], max, na.rm = TRUE),
    sort(unique(plot_dat[[chr_col]]))
  )
  chr_offsets <- cumsum(c(0, as.numeric(chr_sizes[-length(chr_sizes)])))
  names(chr_offsets) <- names(chr_sizes)

  plot_dat$global_bp <- plot_dat[[bp_col]] + chr_offsets[as.character(plot_dat[[chr_col]])]

  chr_mid <- stats::setNames(
    tapply(plot_dat$global_bp, plot_dat[[chr_col]],
           function(x) mean(range(x))),
    sort(unique(plot_dat[[chr_col]]))
  )

  colors <- ifelse(plot_dat[[chr_col]] %% 2 == 0, "#3b82f6", "#111827")
  if (threshold > 0) {
    sig_idx <- which(plot_dat[[p_col]] < threshold)
    colors[sig_idx] <- "#d73027"
  }

  ymax <- max(-log10(plot_dat[[p_col]][plot_dat[[p_col]] > 0]), na.rm = TRUE) * 1.1

  graphics::plot(
    plot_dat$global_bp, -log10(plot_dat[[p_col]]),
    pch = 16, cex = 0.45, col = colors,
    xaxt = "n",
    xlab = "Chromosome",
    ylab = expression(-log[10](P)),
    main = title,
    ylim = c(0, ymax),
    ...
  )
  graphics::axis(1, at = chr_mid, labels = names(chr_mid), cex.axis = 0.75)

  if (!is.null(threshold) && threshold > 0) {
    graphics::abline(h = -log10(threshold), col = "#d73027", lty = 2, lwd = 1.5)
  }

  invisible(plot_dat)
}

#' @rdname plot_manhattan
#' @export
plot_tmrgxe_manhattan <- function(...) {
  plot_manhattan(...)
}

#' Interaction QQ Plot
#'
#' @param results Data frame with P-value column.
#' @param p_col Character. P-value column name.
#' @param title Character. Plot title.
#' @param ... Additional arguments.
#'
#' @export
plot_qq <- function(
  results,
  p_col = "InteractionP",
  title = "Interaction QQ plot",
  ...
) {
  if (!p_col %in% names(results)) {
    for (alt in c("PleioP_MR", "p", "P")) {
      if (alt %in% names(results)) { p_col <- alt; break }
    }
  }

  p <- sort(results[[p_col]][results[[p_col]] > 0], na.last = NA)
  if (length(p) < 10) stop("Too few valid P-values for QQ plot")

  expected <- -log10(stats::ppoints(length(p)))
  observed <- -log10(p)

  chisq <- stats::qchisq(1 - p, 1)
  lambda_gc <- median(chisq, na.rm = TRUE) / stats::qchisq(0.5, 1)

  graphics::plot(
    expected, observed,
    pch = 16, cex = 0.55, col = "#37415199",
    xlab = expression(Expected ~ -log[10](P)),
    ylab = expression(Observed ~ -log[10](P)),
    main = if (is.finite(lambda_gc))
      paste0(title, "  (\U03BB = ", sprintf("%.3f", lambda_gc), ")")
      else title,
    ...
  )
  graphics::abline(0, 1, col = "#d73027", lwd = 2)

  invisible(data.frame(expected = expected, observed = observed, lambda_gc = lambda_gc))
}

#' @rdname plot_qq
#' @export
plot_tmrgxe_qq <- function(...) {
  plot_qq(...)
}

#' Interaction Scatter Plot
#'
#' @param results Data frame with effect columns.
#' @param imrp_result Optional IMRP result (for theta line).
#' @param beta_a_col, beta_b_col Character. Effect column names.
#' @param title Character. Plot title.
#' @param ... Additional arguments.
#'
#' @export
plot_scatter <- function(
  results,
  imrp_result = NULL,
  beta_a_col = NULL,
  beta_b_col = NULL,
  title = "Effect B vs Effect A",
  ...
) {
  # Detect columns
  if (is.null(beta_b_col)) {
    for (c in c("x1", "beta_b", "BETA_GWIS")) {
      if (c %in% names(results)) { beta_b_col <- c; break }
    }
  }
  if (is.null(beta_a_col)) {
    for (c in c("x2", "beta_a", "BETA_GWAS_ALIGNED", "BETA_GWAS")) {
      if (c %in% names(results)) { beta_a_col <- c; break }
    }
  }

  if (is.null(beta_b_col) || is.null(beta_a_col)) {
    stop("Cannot detect effect columns in results")
  }

  x <- results[[beta_b_col]]
  y <- results[[beta_a_col]]

  slope <- NULL
  if (!is.null(imrp_result)) {
    if (is.list(imrp_result) && !is.null(imrp_result$causal_estimate)) {
      slope <- imrp_result$causal_estimate
    }
  }
  if (is.null(slope) || !is.finite(slope)) {
    slope <- stats::lm(y ~ x)$coefficients[2]
  }

  # Colors for significant variants
  sig_col <- if ("is_significant" %in% names(results)) {
    ifelse(results$is_significant, "#d73027", "#4d4d4d55")
  } else if ("is_significant_gxe" %in% names(results)) {
    ifelse(results$is_significant_gxe, "#d73027", "#4d4d4d55")
  } else {
    "#4d4d4d55"
  }

  graphics::plot(
    x, y,
    pch = 16, col = sig_col,
    xlab = "Effect B (exposure)",
    ylab = "Effect A (outcome)",
    main = title,
    ...
  )
  graphics::abline(a = 0, b = slope, col = "#2b6cb0", lwd = 2)
  graphics::legend("topleft",
    legend = c("All variants", "Significant", "Regression"),
    col = c("#4d4d4d55", "#d73027", "#2b6cb0"),
    pch = c(16, 16, NA), lty = c(NA, NA, 1), bty = "n")
}

#' @rdname plot_scatter
#' @export
plot_tmrgxe_scatter <- function(...) {
  plot_scatter(...)
}

#' Conditional Analysis Forest Plot
#'
#' Shows how interaction estimates change before and after conditioning on a
#' key variant.
#'
#' @param conditional_result Object from \code{\link{conditional_analysis}}.
#'
#' @export
plot_conditional <- function(conditional_result) {
  if (!inherits(conditional_result, "mrgxe_conditional")) {
    stop("conditional_result must be from conditional_analysis()")
  }

  pre <- conditional_result$pre_conditioning$top_variants

  if (nrow(pre) == 0) {
    message("No significant variants to plot")
    return(invisible(NULL))
  }

  id_col <- if ("id" %in% names(pre)) "id" else "SNP"

  graphics::plot(
    pre$InteractionBeta, -log10(pre$InteractionP),
    pch = 16, col = "#2563eb",
    xlab = "Interaction Beta",
    ylab = expression(-log[10](P)),
    main = sprintf("Top variants (conditioned on: %s)",
                   conditional_result$condition_variant)
  )
  graphics::text(pre$InteractionBeta, -log10(pre$InteractionP),
        labels = pre[[id_col]], pos = 4, cex = 0.6)
}
