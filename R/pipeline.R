#' Run Interaction Screening Pipeline
#'
#' High-level pipeline that handles any study design (GxE, genome-to-genome,
#' custom) from data import to reporting.
#'
#' @param effect_a_path Character. Path to effect A data (y-axis).
#' @param effect_b_path Character. Path to effect B data (x-axis).
#' @param design_type Character. Study design: \code{"gxe"}, \code{"genome_to_genome"},
#'   or \code{"custom"}. Default \code{"gxe"}.
#' @param design Optional pre-built \code{mrgxe_design} object (overrides design_type).
#' @param id_col, beta_a_col, beta_b_col, se_a_col, se_b_col, n_a_col, n_b_col,
#'   p_a_col, p_b_col Character. Column mappings (see \code{\link{import_effects}}).
#' @param config List or character. Config file path or list of parameters.
#' @param output_dir Character. Output directory.
#' @param run_plots Logical. Generate plots.
#' @param run_report Logical. Generate report.
#' @param n_perm Integer. Permutations for empirical threshold.
#'   Default 0 (skip permutation).
#' @param ... Additional parameters.
#'
#' @return \code{mrgxe_pipeline} object.
#'
#' @export
run_pipeline <- function(
  effect_a_path = NULL,
  effect_b_path = NULL,
  design_type = "gxe",
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
  config = NULL,
  output_dir = "./results",
  run_plots = TRUE,
  run_report = TRUE,
  n_perm = 0,
  ...
) {
  # Resolve configuration
  if (is.character(config)) {
    config <- read_config(config)
  }

  params <- list(...)
  if (is.list(config) && !is.null(config$data)) {
    if (is.null(effect_a_path)) effect_a_path <- config$data$effect_a_path
    if (is.null(effect_b_path)) effect_b_path <- config$data$effect_b_path
    design_type <- config$data$design_type %||% design_type
    params <- modifyList(params, config$params %||% list())
  }

  msg <- function(...) { if (!quiet) message(...) }
  quiet <- params$quiet %||% FALSE

  # Create design
  if (is.null(design)) {
    design <- study_design(design_type)
  }

  # Output dirs
  table_dir <- file.path(output_dir, "tables")
  figure_dir <- file.path(output_dir, "figures")
  for (d in c(output_dir, table_dir, figure_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  # ---- Step 1: Import ----
  msg(sprintf("=== Step 1: Import effects (%s design) ===", design$design_type))
  data <- import_effects(
    data_a = effect_a_path,
    data_b = effect_b_path,
    design = design,
    id_col = id_col,
    beta_a_col = beta_a_col,
    beta_b_col = beta_b_col,
    se_a_col = se_a_col,
    se_b_col = se_b_col,
    n_a_col = n_a_col,
    n_b_col = n_b_col,
    p_a_col = p_a_col,
    p_b_col = p_b_col
  )

  # ---- Step 2: Harmonize ----
  msg("=== Step 2: Harmonize effects ===")
  harm <- harmonize_effects(
    data,
    design = design,
    standardize = design$default_standardize,
    quiet = quiet
  )
  utils::write.csv(harm$qc_counts,
            file.path(table_dir, "qc_counts.csv"), row.names = FALSE)

  # ---- Step 3: Select instruments ----
  msg("=== Step 3: Select instruments ===")
  p_threshold <- params$iv_p_threshold %||% design$default_iv_p_threshold %||% 5e-8
  iv <- select_instruments(
    harm$cleaned,
    p_col = "P",
    p_threshold = p_threshold,
    ld_method = params$ld_method %||% "distance",
    window_kb = params$ld_window_kb %||% 500,
    min_instruments = params$min_instruments %||% 10,
    relax_threshold = params$relax_threshold %||% 5e-6,
    quiet = quiet
  )
  utils::write.csv(iv$instruments,
            file.path(table_dir, "selected_instruments.csv"), row.names = FALSE)

  # ---- Step 4: Estimate rho (if needed) ----
  if (design$requires_rho) {
    msg("=== Step 4: Estimate rho ===")
    rho <- estimate_rho(harm$cleaned,
                        null_p_threshold = params$null_p_threshold %||% 0.05)
    rho_val <- rho$rho
  } else {
    msg("=== Step 4: rho not required for this design; using 0 ===")
    rho_val <- 0
    rho <- list(rho = 0, n_null = 0)
  }

  # ---- Step 5: Screen ----
  msg("=== Step 5: Interaction screening ===")
  screen <- screen_interaction(
    data = harm$cleaned,
    instruments = iv$instruments,
    rho = rho_val,
    standardize = design$default_standardize %||% TRUE,
    method = params$mr_method %||% "IVW",
    signif_threshold = params$screen_threshold %||% 5e-8,
    design = design
  )
  utils::write.csv(screen$results,
            file.path(table_dir, "interaction_screen_all.csv"), row.names = FALSE)
  utils::write.csv(screen$significant_variants,
            file.path(table_dir, "interaction_significant.csv"), row.names = FALSE)

  # ---- Step 6: Permutation (optional) ----
  perm_result <- NULL
  if (n_perm > 0) {
    msg(sprintf("=== Step 6: Permutation testing (%d runs) ===", n_perm))
    perm_result <- permutation_threshold(
      harm$cleaned,
      n_perm = n_perm,
      seed = params$perm_seed %||% NULL,
      trace = !quiet
    )
  } else {
    msg("=== Step 6: Permutation skipped ===")
  }

  # ---- Step 7: Plots ----
  if (run_plots) {
    msg("=== Step 7: Generate plots ===")
    tryCatch({
      grDevices::png(file.path(figure_dir, "manhattan_interaction.png"),
                     width = 1500, height = 800, res = 150)
      plot_manhattan(screen$results, p_col = "InteractionP",
                     threshold = params$screen_threshold %||% 5e-8,
                     title = sprintf("Interaction Manhattan (%s)", design$design_type))
      grDevices::dev.off()
    }, error = function(e) warning("Manhattan plot failed: ", e$message))

    tryCatch({
      grDevices::png(file.path(figure_dir, "qq_interaction.png"),
                     width = 1000, height = 1000, res = 150)
      plot_qq(screen$results, p_col = "InteractionP")
      grDevices::dev.off()
    }, error = function(e) warning("QQ plot failed: ", e$message))

    tryCatch({
      grDevices::png(file.path(figure_dir, "scatter_interaction.png"),
                     width = 1400, height = 1000, res = 150)
      plot_scatter(screen$results, imrp_result = screen$imrp,
                   title = sprintf("Effect B vs Effect A (%s)", design$design_type))
      grDevices::dev.off()
    }, error = function(e) warning("Scatter plot failed: ", e$message))
  }

  # ---- Step 8: Report ----
  if (run_report) {
    msg("=== Step 8: Generate report ===")
    generate_report(
      screen_result = screen,
      harm_result = harm,
      iv_result = iv,
      rho_result = rho,
      perm_result = perm_result,
      output_file = file.path(output_dir, "INTERACTION_REPORT.md")
    )
  }

  result <- list(
    design = design,
    data = data,
    harmonized = harm,
    instruments = iv,
    rho = rho,
    screen = screen,
    permutation = perm_result,
    output_dir = output_dir,
    params = params
  )
  class(result) <- "mrgxe_pipeline"
  msg("=== Pipeline complete ===")
  result
}

#' @export
print.mrgxe_pipeline <- function(x, ...) {
  cat("MRGxE pipeline result\n")
  cat(sprintf("  Design:     %s\n", x$design$design_type))
  cat(sprintf("  Variants:   %d\n", x$harmonized$nrow_cleaned %||%
              nrow(x$harmonized$cleaned %||% x$harmonized$qc_counts)))
  cat(sprintf("  Theta:      %.4f\n", x$screen$imrp$causal_estimate))
  cat(sprintf("  Significant: %d\n", x$screen$n_significant))
  cat(sprintf("  Rho:        %.4f\n", x$rho$rho %||% 0))
  cat(sprintf("  Output:     %s\n", x$output_dir))
}
