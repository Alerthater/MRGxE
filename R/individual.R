#' Compute Effects from Individual-Level Data
#'
#' For studies with individual-level data (not summary statistics), this function
#' computes the genetic effects needed for interaction screening directly.
#' Supports both linear and logistic regression models.
#'
#' @param phenotype Numeric vector or data frame column name. Trait/outcome values.
#' @param genotype_a Matrix or data frame. Genotypes for entity A (e.g., host SNPs).
#' @param genotype_b Vector or data frame column. Genotype for entity B
#'   (e.g., EBV strain indicator or single variant).
#' @param covariate_matrix Matrix or data frame. Covariates (age, sex, PCs, etc.).
#' @param family Character. \code{"gaussian"} for quantitative traits,
#'   \code{"binomial"} for binary traits. Default \code{"gaussian"}.
#' @param model Character. \code{"marginal"} for GWAS-style marginal effect (Y ~ G + covariates),
#'   or \code{"interaction"} for full model with interaction term.
#'
#' @return Data frame with effect estimates (beta, se, P) for each variant in genotype_a.
#'
#' @export
compute_effects_from_individual <- function(
  phenotype,
  genotype_a,
  genotype_b = NULL,
  covariate_matrix = NULL,
  family = c("gaussian", "binomial"),
  model = c("marginal", "interaction")
) {
  family <- match.arg(family)
  model <- match.arg(model)

  n <- length(phenotype)
  if (is.matrix(genotype_a)) {
    p <- ncol(genotype_a)
    var_names <- colnames(genotype_a) %||% paste0("VAR", seq_len(p))
  } else {
    p <- 1
    var_names <- "VAR1"
    genotype_a <- matrix(genotype_a, ncol = 1)
  }

  # Build base formula components
  base_vars <- if (!is.null(covariate_matrix)) {
    as.data.frame(covariate_matrix)
  } else {
    data.frame(intercept = rep(1, n))
  }

  results <- data.frame(
    variant = var_names,
    beta = numeric(p),
    se = numeric(p),
    p = numeric(p),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(p)) {
    g_i <- genotype_a[, i]

    if (model == "marginal") {
      # Y ~ G + covariates
      model_data <- cbind(
        data.frame(y = phenotype, g = g_i),
        base_vars
      )
      fit <- stats::glm(y ~ ., data = model_data, family = family)
      results$beta[i] <- stats::coef(fit)["g"]
      results$se[i] <- sqrt(diag(stats::vcov(fit)))["g"]
      results$p[i] <- 2 * stats::pnorm(-abs(results$beta[i] / results$se[i]))
    } else if (model == "interaction" && !is.null(genotype_b)) {
      # Y ~ G + E + G:E + covariates
      model_data <- cbind(
        data.frame(y = phenotype, g = g_i, e = genotype_b),
        base_vars
      )
      fit <- stats::glm(y ~ g * e + ., data = model_data, family = family)
      # Return the interaction coefficient
      coef_name <- "g:e"
      if (!coef_name %in% names(stats::coef(fit))) {
        # Try alternative naming
        coef_name <- grep(":e$", names(stats::coef(fit)), value = TRUE)[1] %||% "g:e"
      }
      results$beta[i] <- stats::coef(fit)[coef_name]
      results$se[i] <- sqrt(diag(stats::vcov(fit)))[coef_name]
      results$p[i] <- 2 * stats::pnorm(-abs(results$beta[i] / results$se[i]))
    }
  }

  results
}

#' Individual-Level GxE Scan
#'
#' Performs a full genome-wide interaction scan using individual-level data.
#' Computes marginal and main effects for each variant, then runs the
#' MR-based interaction screening.
#'
#' @param phenotype Numeric vector. Trait values.
#' @param genotype_matrix Matrix (N x P). Genotypes for all P variants.
#' @param environment Vector. Environmental variable.
#' @param covariates Matrix. Covariates.
#' @param family Character. Model family.
#' @param ... Additional arguments passed to \code{\link{screen_interaction}}.
#'
#' @return \code{mrgxe_screen} object.
#'
#' @export
individual_gxe_scan <- function(
  phenotype,
  genotype_matrix,
  environment,
  covariates = NULL,
  family = "gaussian",
  ...
) {
  n <- length(phenotype)
  p <- ncol(genotype_matrix)
  var_names <- colnames(genotype_matrix) %||% paste0("rs", seq_len(p))

  message(sprintf("Computing marginal effects for %d variants...", p))

  # Compute marginal effects (GWAS-style, Y ~ G + covariates)
  marginal <- compute_effects_from_individual(
    phenotype, genotype_matrix, NULL, covariates,
    family = family, model = "marginal"
  )

  message("Computing main effects (with E in model)...")

  # Compute main effects (GWIS-style, Y ~ G + E + covariates)
  main_effects <- numeric(p)
  main_se <- numeric(p)
  for (i in seq_len(p)) {
    fit <- stats::glm(
      phenotype ~ genotype_matrix[, i] + environment + .,
      data = if (!is.null(covariates)) {
        cbind(data.frame(covariates))
      } else {
        data.frame(dummy = rep(1, n))
      },
      family = family
    )
    # Extract genetic main effect coefficient
    coef_idx <- 2  # second coefficient (first is intercept, second is G)
    main_effects[i] <- stats::coef(fit)[coef_idx]
    main_se[i] <- sqrt(diag(stats::vcov(fit)))[coef_idx]
  }

  # Build summary-statistics-like data frame
  data <- data.frame(
    id = var_names,
    CHR = rep(1, p), BP = seq_len(p),
    beta_a = marginal$beta, se_a = marginal$se, n_a = n,
    beta_b = main_effects, se_b = main_se, n_b = n,
    P = marginal$p,
    stringsAsFactors = FALSE
  )

  # Run interaction screening
  screen_interaction(data = data, ..., rho = 0)
}
