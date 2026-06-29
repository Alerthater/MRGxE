#' Simulate GxE Summary Statistics
#'
#' Generates synthetic GWAS and GWIS summary statistics with known GxE
#' interaction variants, for testing and demonstration purposes.
#'
#' @param seed Integer. Random seed for reproducibility. Default \code{NULL}.
#' @param n_variants Integer. Number of variants to simulate. Default 1000.
#' @param n_palindromic Integer. Number of variants to make palindromic (A/T or C/G).
#'   Used for testing palindromic SNP filtering. Default 0.
#' @param true_theta Numeric. The true theta (regression slope) value. Default 0.72.
#' @param n_gxe Integer. Number of true GxE interaction variants. Default 50.
#' @param n_gwas Integer. GWAS sample size. Default 5000.
#' @param n_gwis Integer. GWIS sample size. Default 5000.
#'
#' @return A list with:
#'   \describe{
#'     \item{gwis}{Data frame of GWIS main effect summary statistics}
#'     \item{gwas}{Data frame of GWAS marginal effect summary statistics}
#'     \item{true_theta}{The true theta value}
#'     \item{true_gxe_snps}{Integer vector of indices of true GxE SNPs}
#'   }
#' @export
simulate_gxe_summary <- function(
  seed = NULL,
  n_variants = 1000,
  n_palindromic = 0,
  true_theta = 0.72,
  n_gxe = 50,
  n_gwas = 5000,
  n_gwis = 5000
) {
  if (!is.null(seed)) set.seed(seed)

  # Simulate variant positions across 22 autosomal chromosomes
  n_chr <- 22
  chr <- sample(1:n_chr, n_variants, replace = TRUE)
  bp <- unlist(lapply(table(chr), function(n) sort(sample(1:2e8, n))))
  snp <- paste0("rs", seq_len(n_variants))

  # Generate effect alleles (A1) and other alleles (A2)
  a1 <- sample(c("A", "C", "G", "T"), n_variants, replace = TRUE)
  a2 <- sample(c("A", "C", "G", "T"), n_variants, replace = TRUE)

  # Ensure A1 != A2
  same <- a1 == a2
  if (any(same)) {
    swap <- c(A = "C", C = "G", G = "T", T = "A")
    a2[same] <- swap[a1[same]]
  }

  # Make some palindromic (A/T, T/A, C/G, G/C) if requested
  if (n_palindromic > 0) {
    pal_idx <- sample(n_variants, min(n_palindromic, n_variants))
    pal_pairs <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
    for (i in pal_idx) {
      a2[i] <- pal_pairs[a1[i]]
    }
  }

  # Generate allele frequencies
  eaf <- runif(n_variants, 0.05, 0.95)

  # Generate GWIS main effects (beta_1) and GWAS marginal effects (alpha)
  # alpha = theta * beta_1 + GxE_deviation + noise

  # GWIS beta: some null (near 0), some with main effect signal
  beta_b <- rnorm(n_variants, mean = 0, sd = 0.1)

  # Make some variants have stronger main effects (for instrument selection)
  n_iv_candidates <- 100
  iv_idx <- sample(n_variants, n_iv_candidates)
  beta_b[iv_idx] <- rnorm(n_iv_candidates, mean = 0.15, sd = 0.04)

  se_b <- 0.02 + runif(n_variants, 0, 0.01)

  # GWAS effects
  # GxE deviation for interaction variants
  gxe_dev <- rep(0, n_variants)
  true_gxe_idx <- sample(setdiff(seq_len(n_variants), iv_idx), n_gxe)
  gxe_dev[true_gxe_idx] <- rnorm(n_gxe, mean = 0.06, sd = 0.03)

  beta_a <- true_theta * beta_b + gxe_dev + rnorm(n_variants, 0, 0.02)
  se_a <- 0.025 + runif(n_variants, 0, 0.01)

  # Compute P-values (only GWAS has P for instrument selection)
  p_gwas <- 2 * pnorm(-abs(beta_a / se_a))

  # Ensure palindromic SNPs get flagged properly: assign consistent alleles
  if (n_palindromic > 0) {
    # already handled above
  }

  # Build GWIS data frame
  gwis <- data.frame(
    SNP  = snp,
    CHR  = chr,
    BP   = bp,
    A1   = a1,
    A2   = a2,
    BETA = beta_b,
    SE   = se_b,
    N    = n_gwis,
    EAF  = eaf,
    stringsAsFactors = FALSE
  )

  # Build GWAS data frame (includes P column for instrument selection)
  gwas <- data.frame(
    SNP  = snp,
    CHR  = chr,
    BP   = bp,
    A1   = a1,
    A2   = a2,
    BETA = beta_a,
    SE   = se_a,
    N    = n_gwas,
    EAF  = eaf,
    P    = p_gwas,
    stringsAsFactors = FALSE
  )

  list(
    gwis = gwis,
    gwas = gwas,
    true_theta = true_theta,
    true_gxe_snps = true_gxe_idx
  )
}

#' Create Example Data (Alias)
#'
#' Alias for \code{\link{simulate_gxe_summary}}. Provided for backward
#' compatibility and readability.
#'
#' @inheritParams simulate_gxe_summary
#'
#' @return Same as \code{\link{simulate_gxe_summary}}.
#' @export
create_example_data <- function(
  seed = NULL,
  n_variants = 1000,
  n_palindromic = 0,
  true_theta = 0.72,
  n_gxe = 50,
  n_gwas = 5000,
  n_gwis = 5000
) {
  simulate_gxe_summary(
    seed = seed,
    n_variants = n_variants,
    n_palindromic = n_palindromic,
    true_theta = true_theta,
    n_gxe = n_gxe,
    n_gwas = n_gwas,
    n_gwis = n_gwis
  )
}

#' Estimate GxE Interaction Heritability (Alias)
#'
#' Alias for \code{\link{estimate_gxe_heritability}}. Provided for compatibility
#' with external code that expects the \code{estimate_interaction_heritability}
#' name.
#'
#' @inheritParams estimate_gxe_heritability
#'
#' @return Same as \code{\link{estimate_gxe_heritability}}.
#' @export
estimate_interaction_heritability <- function(
  tmrgxe_result,
  method = c("bigsnpr", "ldscR"),
  ld_ref_path = NULL,
  population = "EUR",
  output_dir = NULL
) {
  estimate_gxe_heritability(
    tmrgxe_result = tmrgxe_result,
    method = method,
    ld_ref_path = ld_ref_path,
    population = population,
    output_dir = output_dir
  )
}
