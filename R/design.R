#' Study Design Abstraction
#'
#' Defines the semantics of an interaction screening study. This is the central
#' abstraction that enables the same pipeline code to handle GxE (Zhu 2024),
#' host-pathogen (Chen 2026), and any other pairwise study design.
#'
#' A study design specifies:
#' - What the two effect estimates represent (label, source)
#' - How theta should be interpreted
#' - What QC steps are needed
#' - How the deviation statistic maps to the interaction of interest
#'
#' @param design_type Character. One of:
#'   \code{"gxe"} — gene-environment interaction (Zhu et al. 2024)
#'   \code{"genome_to_genome"} — host-pathogen interaction (Chen et al. 2026)
#'   \code{"custom"} — user-defined design
#' @param label_a Character. Label for the "outcome" effect (y-axis).
#' @param label_b Character. Label for the "exposure" effect (x-axis).
#' @param theta_interpretation Character. What theta represents.
#' @param requires_rho Logical. Whether rho (sample overlap correlation) is needed.
#'   Default TRUE for GxE; FALSE for genome-to-genome.
#' @param requires_individual_data Logical. Whether individual-level data is needed
#'   (as opposed to summary statistics).
#' @param interaction_meaning Character. Description of what the deviation from
#'   regression line means in this design context.
#' @param ... Additional design-specific parameters.
#'
#' @return A list of class \code{"mrgxe_design"} containing all design metadata.
#'
#' @examples
#' # GxE design (Zhu 2024)
#' design_gxe <- study_design("gxe")
#'
#' # Host-pathogen design (Chen 2026)
#' design_hp <- study_design("genome_to_genome",
#'   label_a = "EBV variant overall effect on NPC",
#'   label_b = "EBV-HLA interaction-independent effect")
#'
#' @export
study_design <- function(
  design_type = c("gxe", "genome_to_genome", "custom"),
  label_a = NULL,
  label_b = NULL,
  theta_interpretation = NULL,
  requires_rho = NULL,
  requires_individual_data = FALSE,
  interaction_meaning = NULL,
  ...
) {
  design_type <- match.arg(design_type)

  design <- switch(
    design_type,
    gxe = gxe_design(),
    genome_to_genome = genome_to_genome_design(),
    custom = custom_design(
      label_a = label_a %||% "Effect A",
      label_b = label_b %||% "Effect B",
      theta_interpretation = theta_interpretation %||%
        "Causal effect of exposure on outcome",
      requires_rho = requires_rho %||% TRUE,
      interaction_meaning = interaction_meaning %||%
        "Deviation from regression line"
    )
  )

  # Override with user-supplied values
  if (!is.null(label_a)) design$label_a <- label_a
  if (!is.null(label_b)) design$label_b <- label_b
  if (!is.null(theta_interpretation)) design$theta_interpretation <- theta_interpretation
  if (!is.null(requires_rho)) design$requires_rho <- requires_rho
  if (!is.null(interaction_meaning)) design$interaction_meaning <- interaction_meaning

  design$requires_individual_data <- requires_individual_data
  design$user_params <- list(...)
  class(design) <- "mrgxe_design"
  design
}

#' @rdname study_design
#' @export
gxe_design <- function() {
  list(
    design_type = "gxe",
    label_a = "GWAS marginal effect",
    label_b = "GWIS main effect",
    theta_interpretation = paste(
      "Contribution of GWIS main effect to GWAS marginal effect;",
      "converges to 1 when GWAS and GWIS are identically conducted"
    ),
    requires_rho = TRUE,
    requires_individual_data = FALSE,
    interaction_meaning = paste(
      "Deviation from GWAS~GWIS regression line;",
      "tests combined GxE and mediation effect"
    ),
    default_standardize = TRUE,
    default_iv_p_threshold = 5e-8,
    default_palindromic_removal = TRUE,
    default_eaf_qc = TRUE
  )
}

#' @rdname study_design
#' @export
genome_to_genome_design <- function() {
  list(
    design_type = "genome_to_genome",
    label_a = "Variant overall effect on trait (with interaction)",
    label_b = "Variant interaction-independent effect",
    theta_interpretation = paste(
      "Proportion of overall effect attributable to",
      "interaction-independent pathway"
    ),
    requires_rho = FALSE,
    requires_individual_data = FALSE,
    interaction_meaning = paste(
      "Deviation from regression line in the paired effect space;",
      "tests interaction between the two genomic entities"
    ),
    default_standardize = FALSE,
    default_iv_p_threshold = 5e-6,
    default_palindromic_removal = FALSE,
    default_eaf_qc = FALSE
  )
}

#' @rdname study_design
#' @export
custom_design <- function(
  label_a = "Effect A",
  label_b = "Effect B",
  theta_interpretation = "Causal effect",
  requires_rho = TRUE,
  interaction_meaning = "Interaction detected by deviation from regression"
) {
  list(
    design_type = "custom",
    label_a = label_a,
    label_b = label_b,
    theta_interpretation = theta_interpretation,
    requires_rho = requires_rho,
    requires_individual_data = FALSE,
    interaction_meaning = interaction_meaning,
    default_standardize = TRUE,
    default_iv_p_threshold = 5e-8,
    default_palindromic_removal = TRUE,
    default_eaf_qc = TRUE
  )
}

#' @export
print.mrgxe_design <- function(x, ...) {
  cat(sprintf("MRGxE study design: %s\n", x$design_type))
  cat(sprintf("  Effect A (y): %s\n", x$label_a))
  cat(sprintf("  Effect B (x): %s\n", x$label_b))
  cat(sprintf("  Theta:        %s\n", x$theta_interpretation))
  cat(sprintf("  Interaction:  %s\n", x$interaction_meaning))
  cat(sprintf("  Needs rho:    %s\n", x$requires_rho))
}
