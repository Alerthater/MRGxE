#' MRGxE: Generalized MR-Based Interaction Screening Framework
#'
#' Detects interactions between two genomic entities by identifying variants
#' that deviate from the regression line between two sets of genetic effect
#' estimates. This generalizes the original GxE framework (Zhu et al. 2024) to
#' any pairwise study design (host-pathogen, gene-gene, etc.), as demonstrated
#' in Chen et al. (2026) for EBV-HLA interaction in nasopharyngeal carcinoma.
#'
#' @section Core concept:
#' Given two sets of effect estimates for the same set of variants:
#' \deqn{\text{Effect}_A = \theta \times \text{Effect}_B + \text{deviation}}
#'
#' Most variants follow the regression line; variants with interaction or
#' mediation effects systematically deviate. The deviation is tested via a
#' pleiotropy-like statistic analogous to IMRP.
#'
#' @section Study design abstraction:
#' The package provides a \code{\link{study_design}} object that captures the
#' semantics of each analysis, enabling the same pipeline code to handle:
#' \itemize{
#'   \item \strong{GxE design} (Zhu 2024): Effect_A = GWAS marginal, Effect_B = GWIS main
#'   \item \strong{Genome-to-genome design} (Chen 2026): Effect_A = EBV variant overall effect
#'     on NPC, Effect_B = EBV-HLA interaction-independent effect
#'   \item \strong{Custom design}: any pair of numeric effect vectors
#' }
#'
#' @section Main workflow:
#' \enumerate{
#'   \item Define study design with \code{\link{study_design}()} or presets
#'   \item Import and harmonize effects
#'   \item Select independent instruments and estimate rho
#'   \item Estimate theta via IMRP
#'   \item Screen with interaction deviation test
#'   \item Validate (two-step, conditional analysis, permutation)
#'   \item Visualize and report
#' }
#'
#' @references
#' Zhu X, Yang Y, Lorincz-Comi N, et al. (2024).
#' An approach to identify gene-environment interactions and reveal new
#' biological insight in complex traits. \emph{Nature Communications}, 15:3385.
#'
#' Chen et al. (2026). EBV strain interacts with host HLA to drive nasopharyngeal
#' carcinoma risk. \emph{Nature}.
#'
#' @docType package
#' @name MRGxE
NULL
