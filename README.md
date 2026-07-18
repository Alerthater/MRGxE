# MRGxE: Genome-to-Genome Interaction Screening via Mendelian Randomization

**Version 0.3.0**

MRGxE implements a generalized MR-based interaction screening framework for detecting pairwise interactions between two genomic entities. The core insight is that the relationship between two related sets of genetic effect estimates follows a linear regression under the null, and variants with interaction or mediation effects systematically deviate from this line — analogous to horizontal pleiotropic variants in MR analysis.

## Citation

If you use MRGxE, please cite:

- **GxE framework**: Zhu et al. (2024) *Nature Communications*. "A MR-based approach for gene-environment interaction screening." [doi:10.1038/s41467-024-47806-3](https://doi.org/10.1038/s41467-024-47806-3)
- **Genome-to-genome design**: Chen et al. (2026) *Nature*. "EBV strain interacts with host HLA to drive nasopharyngeal carcinoma risk."

## Supported Study Designs

| Design | Description | Reference |
|--------|-------------|-----------|
| `gxe` | Gene-environment interaction | Zhu et al. 2024 |
| `genome_to_genome` | Host-pathogen, gene-gene interaction | Chen et al. 2026 |
| `custom` | Any pair of effect estimates | User-defined |

## Features

- **IMRP-based screening**: Interaction deviation testing using the IMRP framework
- **Multiple data modes**: Summary statistics or individual-level data
- **Study design abstraction**: Unified interface for GxE, genome-to-genome, and custom designs
- **QC by design**: Automatic QC defaults tailored to each design type (standardization, palindromic removal, EAF checks, etc.)
- **Conditional analysis**: Test whether a lead variant explains a regional signal (e.g., HLA-A*11:01 conditioning)
- **Permutation testing**: Empirical significance thresholds, including BRASS-style clustered permutation
- **Individual-level simulation**: Reproduce Zhu et al. (2024) simulation methodology (Fig 2A–F)
- **Visualization**: Manhattan, QQ, scatter, and conditional plots
- **Pipeline wrapper**: One-call end-to-end analysis with automated reporting

## Installation

```r
# 1. Install IMRP (required dependency)
remotes::install_github("XiaofengZhuCase/IMRP")

# 2. Install MRGxE from GitHub
remotes::install_github("Alerthater/MRGxE")
```

## Quick Start

### GxE Screening (Zhu et al. 2024)

```r
library(MRGxE)

# Step 1: Import and harmonize effects
data <- import_effects("gwas.txt", "gwis.txt", design = gxe_design())
harm <- harmonize_effects(data, design = gxe_design())

# Step 2: Select instruments
iv <- select_instruments(harm$cleaned, p_threshold = 5e-8)

# Step 3: Estimate rho (sample overlap)
rho <- estimate_rho(harm$cleaned)

# Step 4: Run genome-wide interaction screening
screen <- screen_interaction(harm$cleaned, iv$instruments, rho$rho, 
                             design = gxe_design())

# Step 5: Visualize
plot_manhattan(screen$results)
plot_qq(screen$results)
plot_scatter(screen$results, imrp_result = screen$imrp)
```

### Genome-to-Genome Screening (Chen et al. 2026)

```r
library(MRGxE)

# Use genome-to-genome design (no rho needed, no standardization)
design_hp <- study_design("genome_to_genome",
  label_a = "EBV overall effect on NPC",
  label_b = "EBV interaction-independent effect")

data <- import_effects("ebv_overall.txt", "ebv_independent.txt", 
                       design = design_hp)
harm <- harmonize_effects(data, design = design_hp)

# Looser threshold for viral genome (fewer variants)
iv <- select_instruments(harm$cleaned, p_threshold = 5e-6)

# rho = 0 for independent genomes
screen <- screen_interaction(harm$cleaned, iv$instruments, rho = 0,
                             standardize = FALSE, design = design_hp)

# Conditional analysis (e.g., conditioning on lead SNP)
ca <- conditional_analysis(screen, condition_snp = "EBV_85841",
                           label = "EBV SNP 85841")
plot_conditional(ca)
```

### One-Call Pipeline

```r
# GxE example
result <- run_pipeline(
  effect_a_path = "gwas.txt",
  effect_b_path = "gwis.txt",
  design_type = "gxe",
  n_perm = 1000           # Optional: permutation threshold
)

# Genome-to-genome example
result <- run_pipeline(
  effect_a_path = "host_effects.txt",
  effect_b_path = "pathogen_effects.txt",
  design_type = "genome_to_genome",
  run_plots = TRUE
)
```

## Methodology Overview

### The Regression Framework

For two sets of effect estimates (e.g., GWAS marginal effect α̂ and GWIS main effect β̂₁), the relationship is:

$$\hat{\alpha} = \theta \hat{\beta}_1 + \text{deviation}$$

- **θ (theta)**: Estimated via IMRP (IVW after pleiotropy removal). In GxE designs, θ converges to 1 when GWAS and GWIS are identically conducted without interaction effects. In genome-to-genome designs, θ represents the proportion of overall effect attributable to interaction-independent pathways.
- **Deviation**: Variants with significant deviation from regression (TMRGxE test) are candidates for interaction or mediation.
- **ρ (rho)**: Correlation between the two effect estimates under the null, accounting for sample overlap (required for GxE, set to 0 for independent genomes).

### Two-Step Validation (Zhu et al. 2024)

1. **Step 1**: Genome-wide TMRGxE screening (P < 5×10⁻⁸)
2. **Step 2**: Bonferroni-corrected TDirect validation of candidates

### Conditional Analysis (Chen et al. 2026)

After identifying a significant interaction region, condition on the lead variant to determine whether it explains the regional signal — performed in either summary-statistics mode or individual-level mode.

## Comparison: GxE vs. Genome-to-Genome

| Parameter | GxE | Genome-to-Genome |
|-----------|-----|------------------|
| Standardization | β / SE / √N | Raw effects |
| rho (sample overlap) | Required | Not needed (ρ = 0) |
| Palindromic removal | Yes | No |
| EAF QC | Yes | No |
| IV P threshold | 5×10⁻⁸ | 5×10⁻⁶ |
| Output from IMRP | TMRGxE P | T-interaction P |

## Changelog

### v0.3.0 (2026-07)
- Added individual-level GxE simulation functions reproducing Zhu et al. (2024):
  `simulate_gxe_fig2ad_exact()`, `simulate_gxe_fig2ef_exact()`, `simulate_gxe_fig2ad()`,
  `simulate_gxe_fig2ef()`, `compute_gxe_effects_individual()`, `generate_gxe_phenotype()`
- Added integration bridge: `run_full_simulation_pipeline()`, `validate_simulation()`,
  `convert_to_mrgxe_format()`
- Added simulation visualization: `plot_fig2ad()`, `plot_fig2ef_power()`
- Added `run_single_scenario()` for individual Fig 2E-F scenario runs

## Dependencies

- **Required**: `IMRP` (from [github.com/XiaofengZhuCase/IMRP](https://github.com/XiaofengZhuCase/IMRP))
- **Suggested**: `ggplot2`, `bigsnpr`, `data.table`, `testthat`, `knitr`, `rmarkdown`

## Author

**Alerthater** — sole developer of this package. All core development was
performed independently, including:

- Architecture design and module interface definitions
- Core algorithm implementation (IMRP integration, interaction screening,
  conditional analysis, permutation testing, QC pipelines)
- Statistical verification and simulation testing against published results
- Documentation (roxygen, vignettes, README)

## Development & AI Acknowledgements

This package was developed by Alerthater. Core methodological design,
architecture decisions, statistical verification, and QC logic were implemented
and validated independently.

Claude Code was used to accelerate boilerplate generation (roxygen docs,
vignette scaffolding, plot templates) under the author's direction and review.

## Contributing

Issues and pull requests are welcome. When reporting bugs, please include a minimal reproducible example.

## License

GPL-3 © Alerthater
