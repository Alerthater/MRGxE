# MRGxE: Genome-to-Genome Interaction Screening via Mendelian Randomization

**Version 0.2.0**

MRGxE implements a generalized MR-based interaction screening framework for detecting pairwise interactions between two genomic entities. The core insight is that the relationship between two related sets of genetic effect estimates follows a linear regression under the null, and variants with interaction or mediation effects systematically deviate from this line — analogous to horizontal pleiotropic variants in MR analysis.

## Citation

If you use MRGxE, please cite:

- **GxE framework**: Zhu et al. (2024) *Nature Communications*. "A MR-based approach for gene-environment interaction screening." [doi:10.1038/s41467-024-47806-3](https://doi.org/10.1038/s41467-024-47806-3)
- **Genome-to-genome design**: Chen et al. (2026) *Nature*. "EBV strain interacts with host HLA to drive nasopharyngeal carcinoma risk."
- **MRGxE package**: Yang et al. (2026). "MRGxE: An R package for genome-to-genome interaction screening."

## Supported Study Designs

| Design | Description | Reference |
|--------|-------------|-----------|
| `gxe` | Gene-environment interaction | Zhu et al. 2024 |
| `genome_to_genome` | Host-pathogen, gene-gene interaction | Chen et al. 2026 |
| `custom` | Any pair of effect estimates | User-defined |

## Features

- **Study design abstraction**: Unified interface for GxE, genome-to-genome, and custom designs
- **IMRP-based screening**: Iterative MR and pleiotropy testing using the TMRGxE/T-interaction statistic
- **Multiple data modes**: Summary statistics or individual-level data
- **QC by design**: Automatic QC defaults tailored to each design type (standardization, palindromic removal, EAF checks, etc.)
- **Conditional analysis**: Test whether a lead variant explains a regional signal (e.g., HLA-A\*11:01 conditioning)
- **Permutation testing**: Empirical significance thresholds, including BRASS-style clustered permutation for structured samples
- **Visualization**: Manhattan, QQ, scatter, and conditional forest plots
- **Pipeline wrapper**: One-call end-to-end analysis with automated reporting
- **Backward compatible**: All v0.1.0 functions (`import_summary_stats`, `harmonize_gwas_gwis`, `run_tmrgxe`, etc.) continue to work

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

### v0.2.0 (2026-06)
- Added `study_design()` abstraction with `gxe`, `genome_to_genome`, `custom` presets
- Added `screen_interaction()` — generalized screening for any design
- Added `conditional_analysis()` — summary-statistics and individual-level modes
- Added `permutation_threshold()` and `brass_permutation()`
- Added `compute_effects_from_individual()` and `individual_gxe_scan()`
- All v0.1 functions maintained with backward compatibility

### v0.1.0 (2026-06)
- Initial release: GxE screening with IMRP
- Core functions: `import_summary_stats()`, `harmonize_gwas_gwis()`, `run_tmrgxe()`

## Dependencies

- **Required**: `IMRP` (from [github.com/XiaofengZhuCase/IMRP](https://github.com/XiaofengZhuCase/IMRP))
- **Suggested**: `ggplot2`, `bigsnpr`, `data.table`, `testthat`, `knitr`, `rmarkdown`

## Contributing

Issues and pull requests are welcome. When reporting bugs, please include a minimal reproducible example.

## License

GPL-3 © Yang et al.
