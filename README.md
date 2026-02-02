# CoxFP: Resampling-Based Multiple Testing Adjustment in Cox Models with Fractional Polynomial Transformations


**Reference:** Liquet, B., Roux, M., & Riou, J. (2026). Resampling-Based Multiple Testing Adjustment for Fractional Polynomial Cox Models.

![Made with R](https://img.shields.io/badge/Made%20with-R-276DC3.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)
![Status](https://img.shields.io/badge/Status-Research%20Prototype-blue.svg)
---


### What are Fractional Polynomials?

Fractional polynomials (FP) are a flexible class of parametric functions used to model non-linear relationships between a continuous predictor and an outcome. Unlike standard polynomials ($x$, $x^²$, $x^³$), FPs use powers from a restricted set: {-2, -1, -0.5, 0, 0.5, 1, 2, 3}.

**Key advantages:**
- More flexible than linear models
- More interpretable than splines
- Can capture various non-linear shapes (U-shaped, J-shaped, monotonic, etc.)

### Why Multiple Testing Correction?

When testing multiple FP transformations (e.g., FP1 with powers -2, -1, 0, 1, 2), we perform multiple hypothesis tests. Without correction, the probability of false positives (Type I error) inflates dramatically. This package implements several correction methods to control the Family-Wise Error Rate (FWER).
This repository implements several multiple‑testing correction strategies:
- **Bonferroni correction**
- **Royston–Sauerbrei MFP procedure**
- **Parametric bootstrap under H₀**
- **Semi‑parametric bootstrap (weighted Cox model)**
- **Wild bootstrap (martingale perturbation)**
- **Permutation of martingale residuals**
  
These methods allow valid inference when selecting the functional form of a continuous covariate in survival models.

### Typical Use Cases

- **Biomarker research**: Finding optimal functional form for biomarker-outcome relationships
- **Epidemiology**: Modeling dose-response curves
- **Clinical research**: Identifying non-linear risk factors in survival analysis
- **Pharmacology**: Modeling drug concentration-effect relationships



---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Function Arguments](#function-arguments)
- [Correction Methods](#correction-methods)
- [Output Structure](#output-structure)
- [Interpretation Guide](#interpretation-guide)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Installation

```r
# Install dependencies
if (!require("survival")) install.packages("survival")

# Clone and source
git clone https://github.com/JeremieRiou/Multiple_Testing_MFP.git
source("cox_fp_corrected.R")
```

**Requirements:** R > 3.6.0, `survival` package

---

## Quick Start

```r
library(survival)
source("cox_fp_functions.R")

lung
lung$status <- as.numeric(paste(ifelse(lung$status=="1","0","1")))

# Define FP transformations to test
fp_powers <- list(
  c(-2), c(-1), c(0), c(1), c(2),      # FP1 models
  c(-1, -1), c(0, 1), c(1, 2)          # FP2 models
)

# Run analysis
result <- cox_fp_analysis(
  formula = Surv(time, status) ~ age + sex,
  data = lung,
  var_interest = "ph.ecog",
  fp_powers = fp_powers,
  method = "bonferroni",
  alpha = 0.05
)

# View results
print(result$p_corrected)    # Adjusted p-value
print(result$best_power)     # Selected transformation
```

---

## Function Arguments

### `cox_fp_analysis()`

| Argument | Type | Default | Description |
|----------|------|---------|-------------|
| `formula` | formula | - | `Surv(time, status) ~ covariates` (exclude `var_interest`) |
| `data` | data.frame | - | Dataset with all variables |
| `var_interest` | character | - | Variable to transform with FP |
| `fp_powers` | list | - | List of FP transformations to test |
| `method` | character | `"bonferroni"` | Correction method (see below) |
| `B` | integer | `500` | Bootstrap replications (resampling methods) |
| `alpha` | numeric | `0.05` | Significance level |
| `offset` | numeric | `0.1` | Offset for positivity |
| `verbose` | logical | `TRUE` | Display progress |

### FP Powers Examples

```r
# FP1 (single power)
c(-2), c(-1), c(0), c(1), c(2)

# FP2 (two powers)
c(-1, -1)   # Repeated: x^(-1) and x^(-1)*log(x)
c(0, 1)     # Different: log(x) and x
c(1, 2)     # Different: x and x^2
```

---

## Correction Methods

| Method | String | Description | Use When | Speed |
|--------|--------|-------------|----------|-------|
| **Naive** | `"naive"` | No correction | Not recommended | Instant |
| **Bonferroni** | `"bonferroni"` | Conservative FWER control | Quick analysis, any n | Instant |
| **Royston-Sauerbrei** | `"royston"` | Sequential testing (MFP) | FP1+FP2, standard analysis | Fast |
| **Parametric Bootstrap** | `"parametric_bootstrap"` | Model-based resampling | n > 100 | Moderate |
| **Semi-parametric Bootstrap** | `"semiparametric_bootstrap"` | Weighted resampling | Robust, n > 100 | Moderate |
| **Wild Bootstrap** | `"wild_bootstrap"` | Martingale perturbation | Small n, uncertain assumptions | Moderate |
| **Martingale Permutation** | `"perm_martingale"` | Residual permutation | Non-parametric | Moderate |

**Recommendation:**
- **Quick/exploratory:** `"bonferroni"`
- **Standard analysis:** `"royston"` (requires FP1 and FP2)
- **Robust inference:** `"parametric_bootstrap"` with B=1000

---

## Output Structure

```r
result <- cox_fp_analysis(...)

# Key results
result$p_corrected          # Adjusted p-value
result$best_pval_raw        # Raw p-value of best model
result$best_power           # Selected FP powers
result$best_fit             # coxph object
result$best_aic             # AIC

# All models
result$all_pvalues_raw      # All raw p-values
result$all_results          # Full details per model

# Metadata
result$correction_method    # Method used
result$n_obs                # Sample size
result$n_tests              # Number of FP models (K)
result$additional_info      # Method-specific info (e.g., bootstrap success rate)
```

---

## Interpretation Guide

### Decision Rule

```r
if (result$p_corrected < alpha) {
  # Significant: Use result$best_power transformation
} else {
  # Not significant: Consider linear model or exclusion
}
```

### Common FP Transformations

| Powers | Shape | Interpretation |
|--------|-------|----------------|
| `c(1)` | Linear | Proportional effect |
| `c(0)` | Log | Diminishing effect |
| `c(2)` | Quadratic | U-shaped |
| `c(-1)` | Inverse | Asymptotic decay |
| `c(0, 1)` | Log-linear | Flexible monotonic |

### Example Interpretation

```r
result$best_power = c(2)
result$p_corrected = 0.00020
```
*"A log-linear FP2 transformation provides the best fit (adjusted p = 0.008), suggesting the effect plateaus at higher values."*

---

## Troubleshooting

### Common Errors

**"Sample too small"**
- Need: n > 50 complete cases, > 10 events
- Solution: Check with `table(complete.cases(data))` and `table(data$status)`

**"No FP model could be adjusted"**
- Cause: Convergence failures
- Solutions:
  ```r
  # Increase offset
  result <- cox_fp_analysis(..., offset = 1.0)
  
  # Reduce FP models
  fp_powers <- list(c(-1), c(0), c(1), c(2))
  ```

**"Too few valid bootstraps"**
- Cause: Bootstrap failures
- Solutions:
  ```r
  # Increase B
  result <- cox_fp_analysis(..., B = 1000)
  
  # Or switch method
  result <- cox_fp_analysis(..., method = "bonferroni")
  ```

### Data Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| Sample size | 50 | 100+ |
| Events | 10 | 30+ |
| Missing data | Complete cases | < 10% missing |
| Variable range | Positive (offset applied) | Wide range |

---

## Limitations

- The computational cost may be high for large datasets or large FP grids, especially with bootstrap-based methods.
- The bootstrap procedures require a sufficient number of events; results may be unstable when the number of events is very small.
- The variable of interest must be strictly positive for FP transformations; an offset is added automatically when needed.
- The method does not currently support time‑dependent covariates or stratified Cox models.
- Only one continuous variable is tested at a time; multivariable FP selection is not implemented.


## Detailed Example

See [Motivating_Example.md](https://github.com/JeremieRiou/Multiple_Testing_MFP/blob/main/Motivating_Example.md) for a comprehensive real-world application.

**Quick comparison of methods:**

```r
# Conservative
result_bonf <- cox_fp_analysis(..., method = "bonferroni")

# Standard
result_roy <- cox_fp_analysis(..., method = "royston")

# Robust (slower)
result_boot <- cox_fp_analysis(..., method = "parametric_bootstrap", B = 1000)

# Compare
data.frame(
  Method = c("Bonferroni", "Royston", "Bootstrap"),
  P_adj = c(result_bonf$p_corrected, result_roy$p_corrected, result_boot$p_corrected)
)
```

---

## Citation

```bibtex
@article{liquet2026coxfp,
  title={Resampling-Based Multiple Testing Adjustment for Fractional Polynomial Cox Models},
  author={Liquet, Benoit and Roux, Mathieu and Riou, J{\'e}r{\'e}mie},
  year={2026},
  note={In preparation}
}
```

---

## Additional Resources

- [Motivating Example](https://github.com/JeremieRiou/Multiple_Testing_MFP/blob/main/Motivating_Example.md)
- Royston & Sauerbrei (2008) - Original MFP methodology

---

## Important Notes

- **Variable positivity:** `var_interest` must be positive (offset applied automatically if needed)
- **Power 0:** Represents log-transformation by FP convention
- **Repeated powers:** `c(1,1)` creates x and xÃlog(x) terms
- **Report adjusted p-values:** Always use `p_corrected` to maintain valid Type I error control
- **Computational time:** Bootstrap methods take 1-5 min depending on n and K

---

**Contact:** [GitHub Issues](https://github.com/JeremieRiou/Multiple_Testing_MFP/issues)
