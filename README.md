# VMS Speed Profile Mixture Model

Bayesian hierarchical mixture model for classifying fishing activity from Vessel Monitoring System (VMS) speed data, using `brms` in R.

## Overview

VMS data records vessel positions and speeds at regular intervals. Different activities — harbour/moored, slow manoeuvring, fishing, and steaming — produce distinct speed profiles. This script fits a **4-component Gamma mixture model** with hierarchical (random) effects across fishing métiers, allowing speed thresholds to vary by gear type and method.

The model estimates:

- **Component means (μ₁–μ₄):** typical speeds for each activity class, with random effects per métier (`LE_L5MET`)
- **Component shapes (shape₁–shape₄):** spread of each speed distribution (fixed across métiers)
- **Mixing weights (θ₁–θ₄):** proportion of observations in each component
- **Speed thresholds:** intersection points between adjacent components, used to classify fishing vs non-fishing activity

## Requirements

- **R** ≥ 4.1
- **R packages:**
  - `brms` — Bayesian regression modelling via Stan
  - `dplyr`, `tidyr`, `purrr` — data manipulation
  - `tidybayes` — posterior summaries
  - `ggplot2`, `patchwork` — plotting
- **Optional:** `cmdstanr` for faster sampling (requires [CmdStan](https://mc-stan.org/cmdstanr/) installation)

Install packages:

```r
install.packages(c("brms", "dplyr", "tidyr", "purrr", "tidybayes", "ggplot2", "patchwork"))
```

## Input Data

The script expects a file `Data/tacsatEflalo2024.RData` containing a data frame `tacsatEflalo` with at minimum:

| Column | Description |
|---|---|
| `SI_SP` | Vessel speed (knots) |
| `LE_GEAR` | Gear type code (Level 4, e.g. `OTB`, `DRB`) |
| `LE_L5MET` | Métier code (Level 5, e.g. `OTB_DEF_>=70_0`) |

## Usage

1. Place your `tacsatEflalo2024.RData` file in a `Data/` subdirectory.
2. Open `vms_mixture_model.R` and adjust settings as needed:
   - **Subsampling:** Uncomment the `slice_sample` block for faster development runs.
   - **Backend:** Set `backend = "cmdstanr"` if CmdStan is installed and working, otherwise leave it as `rstan`.
   - **Number of components (K):** Default is 4. Reduce to 3 if a simpler model is more appropriate for your data.
3. Run the script.

## Model Structure

```
SI_SP ~ Gamma mixture (K=4, ordered)

Component means:    mu_k ~ Intercept + (1 | LE_L5MET)    for k = 1..4
Component shapes:   shape_k ~ Intercept                   for k = 1..4
```

The `order = TRUE` constraint ensures μ₁ < μ₂ < μ₃ < μ₄, so components map consistently to increasing speed activities.

Random effects on means allow each métier to have its own characteristic speeds, while shape parameters are shared globally to keep the model tractable.

## Outputs

- **Fitted model object** cached as `brms_vms_mix_k4_gamma_simplified.rds`
- **Threshold table** printed to console — speed boundaries between activity classes per métier
- **Diagnostic plots** — trace plots, posterior predictive checks
- **Threshold plots** saved as `Threshold_Plot_{GEAR}.png` — histograms of observed speeds overlaid with fitted mixture components and threshold lines

## Performance Notes

This script is optimised for speed compared to a fully nested model with random effects on all parameters. Key simplifications:

- Random effects on component means only (not shapes)
- Single grouping level (`LE_L5MET`) rather than nested (`LE_GEAR/LE_L5MET`)
- Optional subsampling for model development

For large datasets (>1M rows), consider subsampling during development and running the full data overnight.

## Licence

This code is provided as-is for research purposes.
