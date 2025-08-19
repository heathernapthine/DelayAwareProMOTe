# ProMOTe

## Overview

This repository implements Variational Bayes (VB) methods for latent-class modelling of **long‑term conditions (LTCs)** using presence indicators and diagnosis times. It contains:
- **Simulation code** for three synthetic datasets (presence‑driven clusters, onset‑driven clusters, and onset‑driven with mixed diagnosis delays).
- **Two inference stacks**:
  1) a **delay‑aware** VB model that treats observed diagnosis time as *onset + delay* with condition‑level delay priors, and  
  2) a **baseline (no‑delay)** VB model that treats diagnosis times as onsets.
- **Prediction and evaluation scripts** (cluster recovery, forward prediction, ablations, identifiability checks) and plotting utilities.

> **Naming disclaimer (synthetic data).** LTC names are taken from UK Biobank phenotype lists purely as labels. In these simulations, the assignment of names to presence/onset/delay mechanisms is arbitrary and does **not** reflect real clinical patterns.


## Dependencies

Core packages used across generation, inference, prediction, and plotting:

- Data & maths: `dplyr`, `tidyr`, `purrr`, `stringr`, `MASS`, `extraDistr`, `truncnorm`
- VB utilities & clustering: `MCMCpack`, `clue`, `mclust`, `aricode`
- Metrics & ROC: `pROC`
- Plotting: `ggplot2`, `scales`, `patchwork`, `cowplot`, `ggridges`, `grid`
- I/O: `readr`

Install with:
```r
install.packages(c(
  "dplyr","tidyr","purrr","stringr","MASS","extraDistr","truncnorm",
  "MCMCpack","clue","mclust","aricode","pROC",
  "ggplot2","scales","patchwork","cowplot","ggridges","grid","readr"
))
```


## Data Structure

Each dataset is saved as a `.rds` containing a named list. Common fields:

- `d` (N×M matrix): Presence (1/0) of condition *m* for person *i*.
- `t` (N×M matrix): **Observed diagnosis time** (years). `NA` if absent/unknown.
- `rho` (length N): Baseline age (study entry) per person.
- `tau` (length N): Current age (end of follow‑up); typically `rho + 30`.
- `iota` (length N): Status at `tau` (encoding used by scripts; right‑censor logic in code).
- `N`, `K`, `M` (ints): Numbers of people, clusters, and conditions.
- `cond_list` (length M): Condition names (labels only; see disclaimer).
- Optional convenience fields (depending on generation script):
  - `onset` (N×M): Latent onset ages used to simulate `t`.
  - `delay` (N×M): Simulated delays (when a delay model was used).
  - `pi_mk` (M×K): True presence probabilities used in simulation.
  - `mu_mk`, `sigma2_mk` (M×K): True onset parameters per cluster.
  - `sex`, `birth_conds`, `male_conds`, `female_conds`: Masks for lifetime/sex‑specific conditions (may be `NULL`).
  - **Mixed‑delay datasets only:**  
    `delay_prior_df` (data frame with per‑condition delay family and parameters),  
    `delay_dist_cond`, `delay_mu_cond`, `delay_sd_cond` and mixture/ uniform parameters.

**Train/Test split.** Scripts follow an 80/20 split (≈160k/40k when `N=200000`) with `set.seed(42)` for reproducibility.


## Functions

Below we summarise the main exported functions grouped by model family. Each item lists **Purpose**, **Inputs**, and **Outputs**.

### A. Delay‑aware model (diagnosis = onset + delay)

#### 1. `VB_gaussian_update_d`
**Purpose:**  
Variational Bayes updates for the Gaussian latent‑class model with fixed `K`, left/right censoring, and **condition‑level delay priors** (mean `mu0`, variance `sigma20`, truncated at 0).

**Inputs:**
- `d`, `t_obs`, `rho`, `tau`, `iota`: See *Data Structure* (`t_obs` is the observed diagnosis matrix).
- `hyperparameters`: List of prior parameters as matrices/vectors:  
  `theta`, `a`, `b`, `u`, `v`, `alpha`, `beta`.
- Initial values: `initial_Cstar` (N×K), `initial_Dstar`, `initial_pstar`, `initial_qstar`, `initial_rstar` (N×M).
- Sizes & control: `N`, `M`, `K`, `epsilon`.
- **Delay priors:** `mu0` (length M), `sigma20` (length M).
- Optional masks: `sex`, `birth_conds`, `male_conds`, `female_conds`.
- `cond_list`: Condition names (for outputs).

**Outputs:**
- Posterior globals:  
  `theta_star`, `pi_a`/`pi_b` (or `a_star`/`b_star`), `mu_u` (`u_star`), `v_star`,  
  `mu_alpha`/`mu_beta` (`alpha_star`/`beta_star`).
- Patient‑level latents: `C_star`, `p_star`, `q_star`, `r_star`, `D_star`.
- Delay diagnostics: `gap_mu_star`, `gap_sigma2_star`.
- Expectations: `expected_t` (onset proxy inside window), `expected_d`.
- ELBO traces: `elbo`, `param_diffs`; iterations: `n_steps`.


#### 2. `VB_gaussian_predictive_density_d`
**Purpose:**  
Posterior predictive for a **new person** using a partial window `[rho, tau]`. Returns cluster weights and condition‑wise quantities, accounting for truncation/censoring and **delay** moments.

**Inputs:**
- `hyperparameters`: Posterior globals (`theta`, `a`, `b`, `u`, `v`, `alpha`, `beta`).
- Per‑person evidence split:  
  `M_obs` (fully observed in window), `M_part` (left‑censored), `M_unobs` (unobserved),  
  with `d_obs`, `t_obs`, `d_part`.
- Window: `rho`, `tau`; condition count: `M`.
- **Delay priors:** `mu0`, `sigma20` (length M).

**Outputs:**
- `phi` (length K): Posterior cluster probabilities.
- `eta` (M×K): `P(event after tau | cluster, condition)`.
- `varpi` (M×K): `P(present | cluster)`.


#### 3. `probability_LTHC_by_T_d`
**Purpose:**  
Compute, for each condition, `P(event in (tau, T] | evidence up to tau)` by mixing `eta` over clusters with `phi` and using the **delay‑aware** Student‑t predictive for observed times.

**Inputs:**
- `parameters`: Output of `VB_gaussian_predictive_density_d` (uses `phi`, `eta`).
- `hyperparameters`: Posterior globals (`u`, `v`, `alpha`, `beta`).
- Window: `T` (forecast horizon), `tau` (cut age), `M`.
- **Delay priors:** `mu0`, `sigma20`.

**Outputs:**
- Numeric vector (length M): probabilities by condition.


#### 4. `expected_LTHC_t_after_tau_d`
**Purpose:**  
Compute `E[Y | Y > tau]` (expected diagnosis time after `tau`) per condition, mixing over clusters with weights proportional to `phi .* eta`.

**Inputs:**
- `parameters`: Output of `VB_gaussian_predictive_density_d`.
- `hyperparameters`: Posterior globals (`u`, `v`, `alpha`, `beta`).
- `tau`: Cut age; `M`: number of conditions.
- **Delay priors:** `mu0`, `sigma20`.

**Outputs:**
- Numeric vector (length M): expected diagnosis times after `tau`.


#### 5. `expected_lst_lefttrunc_d`
**Purpose:**  
Expected value of a **left‑truncated** location‑scale Student‑t distribution; used by the delay‑aware predictive to get `E[Y | Y > tau]`.

**Inputs:** `df`, `mu`, `sigma`, `tau`.

**Outputs:** Numeric matrix/vector of expectations (shape conforms to inputs).


### B. Baseline model (no explicit delay)

#### 6. `VB_gaussian_update`
**Purpose:**  
Variational Bayes updates for the Gaussian latent‑class model with censoring **without** modelling delay (diagnosis ≈ onset in‑window).

**Inputs:**  
`d`, `t`, `rho`, `tau`, `iota`, `hyperparameters`, initial `*_star`, sizes `N/M/K`, `epsilon`, optional `sex`/`birth_conds`/`male_conds`/`female_conds`, `cond_list`.

**Outputs:**  
`theta_star`, `a_star`, `b_star`, `u_star`, `v_star`, `alpha_star`, `beta_star`,  
`C_star`, `p_star`, `q_star`, `r_star`, `D_star`, `expected_t`, `n_steps`, ELBO traces.


#### 7. `VB_gaussian_predictive_density`
**Purpose:**  
Posterior predictive for a new person on `[rho, tau]` ignoring delay.

**Inputs:**  
Posterior globals + per‑person `M_obs/M_part/M_unobs`, `d_obs/t_obs/d_part`, `rho`, `tau`, `M`.

**Outputs:**  
`phi` (clusters), `eta` (condition after‑tau probability per cluster), `varpi`.


#### 8. `probability_LTHC_by_T`
**Purpose:**  
Probability that each condition occurs in `(tau, T]` under the **baseline** predictive.

**Inputs:** `parameters`, `hyperparameters`, `T`, `tau`, `M`.  
**Outputs:** Vector of probabilities (length M).


#### 9. `expected_LTHC_t_after_tau`
**Purpose:**  
Expected diagnosis time after `tau` under the **baseline** predictive.

**Inputs:** `parameters`, `hyperparameters`, `tau`, `M`.  
**Outputs:** Vector of expected times (length M).


## Repository Layout

### `data/` - dataset generators (synthetic)
- `generate_presence_based_clusters.R`  
  **Purpose:** Presence‑driven clusters; diagnosis = onset + truncated‑normal delay (cluster‑independent).  
  **Output:** `data/generated_presence_promote_style.rds`.
- `generate_onset_based_clusters.R`  
  **Purpose:** Onset‑driven clusters with truncated‑normal delays; presence uninformative.  
  **Output:** `data/generated_onset_promote_style.rds`.
- `generate_mixed_delays.R`  
  **Purpose:** Onset‑driven clusters with **heterogeneous delay families** (gaussian, uniform, two‑component mixture); stores per‑condition delay priors (`delay_prior_df`).  
  **Output:** `data/generated_promote_style_mixed_delays.rds`.
- `generate_mixed_delays_bimodal_flag.R`  
  **Purpose:** As above, plus a per‑condition **early/late (bimodal) flag** for analysis.  
  **Output:** *(see script; same `.rds` shape with extra flags).*

> Train/test splits in downstream scripts use `set.seed(42)` and an 80/20 split (≈160k/40k for `N=200000`).

### `functions/` - baseline (no‑delay)
- `ProMOTe_VB.R` - `VB_gaussian_update`.
- `ProMOTe_Predictive.R` - `VB_gaussian_predictive_density`.
- `ProMOTe_LTCby.R` - `probability_LTHC_by_T`.
- `ProMOTe_LTCt.R` - `expected_LTHC_t_after_tau`.
- `ProMOTe_utility.R` - `expected_lst_lefttrunc`

### `functionswithdelay/` - delay‑aware
- `delay_VB.R` - `VB_gaussian_update_d`.
- `ProMOTe_Predictive_delay.R` - `VB_gaussian_predictive_density_d`.
- `ProMOTe_LTCby_delay.R` - `probability_LTHC_by_T_d`.
- `ProMOTe_LTCt_delay.R` - `expected_LTHC_t_after_tau_d`.
- `ProMOTe_utility_delay.R` - `expected_lst_lefttrunc_d`

### `scripts/` - training, prediction, analysis, plots
- `prediction.R`  
  **Purpose:** End‑to‑end baseline vs delay‑aware comparison on mixed‑delay data; cluster recovery, forward prediction, MAE.
- `prediction_mixed_priors.R`  
  **Purpose:** Predict with stored condition‑level delay priors; saves per‑patient probabilities/expectations.
- `ablation_test.R` / `ablation_plot.R`  
  **Purpose:** Prior‑strength ablation over delay variance scales; writes CSVs and plots (cluster recovery, shrinkage).
- `identifiability_check.R`  
  **Purpose:** Swap‑correlation and invariance diagnostics across different delay‑prior scales.
- `censoring_percentages.R`  
  **Purpose:** Per‑cluster left/right censoring summaries from `d`, `t`, `rho`, `tau`, `iota`.
- `bimodal_in_mixed.R`  
  **Purpose:** Bimodal‑flag analysis on mixed‑delay dataset; chi‑square early/late vs clusters; MAE on late subgroup.
- `plot_elbo.R`  
  **Purpose:** Compare ELBO trajectories for delay‑aware vs no‑delay models and plot parameter differences against elbo.
- `plot_risk_profiles.R`  
  **Purpose:** Visualise cluster‑specific trajectories from the delay‑aware posterior.
- `bubble_plot.R`  
  **Purpose:** TEST‑split bubble chart of inferred disease presence by cluster using the delay‑aware posterior.

### `src/plots/`
Saved figures from ablations and cluster trajectory visualisations (e.g., prior‑scale vs ARI/NMI, shrinkage boxplots, onset/diagnosis overlays).

## Notes on Censoring and Masks

- **Left‑censored:** `d==1` with `t < rho` (diagnosis before baseline).  
- **Right‑censored:** `d==0` with survival state indicating no event by `tau`.  
- Birth‑only and sex‑specific conditions are hard‑coded via `birth_conds`, `male_conds`, `female_conds` where provided.


## Reproducibility

- Global random seeds are set in generation and analysis scripts (commonly `set.seed(42)`).
- All splits, priors, and hyperparameters used in ablations are recorded in the saved RDS/CSV outputs.

