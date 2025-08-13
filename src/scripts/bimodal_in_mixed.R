# ============================================
# Mixed-delays evaluation with chi-square test
# Recomputes predictions from saved posteriors
# ============================================
set.seed(42)

suppressPackageStartupMessages({
  library(mclust)
  library(aricode)
  library(pROC)
  library(ggplot2)
  library(truncnorm)
  library(dplyr)
  library(tidyr)
  library(clue)
  library(MCMCpack)
})

# ---- Source model code ----
# Delay-aware predictive utilities
source("src/functionswithdelay/delay_VB.R")                # defines VB structures (training not used here)
source("src/functionswithdelay/ProMOTe_Predictive_delay.R")
source("src/functionswithdelay/ProMOTe_LTCby_delay.R")
source("src/functionswithdelay/ProMOTe_LTCt_delay.R")
source("src/functionswithdelay/ProMOTe_utility_delay.R")

# Source the original no-delay implementation.
source("src/functions/ProMOTe_VB.R")
source("src/functions/ProMOTe_Predictive.R")
source("src/functions/ProMOTe_LTCby.R")
source("src/functions/ProMOTe_LTCt.R")
source("src/functions/ProMOTe_utility.R")

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ---- Load data & trained posteriors ----
data_all        <- readRDS("data/generated_promote_style_mixed_delays.rds")
posterior_delay <- readRDS("src/resultsmixeddatatwelth/posterior_val_delay_train.rds")
posterior_base  <- readRDS("src/resultsmixeddatatwelth/posterior_val_no_delay_train.rds")

# ---- Train/test split (reproduced) ----
n_total   <- as.integer(nrow(data_all$d))   # ensure integer
train_prop <- 0.8
train_idx <- sample(seq_len(n_total), size = floor(train_prop * n_total))
test_idx  <- setdiff(seq_len(n_total), train_idx)

cat("Total patients:", n_total, "\n")
cat("Training patients:", length(train_idx), "\n")
cat("Test patients:", length(test_idx), "\n\n")

subset_data <- function(data, idx) {
  list(
    d = data$d[idx, , drop = FALSE],
    t = data$t[idx, , drop = FALSE],
    rho = data$rho[idx],
    tau = data$tau[idx],
    iota = data$iota[idx],
    onset = data$onset[idx, , drop = FALSE],
    N = length(idx),
    M = data$M,
    sex = if (!is.null(data$sex)) data$sex[idx] else NULL,
    birth_conds = data$birth_conds,
    male_conds  = data$male_conds,
    female_conds= data$female_conds,
    cond_list   = data$cond_list,
    delay_prior_df  = if (!is.null(data$delay_prior_df))  data$delay_prior_df  else NULL,
    delay_dist_cond = if (!is.null(data$delay_dist_cond)) data$delay_dist_cond else NULL,
    delay_mu_cond   = if (!is.null(data$delay_mu_cond))   data$delay_mu_cond   else NULL,
    delay_sd_cond   = if (!is.null(data$delay_sd_cond))   data$delay_sd_cond   else NULL
  )
}

train_data <- subset_data(data_all, train_idx)
test_data  <- subset_data(data_all,  test_idx)

M <- test_data$M
cond_list <- test_data$cond_list

# ---- Determine K from delay-aware posterior ----
K <- ncol(posterior_delay$posterior.parameters$pi_a)

# ---- Build delay prior (mu0, sigma20) per condition ----
df <- test_data$delay_prior_df
stopifnot(!is.null(df))
df <- df[match(cond_list, df$condition), ]

is_g <- df$delay_dist == "gaussian"
is_u <- df$delay_dist == "uniform"
is_m <- df$delay_dist == "mixture2"

# Use the flag provided by the generator (per condition, length M)
is_bimodal <- as.logical(data_all$is_bimodal)

mu0     <- numeric(M)
sigma20 <- numeric(M)

# Gaussian
mu0[is_g]     <- df$delay_mu[is_g]
sigma20[is_g] <- (df$delay_sd[is_g])^2

# Uniform -> mean/variance
mu0[is_u]     <- (df$unif_a[is_u] + df$unif_b[is_u]) / 2
sigma20[is_u] <- (df$unif_b[is_u] - df$unif_a[is_u])^2 / 12

# Mixture(2) -> collapsed mean/variance
mix_mean <- df$mix_w1*df$mix_mu1 + (1 - df$mix_w1)*df$mix_mu2
mix_var  <- df$mix_w1*(df$mix_sd1^2 + (df$mix_mu1 - mix_mean)^2) +
            (1 - df$mix_w1)*(df$mix_sd2^2 + (df$mix_mu2 - mix_mean)^2)
mu0[is_m]     <- mix_mean[is_m]
sigma20[is_m] <- mix_var[is_m]

# ---- Assemble hyperparameters from saved posteriors ----
# Delay-aware posterior -> predictive hyperparameters
ppd <- posterior_delay$posterior.parameters
theta_post_d <- if (!is.null(ppd$gamma_alpha)) ppd$gamma_alpha / sum(ppd$gamma_alpha) else rep(1 / K, K)
a_post_d     <- ppd$pi_a
b_post_d     <- ppd$pi_b
u_post_d     <- ppd$mu_u
v_post_d     <- ppd$v_star
alpha_post_d <- ppd$mu_alpha
beta_post_d  <- ppd$mu_beta
hyper_post_delay <- list(theta_post_d, a_post_d, b_post_d, u_post_d, v_post_d, alpha_post_d, beta_post_d)

# Baseline posterior -> predictive hyperparameters
ppb <- posterior_base$posterior.parameters
theta_post_b <- if (!is.null(ppb$theta_star)) ppb$theta_star / sum(ppb$theta_star) else rep(1 / K, K)
a_post_b     <- ppb$a_star
b_post_b     <- ppb$b_star
u_post_b     <- ppb$u_star
v_post_b     <- ppb$v_star
alpha_post_b <- ppb$alpha_star
beta_post_b  <- ppb$beta_star
hyper_post_base <- list(theta_post_b, a_post_b, b_post_b, u_post_b, v_post_b, alpha_post_b, beta_post_b)

# ---- Student-t helpers expected by predictive code ----
plst <- function(x, df, mu, sigma) {
  sigma <- pmax(sigma, 1e-12)
  pt((x - mu) / sigma, df = df)
}
dlst <- function(x, df, mu, sigma) {
  sigma <- pmax(sigma, 1e-12)
  dt((x - mu) / sigma, df = df) / sigma
}

# ---- Helper: build per-patient view for a time window [rho_i, tau_i] ----
make_view <- function(d_row, t_row, rho_i, tau_i, M) {
  M_obs_idx  <- which(d_row == 1 & !is.na(t_row) & t_row >= rho_i & t_row <= tau_i)  # fully observed
  M_part_idx <- which(d_row == 1 & !is.na(t_row) & t_row <  rho_i)                   # left-censored
  all_idx    <- seq_len(M)
  M_unobs_idx <- setdiff(all_idx, union(M_obs_idx, M_part_idx))
  list(
    M_obs  = M_obs_idx,
    M_part = M_part_idx,
    M_unobs= M_unobs_idx,
    d_obs  = rep.int(1L, length(M_obs_idx)),
    t_obs  = if (length(M_obs_idx)) t_row[M_obs_idx] else numeric(0),
    d_part = rep.int(1L, length(M_part_idx))
  )
}

# ---- Compute raw cluster predictions (un-aligned) on full window ----
N_test <- test_data$N
phi_mat_delay <- matrix(NA_real_, N_test, K)
phi_mat_base  <- matrix(NA_real_, N_test, K)

cat("Predicting raw clusters on full window ...\n")
for (i in seq_len(N_test)) {
  vview <- make_view(test_data$d[i, ], test_data$t[i, ], test_data$rho[i], test_data$tau[i], M)

  # Delay-aware predictive (needs mu0, sigma20)
  pred_d <- VB_gaussian_predictive_density_d(
    hyperparameters = hyper_post_delay,
    M_obs = vview$M_obs, M_part = vview$M_part, M_unobs = vview$M_unobs,
    d_obs = vview$d_obs, t_obs = vview$t_obs, d_part = vview$d_part,
    rho = test_data$rho[i], tau = test_data$tau[i], M = M,
    mu0 = mu0, sigma20 = sigma20
  )
  # If your codebase doesn't have *_density_d(), use VB_gaussian_predictive_density(...) with mu0/sigma20.
  phi_mat_delay[i, ] <- pred_d$phi

  # No-delay predictive
  pred_b <- VB_gaussian_predictive_density(
    hyperparameters = hyper_post_base,
    M_obs = vview$M_obs, M_part = vview$M_part, M_unobs = vview$M_unobs,
    d_obs = vview$d_obs, t_obs = vview$t_obs, d_part = vview$d_part,
    rho = test_data$rho[i], tau = test_data$tau[i], M = M
  )
  phi_mat_base[i, ] <- pred_b$phi
}

raw_pred_delay <- max.col(phi_mat_delay, ties.method = "first")
raw_pred_base  <- max.col(phi_mat_base,  ties.method = "first")

# ---- After-cut expected diagnosis ages (recomputed) ----
set.seed(42)
cut_ages <- runif(N_test, min = 50, max = 90)
cut_mat  <- matrix(cut_ages, nrow = N_test, ncol = M, byrow = FALSE)

E_t_after_delay <- matrix(NA_real_, N_test, M)
E_t_after_base  <- matrix(NA_real_, N_test, M)

cat("Computing after-cut expected diagnosis times ...\n")
for (i in seq_len(N_test)) {
  cut_i <- cut_ages[i]
  vview <- make_view(test_data$d[i, ], test_data$t[i, ], test_data$rho[i], cut_i, M)

  pred_d <- VB_gaussian_predictive_density_d(
    hyperparameters = hyper_post_delay,
    M_obs = vview$M_obs, M_part = vview$M_part, M_unobs = vview$M_unobs,
    d_obs = vview$d_obs, t_obs = vview$t_obs, d_part = vview$d_part,
    rho = test_data$rho[i], tau = cut_i, M = M,
    mu0 = mu0, sigma20 = sigma20
  )
  pred_b <- VB_gaussian_predictive_density(
    hyperparameters = hyper_post_base,
    M_obs = vview$M_obs, M_part = vview$M_part, M_unobs = vview$M_unobs,
    d_obs = vview$d_obs, t_obs = vview$t_obs, d_part = vview$d_part,
    rho = test_data$rho[i], tau = cut_i, M = M
  )

  # Expected time after cut (per condition)
  E_t_after_delay[i, ] <- expected_LTHC_t_after_tau_d(
    parameters = pred_d, hyperparameters = hyper_post_delay, tau = cut_i, M = M,
    mu0 = mu0, sigma20 = sigma20
  )
  # If *_after_tau_d() doesn't exist, use expected_LTHC_t_after_tau(...) with mu0/sigma20.
  E_t_after_base[i, ] <- expected_LTHC_t_after_tau(
    parameters = pred_b, hyperparameters = hyper_post_base, tau = cut_i, M = M
  )
}

# ---- Ground truth + evaluation masks ----
d_true   <- test_data$d
t_true   <- test_data$t
tau_true <- test_data$tau

is_pos <- (d_true == 1) & (t_true > cut_mat) & (t_true <= tau_true)   # observed after-cut events
is_neg <- (d_true == 0) | ((d_true == 1) & (t_true <= cut_mat))       # not after cut
is_unk <- (d_true == 1) & (t_true > tau_true)                         # right censored after cut
eval_mask <- (is_pos | is_neg) & !is_unk

# ---- Early/Late labels (from generator) for test split ----
delay_component <- as.matrix(data_all$delay_component[test_idx, , drop = FALSE])
stopifnot(nrow(delay_component) == N_test, ncol(delay_component) == M)

# ---- MAE for late subgroup (observed after-cut events) ----
mae_mask_mat <- is_pos & !is.na(E_t_after_delay) & !is.na(E_t_after_base)
err_delay <- abs(E_t_after_delay - t_true)
err_base  <- abs(E_t_after_base  - t_true)

is_late <- (delay_component == "late")

get_mae_by_condition_late <- function(is_late, mae_mask_mat, err_mat) {
  sapply(seq_len(M), function(m) {
    idx <- which(is_late[, m] & mae_mask_mat[, m])
    if (length(idx) > 0) {
      mean(err_mat[cbind(idx, rep(m, length(idx)))])
    } else NA_real_
  })
}

mae_delay_late <- get_mae_by_condition_late(is_late, mae_mask_mat, err_delay)
mae_base_late  <- get_mae_by_condition_late(is_late, mae_mask_mat, err_base)
mae_diff_late  <- mae_delay_late - mae_base_late

# ---- Chi-square: Early vs Late vs predicted clusters (per condition) ----
chi_pvals <- function(group_matrix, pred_clusters, is_bimodal, cond_names) {
  out <- tibble(condition = cond_names, p_value = as.numeric(NA))
  for (m in seq_along(cond_names)) {
    if (!is_bimodal[m]) next
    grp <- group_matrix[, m]
    valid <- (grp == "early" | grp == "late") & !is.na(grp)
    if (!any(valid)) next
    tab <- table(grp[valid], pred_clusters[valid])
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    if (ncol(tab) < 2 || sum(valid) < 10) next
    out$p_value[m] <- suppressWarnings(tryCatch(chisq.test(tab)$p.value, error = function(e) NA_real_))
  }
  out
}

chi_delay <- chi_pvals(delay_component, raw_pred_delay, is_bimodal, cond_list)
chi_base  <- chi_pvals(delay_component, raw_pred_base,  is_bimodal, cond_list)

# ---- Summary table ----
summary_stats <- tibble(
  condition       = cond_list,
  is_bimodal      = is_bimodal,
  mae_delay_late  = mae_delay_late,
  mae_base_late   = mae_base_late,
  mae_diff_late   = mae_delay_late - mae_base_late,
  p_base          = chi_base$p_value,
  p_delay         = chi_delay$p_value
) %>%
  mutate(
    logp_base  = ifelse(is.na(p_base),  NA_real_, -log10(pmax(p_base,  1e-300))),
    logp_delay = ifelse(is.na(p_delay), NA_real_, -log10(pmax(p_delay, 1e-300)))
  )

write.csv(summary_stats, file = "results/summary_stats_mixed_bimodal.csv", row.names = FALSE)

# ---- Plots ----
light_pink <- "#fbaee1"
light_blue <- "#A6CEE3"
purple     <- "#bd5bd3"

# 1) MAE difference for late patients (bimodal conditions only)
p1 <- summary_stats %>%
  filter(is_bimodal, !is.na(mae_diff_late)) %>%
  ggplot(aes(x = reorder(condition, mae_diff_late), y = mae_diff_late)) +
  geom_bar(stat = "identity", fill = light_blue) +
  coord_flip() +
  labs(
    title = "MAE(Delay-Aware) - MAE(Baseline) for Late Patients",
    x = "Condition", y = "MAE Difference (years)"
  ) +
  theme_classic(base_size = 12)

ggsave("results/mae_diff_late_barplot_mixed.png", p1, width = 9, height = 7, dpi = 300)

# 2) Chi-square -log10(p) for early/late vs cluster separation
# 2) Chi-square -log10(p) for early/late vs cluster separation
df_long_chi <- summary_stats %>%
  dplyr::filter(is_bimodal, !is.na(logp_base) | !is.na(logp_delay)) %>%
  dplyr::select(condition, logp_base, logp_delay) %>%
  tidyr::pivot_longer(cols = c("logp_base", "logp_delay"),
                      names_to = "model", values_to = "logp") %>%
  dplyr::mutate(model = dplyr::recode(model,
                                      logp_base = "Baseline",
                                      logp_delay = "Delay-Aware"))

p2 <- ggplot(df_long_chi, aes(x = reorder(condition, logp, na.rm = TRUE), y = logp, fill = model)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  annotate("text", x = 1, y = -log10(0.05) + 0.1, label = "p = 0.05", hjust = 0, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("Baseline" = light_pink, "Delay-Aware" = light_blue)) +
  labs(x = "Condition", y = "-log10(p-value)") +
  theme_classic(base_size = 12)

ggsave("results/chi2_logp_barplot_mixed.png", p2, width = 9, height = 7, dpi = 300)

# ---- Save generated predictions for reuse ----
saveRDS(raw_pred_delay,  "results/raw_pred_delay.rds")
saveRDS(raw_pred_base,   "results/raw_pred_base.rds")
saveRDS(E_t_after_delay, "results/E_t_after_delay.rds")
saveRDS(E_t_after_base,  "results/E_t_after_base.rds")
saveRDS(cut_ages,        "results/cut_ages.rds")

cat("Saved:\n",
    "  results/summary_stats_mixed_bimodal.csv\n",
    "  results/mae_diff_late_barplot_mixed.png\n",
    "  results/chi2_logp_barplot_mixed.png\n",
    "  results/raw_pred_delay.rds\n",
    "  results/raw_pred_base.rds\n",
    "  results/E_t_after_delay.rds\n",
    "  results/E_t_after_base.rds\n",
    "  results/cut_ages.rds\n")
