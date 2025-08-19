# TEST-split bubble chart of inferred disease presence by cluster using the delay-aware predictive;
# bubble size = prevalence, opacity = log-relative-risk vs baseline, with a right-hand bar chart of
# ground-truth cluster share (test), aligned to predicted clusters, sorted decreasing.

library(ggplot2)
library(ggforestplot)
library(clue)
library(dplyr)
library(stringr)
library(scales)
library(patchwork)

set.seed(42)

# Load data & posterior
data_all <- readRDS("data/generated_promote_style_mixed_delays.rds")
posterior_train <- readRDS("src/resultsmixeddata/posterior_val_delay_train.rds")

source("src/functionswithdelay/ProMOTe_Predictive_delay.R")  # VB_gaussian_predictive_density_d
source("src/functionswithdelay/ProMOTe_LTCby_delay.R")       # probability_LTHC_by_T_d
source("src/functionswithdelay/ProMOTe_utility_delay.R")     # helpers etc.

n_total    <- nrow(data_all$d)
train_prop <- 0.8
train_idx  <- sample(seq_len(n_total), size = floor(train_prop * n_total))
test_idx   <- setdiff(seq_len(n_total), train_idx)

subset_data <- function(data, idx) {
  list(
    d = data$d[idx, , drop = FALSE],
    t = data$t[idx, , drop = FALSE],
    onset = data$onset[idx, , drop = FALSE],
    z = data$z[idx],
    rho = data$rho[idx],
    tau = data$tau[idx],
    iota = data$iota[idx],
    N = length(idx),
    M = data$M,
    K = data$K,
    cond_list = data$cond_list,
    delay_prior_df  = data$delay_prior_df,
    delay_dist_cond = data$delay_dist_cond,
    delay_mu_cond   = data$delay_mu_cond,
    delay_sd_cond   = data$delay_sd_cond
  )
}

test_data <- subset_data(data_all, test_idx)
M <- test_data$M
conds_full  <- if (!is.null(test_data$cond_list)) test_data$cond_list else as.character(seq_len(M))
cond_labels <- str_trunc(as.character(conds_full), width = 15, side = "right", ellipsis = "…")

# Build predictive (test-time) hyperparameters from trained posterior 
pp <- posterior_train$posterior.parameters
a_post     <- pp$pi_a
b_post     <- pp$pi_b
u_post     <- pp$mu_u
v_post     <- pp$v_star
alpha_post <- pp$mu_alpha
beta_post  <- pp$mu_beta
theta_post <- if (!is.null(pp$gamma_alpha)) pp$gamma_alpha / sum(pp$gamma_alpha) else rep(1 / ncol(a_post), ncol(a_post))
hyper_post <- list(theta_post, a_post, b_post, u_post, v_post, alpha_post, beta_post)
K <- length(theta_post)

# Condition-specific delay priors (mu0, sigma20) 
df <- test_data$delay_prior_df
stopifnot(!is.null(df))
df <- df[match(conds_full, df$condition), ]

is_g <- df$delay_dist == "gaussian"
is_u <- df$delay_dist == "uniform"
is_m <- df$delay_dist == "mixture2"

mu0     <- numeric(M)
sigma20 <- numeric(M)

mu0[is_g]     <- df$delay_mu[is_g]
sigma20[is_g] <- (df$delay_sd[is_g])^2

mu0[is_u]     <- (df$unif_a[is_u] + df$unif_b[is_u]) / 2
sigma20[is_u] <- (df$unif_b[is_u] - df$unif_a[is_u])^2 / 12

mix_mean <- df$mix_w1*df$mix_mu1 + (1 - df$mix_w1)*df$mix_mu2
mix_var  <- df$mix_w1*(df$mix_sd1^2 + (df$mix_mu1 - mix_mean)^2) +
            (1 - df$mix_w1)*(df$mix_sd2^2 + (df$mix_mu2 - mix_mean)^2)
mu0[is_m]     <- mix_mean[is_m]
sigma20[is_m] <- mix_var[is_m]

# Test-time responsibilities & expected presence 
N_test <- test_data$N
phi_mat  <- matrix(NA_real_, N_test, K) # responsibilities per patient
E_d_test <- matrix(NA_real_, N_test, M) # E[d_ij] over [rho_i, tau_i]

make_view <- function(d_row, t_row, rho_i, tau_i, M) {
  M_obs_idx  <- which(d_row == 1 & !is.na(t_row) & t_row >= rho_i & t_row <= tau_i)
  M_part_idx <- which(d_row == 1 & !is.na(t_row) & t_row <  rho_i)
  all_idx    <- seq_len(M)
  M_unobs_idx<- setdiff(all_idx, union(M_obs_idx, M_part_idx))
  list(
    M_obs  = M_obs_idx,
    M_part = M_part_idx,
    M_unobs= M_unobs_idx,
    d_obs  = rep.int(1L, length(M_obs_idx)),
    t_obs  = if (length(M_obs_idx)) t_row[M_obs_idx] else numeric(0),
    d_part = rep.int(1L, length(M_part_idx))
  )
}

for (i in seq_len(N_test)) {
  vw <- make_view(test_data$d[i, ], test_data$t[i, ], test_data$rho[i], test_data$tau[i], M)

  pred <- VB_gaussian_predictive_density_d(
    hyperparameters = hyper_post,
    M_obs  = vw$M_obs,  M_part = vw$M_part,  M_unobs= vw$M_unobs,
    d_obs  = vw$d_obs,  t_obs  = vw$t_obs,   d_part = vw$d_part,
    rho    = test_data$rho[i],
    tau    = test_data$tau[i],
    M      = M,
    mu0    = mu0,
    sigma20 = sigma20
  )

  # Responsibilities.
  phi_row <- as.numeric(pred$phi)
  if (!all(is.finite(phi_row)) || sum(phi_row) <= 0) {
    phi_row <- rep(1 / K, K)
  } else {
    phi_row[!is.finite(phi_row)] <- 0
    s <- sum(phi_row); if (!is.finite(s) || s <= 0) s <- 1
    phi_row <- pmax(phi_row, 0) / s
  }
  phi_mat[i, ] <- phi_row

  # Expected presence by end of window (upper bound = tau_i).
  E_d_row <- probability_LTHC_by_T_d(
    parameters = pred, hyperparameters = hyper_post,
    T = test_data$tau[i], tau = test_data$rho[i],
    M = M, mu0 = mu0, sigma20 = sigma20
  )
  if (length(E_d_row) != M || !all(is.finite(E_d_row))) {
    tmp <- rep(0, M); tmp[vw$M_obs] <- 1
    E_d_row <- tmp
  }
  E_d_test[i, ] <- pmin(pmax(E_d_row, 0), 1)
}

# Responsibilities & inferred prevalence on test.
row_sums <- rowSums(phi_mat); row_sums[!is.finite(row_sums) | row_sums <= 0] <- 1
R_test <- sweep(phi_mat, 1, row_sums, "/")

num_test   <- t(E_d_test) %*% R_test                        # M x K
denom_test <- matrix(colSums(R_test), nrow = M, ncol = K, byrow = TRUE)
prevalence_test <- num_test / denom_test                    # M x K
prevalence_test[!is.finite(prevalence_test)] <- 0

# Align predicted with true and ORDER by decreasing GT share.
pred_hard <- max.col(R_test, ties.method = "first")
z_true_f  <- factor(as.integer(test_data$z), levels = seq_len(K))
pred_f    <- factor(pred_hard,               levels = seq_len(K))
cm_mat    <- as.matrix(table(z_true_f, pred_f))            

perm_true_to_pred <- as.integer(clue::solve_LSAP(cm_mat, maximum = TRUE))
perm_pred_to_true <- match(seq_len(K), perm_true_to_pred)

gt_counts <- tabulate(as.integer(test_data$z), nbins = K)
gt_share  <- gt_counts / sum(gt_counts)

gt_share_aligned <- gt_share[perm_pred_to_true] # predicted-cluster index
ord <- order(gt_share_aligned, decreasing = TRUE)           

# Apply this order to everything we plot by cluster (cols)
prevalence_test <- prevalence_test[, ord, drop = FALSE]
mass_test       <- colSums(R_test)[ord]                      
K_plot <- ncol(prevalence_test)
cluster_labels <- as.character(seq_len(K_plot))

# Relative risk vs baseline.
baseline_test <- colMeans(test_data$d)
baseline_safe <- pmax(baseline_test, 1e-8)
relative_risk_test <- sweep(prevalence_test, 1, baseline_safe, "/")

# Map log-RR to opacity.
eps <- 1e-8
log_rr_cap <- 2
log_rr <- log(pmax(relative_risk_test, eps))
log_rr <- pmax(pmin(log_rr, log_rr_cap), -log_rr_cap)

prevalence_df <- data.frame(
  cond_id    = rep(seq_len(M), times = K_plot),
  cluster    = factor(rep(cluster_labels, each = M), levels = cluster_labels),
  condition  = rep(cond_labels, times = K_plot),
  prevalence = as.vector(pmin(pmax(prevalence_test, 0), 1)),
  log_rr     = as.vector(log_rr)
)

# Keep same half-condition selection.
set.seed(42)
cond_keep_idx     <- sort(sample(seq_len(M), size = floor(M/2)))
cond_keep_labels  <- cond_labels[cond_keep_idx]

prevalence_df <- prevalence_df %>%
  filter(cond_id %in% cond_keep_idx) %>%
  mutate(condition = factor(condition, levels = cond_keep_labels))

pal <- scales::hue_pal()(K_plot)

# Bubble plot.
alpha_min <- 0.2
alpha_max <- 1.00
log_rr_cap <- 1   


p_bubbles <- ggplot(prevalence_df, aes(x = condition, y = cluster)) +
  ggforestplot::geom_stripes(
    aes(y = as.integer(cluster)), odd = "#33333322", even = "#00000000", inherit.aes = FALSE
  ) +
  geom_point(
    aes(size = prevalence, fill = cluster, alpha = log_rr),
    shape = 21, colour = "black", stroke = 0.15
  ) +
  scale_size_continuous(
    name = "Disease Prevalence",
    range = c(0.8, 7), limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = c("0.00","0.25","0.50","0.75","1.00")
  ) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_alpha_continuous(
    name   = "Relative Risk",
    limits = c(-log_rr_cap, log_rr_cap),
    range  = c(alpha_min, alpha_max),
    breaks = log(c(0.5, 1, 2)),
    labels = c("Less than Half", "Baseline", "Double"),
    oob    = scales::squish
  ) +
  scale_y_discrete(limits = rev(levels(prevalence_df$cluster))) +
  scale_x_discrete(limits = levels(prevalence_df$condition)) +
  labs(x = "Condition", y = "Cluster", title = "Posterior Disease Presence") +
  coord_cartesian(clip = "off") +
  (theme_forest() +
     theme(
       panel.background = element_rect(fill = "white", colour = "white"),
       plot.background  = element_rect(fill = "white", colour = "white"),
       text             = element_text(size = 16),
       legend.title     = element_text(size = 14),
       legend.text      = element_text(size = 12),
       legend.position  = "left",
       axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
       axis.text.y      = element_text(size = 10),
       plot.title       = element_text(hjust = 0, face = "bold", size = 18),
       panel.grid.major.x = element_blank(),
       panel.grid.minor   = element_blank()
     )) +
  guides(
    size  = guide_legend(order = 1, override.aes = list(shape = 21, fill = "grey85")),
    alpha = guide_legend(order = 2, override.aes = list(shape = 21, fill = "grey50", size = 5))
  )

# Right-hand bars: ground-truth share (sorted decreasing) 
gt_share_plot <- gt_share_aligned[ord]   

bar_df <- data.frame(
  cluster    = factor(cluster_labels, levels = cluster_labels),
  mass_share = gt_share_plot
)

p_bar_right <- ggplot(bar_df, aes(y = cluster, x = mass_share, fill = cluster)) +
  geom_col(width = 0.9) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_discrete(limits = rev(levels(prevalence_df$cluster))) +
  scale_x_continuous(
    limits = c(0, 0.25), # cap at 50%
    labels = scales::label_percent(accuracy = 1),
    breaks = c(0, 0.25, 0.5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 0)
  )

# Combine & save.
combo <- p_bubbles + p_bar_right + patchwork::plot_layout(widths = c(8, 1.2))
ggsave("src/plots/test_inferred_prevalence.png", combo, width = 12, height = 6, dpi = 300)


