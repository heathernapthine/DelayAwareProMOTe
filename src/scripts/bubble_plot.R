library(ggforestplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(clue)
library(scales)     # percent_format()

set.seed(42)

# --- Load & split as you had it ------------------------------------------------
data_all <- readRDS("data/generated_promote_style_mixed_delays.rds")

n_total   <- nrow(data_all$d)
train_prop <- 0.8
train_idx <- sample(1:n_total, size = floor(train_prop * n_total))
test_idx  <- setdiff(1:n_total, train_idx)

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
    male_conds = data$male_conds,
    female_conds = data$female_conds,
    cond_list = data$cond_list,
    delay_prior_df  = if (!is.null(data$delay_prior_df))  data$delay_prior_df  else NULL,
    delay_dist_cond = if (!is.null(data$delay_dist_cond)) data$delay_dist_cond else NULL,
    delay_mu_cond   = if (!is.null(data$delay_mu_cond))   data$delay_mu_cond   else NULL,
    delay_sd_cond   = if (!is.null(data$delay_sd_cond))   data$delay_sd_cond   else NULL
  )
}

train_data <- subset_data(data_all, train_idx)
test_data  <- subset_data(data_all,  test_idx)

# --- Pull inferred prevalence --------------------------------------------------
posterior_train <- readRDS("src/elbonew/posterior_val_delay_train.rds")
prevalence <- posterior_train$posterior.parameters$expected_d   # dimension: M x K (assumed)

# Derive M, K from the objects we actually have
M <- ncol(train_data$d)
K <- ncol(prevalence)

# Condition labels (use your names if present)
conds <- if (!is.null(train_data$cond_list)) train_data$cond_list else as.character(seq_len(M))

# Baseline prevalence across all patients (to compute RR; protect zeros)
baseline <- colMeans(train_data$d)
eps <- 1e-8
baseline_adj <- pmax(baseline, eps)

# --- True prevalence by cluster (from training iota) ---------------------------
true_prevalence <- matrix(0, nrow = M, ncol = K)
for (k in 1:K) {
  idx_k <- which(train_data$iota == k)
  if (length(idx_k) > 0) {
    true_prevalence[, k] <- colMeans(train_data$d[idx_k, , drop = FALSE])
  }
}

# Align inferred cluster columns to true clusters using LSAP on prevalence profiles
cost_matrix <- matrix(0, nrow = K, ncol = K)
for (i in 1:K) {
  for (j in 1:K) {
    cost_matrix[i, j] <- sum((prevalence[, i] - true_prevalence[, j])^2)
  }
}
reorder <- solve_LSAP(cost_matrix)
true_prevalence <- true_prevalence[, reorder, drop = FALSE]

# Relative risk (optional, used only for alpha)
relative_risk      <- prevalence      / baseline_adj
true_relative_risk <- true_prevalence / baseline_adj

# --- Long data frame -----------------------------------------------------------
prevalencedf <- tibble(
  cluster          = factor(rep(seq_len(K), each = M), levels = seq_len(K)),
  condition        = factor(rep(conds, times = K), levels = conds),
  prevalence       = as.vector(prevalence),
  true_prevalence  = as.vector(true_prevalence),
  rr               = as.vector(relative_risk),
  true_rr          = as.vector(true_relative_risk)
)

# Keep the top N conditions by baseline prevalence (makes x-axis readable)
top_n <- 40
keep_conds <- names(sort(baseline, decreasing = TRUE))[seq_len(min(top_n, length(baseline)))]
prevalencedf <- prevalencedf %>%
  filter(condition %in% keep_conds) %>%
  mutate(condition = factor(condition, levels = keep_conds))  # order by overall prevalence

# --- A reusable plotting function ---------------------------------------------
dotplot <- function(dat, value_col, rr_col, title = NULL) {
  ggplot(dat, aes(x = condition, y = cluster)) +
    # alternate grey stripes by CLUSTER rows (y-axis)
    ggforestplot::geom_stripes(
      aes(y = as.numeric(cluster)),
      odd = "#00000008", even = "#00000000",
      inherit.aes = FALSE
    ) +
    # bubbles: size & fill by prevalence, subtle outline
    geom_point(
      aes(size = {{ value_col }}, fill = {{ value_col }}, alpha = {{ rr_col }}),
      shape = 21, colour = "black", stroke = 0.15
    ) +
    scale_size_area(
      name = "Prevalence",
      max_size = 6,
      breaks = c(0.01, 0.05, 0.10),
      labels = percent_format(accuracy = 1)
    ) +
    scale_fill_viridis_c(name = "Prevalence", direction = -1) +
    scale_alpha_continuous(range = c(0.35, 1), guide = "none") +
    labs(x = "Condition", y = "Cluster", title = title) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.margin = margin(10, 20, 10, 10)
    )
}

p_inferred <- dotplot(prevalencedf, value_col = prevalence,      rr_col = rr,      title = "Inferred prevalence")
p_true     <- dotplot(prevalencedf, value_col = true_prevalence, rr_col = true_rr, title = "True prevalence")

ggsave("inferred_prevalence.png", p_inferred, width = 14, height = 8, dpi = 300)
ggsave("true_prevalence.png",     p_true,     width = 14, height = 8, dpi = 300)
