library(ggforestplot)

set.seed(42)

library(mclust)     
library(aricode)    
library(pROC)       
library(ggplot2)    
library(truncnorm)
library(dplyr)
library(tidyr)
library(clue)        
library(MCMCpack) 


data_all <- readRDS("data/generated_promote_style_mixed_delays.rds")

# Create an 80/20 train test split.
n_total   <- nrow(data_all$d)
train_prop <- 0.8

train_idx <- sample(1:n_total, size = floor(train_prop * n_total))
test_idx  <- setdiff(1:n_total, train_idx)

cat("Total patients:", n_total, "\n")
cat("Training patients:", length(train_idx), "\n")
cat("Test patients:", length(test_idx), "\n\n")

# Define a helper to slice the dataset by index.
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

    # Delay priors from data file are optional and can be added here.
    delay_prior_df  = if (!is.null(data$delay_prior_df))  data$delay_prior_df  else NULL,
    delay_dist_cond = if (!is.null(data$delay_dist_cond)) data$delay_dist_cond else NULL,
    delay_mu_cond   = if (!is.null(data$delay_mu_cond))   data$delay_mu_cond   else NULL,
    delay_sd_cond   = if (!is.null(data$delay_sd_cond))   data$delay_sd_cond   else NULL
  )
}

# Build train and test splits.
train_data <- subset_data(data_all, train_idx)
test_data  <- subset_data(data_all,  test_idx)

N <- 200000
M <- 75
K <- 10

posterior_train <- readRDS("src/elbonew/posterior_val_delay_train.rds")
prevalence <- posterior_train$posterior.parameters$expected_d
conds <- as.character(1:M)
baseline <- colMeans(train_data$d) 
# true prevalence: proportion of patients in each cluster who have each condition
true_prevalence <- matrix(0, nrow = M, ncol = K)

for (k in 1:K) {
  idx_k <- which(train_data$iota == k)
  if (length(idx_k) > 0) {
    true_prevalence[, k] <- colMeans(train_data$d[idx_k, , drop = FALSE])
  }
}

cost_matrix <- matrix(0, nrow = K, ncol = K)
for (i in 1:K) {
  for (j in 1:K) {
    cost_matrix[i, j] <- sum((prevalence[, i] - true_prevalence[, j])^2)
  }
}

reorder <- solve_LSAP(cost_matrix)
true_prevalence <- true_prevalence[,reorder]

relative_risk <- prevalence/baseline
true_relative_risk <- true_prevalence/baseline

prevalencedf <- data.frame(cluster = as.factor(rep(1:K, each = M)), condition = rep(conds, K), prevalence = c(prevalence), relative_risk = c(relative_risk), true_prevalence = c(true_prevalence), true_relative_risk = c(true_relative_risk))
prevalencedf$relative_risk_discretised <- as.factor(ifelse(prevalencedf$relative_risk>2, 3, ifelse(prevalencedf$relative_risk<0.5, 1, 2)))
prevalencedf$true_relative_risk_discretised <- as.factor(ifelse(prevalencedf$true_relative_risk>2, 3, ifelse(prevalencedf$true_relative_risk<0.5, 1, 2)))

prevalencedf <- prevalencedf %>% filter(condition %in% as.character(1:10))


# Create first plot and assign to variable
plot1 <- ggplot(prevalencedf, aes(x = cluster, y = condition)) + 
  geom_point(aes(size = prevalence, col = cluster, alpha = relative_risk_discretised)) + 
  scale_y_discrete(limits = conds) + 
  scale_color_discrete(guide = "none") + 
  scale_size_continuous(guide = "none") + 
  scale_alpha_discrete(guide = "none") + 
  geom_point(aes(size = prevalence, col = cluster), shape = 1) + 
  ylab("Condition") + 
  xlab("Cluster") + 
  theme_forest() %+replace% theme(panel.background = element_rect("white", "white"), 
                                  plot.background = element_rect(fill = "white", color = "white"), 
                                  text = element_text(size = 24)) + 
  geom_stripes(aes(y=unclass(condition)), odd = "#33333333", even = "#00000000", inherit.aes = FALSE)

# Create second plot and assign to variable
plot2 <- ggplot(prevalencedf, aes(x = cluster, y = condition)) + 
  geom_point(aes(size = true_prevalence, col = cluster, alpha = true_relative_risk_discretised)) + 
  scale_y_discrete(limits = conds) + 
  scale_color_discrete(guide = "none") + 
  scale_size_continuous(guide = "none") + 
  scale_alpha_discrete(guide = "none") + 
  geom_point(aes(size = true_prevalence, col = cluster), shape = 1) + 
  ylab("Condition") + 
  xlab("Cluster") + 
  theme_forest() %+replace% theme(panel.background = element_rect("white", "white"), 
                                  plot.background = element_rect(fill = "white", color = "white"), 
                                  text = element_text(size = 24)) + 
  geom_stripes(aes(y=unclass(condition)), odd = "#33333333", even = "#00000000", inherit.aes = FALSE)

# Save plots as PNG files
ggsave("prevalence_plot.png", plot1, width = 12, height = 8, dpi = 300)
ggsave("true_prevalence_plot.png", plot2, width = 12, height = 8, dpi = 300)

cat("Plots saved as:\n")
cat("- prevalence_plot.png\n")
cat("- true_prevalence_plot.png\n")