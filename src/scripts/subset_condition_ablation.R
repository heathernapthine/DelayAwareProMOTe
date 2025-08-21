# Loads synthetic multimorbidity dataset with condition-specific delay priors 
# Fits the delay-aware ProMOTe model under different prior-variance scalings
# Evaluates clustering (ACC/ARI/NMI) and forward prediction (Accuracy, AUROC, MAE) 
# restricting the predictive density calculation (VB_gaussian_predictive_density_d)
# to only a chosen subset of conditions (columns).


set.seed(42)

library(mclust); library(aricode); library(pROC)
library(dplyr);  library(tidyr);   library(stringr); library(clue)
library(truncnorm); library(MCMCpack); library(purrr)

source("src/functionswithdelay/delay_VB.R")                
source("src/functionswithdelay/ProMOTe_LTCby_delay.R")      
source("src/functionswithdelay/ProMOTe_LTCt_delay.R")       
source("src/functionswithdelay/ProMOTe_Predictive_delay.R") 
source("src/functionswithdelay/ProMOTe_utility_delay.R")    

data_all <- readRDS("data/generated_promote_style_mixed_delays.rds")
n_total <- nrow(data_all$d); train_prop <- 0.8
train_idx <- sample.int(n_total, size = floor(train_prop * n_total))
test_idx  <- setdiff(seq_len(n_total), train_idx)

subset_data <- function(data, idx){
  list(
    d = data$d[idx,, drop=FALSE], t = data$t[idx,, drop=FALSE],
    rho = data$rho[idx], tau = data$tau[idx], iota = data$iota[idx],
    N = length(idx), M = data$M, cond_list = data$cond_list,
    sex = if (!is.null(data$sex)) data$sex[idx] else NULL,
    birth_conds = data$birth_conds, male_conds = data$male_conds,
    female_conds = data$female_conds, delay_prior_df = data$delay_prior_df
  )
}
train_data <- subset_data(data_all, train_idx)
test_data  <- subset_data(data_all, test_idx)

K <- 10; epsilon <- 0.1; M <- train_data$M; cond_list <- train_data$cond_list

# Delay priors (base)
df <- train_data$delay_prior_df; df <- df[match(cond_list, df$condition), ]
delay_family <- df$delay_dist

mu0_base     <- numeric(M); sigma20_base <- numeric(M)
is_g <- delay_family == "gaussian"; is_u <- delay_family == "uniform"; is_m <- delay_family == "mixture2"

mu0_base[is_g]     <- df$delay_mu[is_g]
sigma20_base[is_g] <- (df$delay_sd[is_g])^2

mu0_base[is_u]     <- (df$unif_a[is_u] + df$unif_b[is_u]) / 2
sigma20_base[is_u] <- (df$unif_b[is_u] - df$unif_a[is_u])^2 / 12

mix_mean <- df$mix_w1*df$mix_mu1 + (1 - df$mix_w1)*df$mix_mu2
mix_var  <- df$mix_w1*(df$mix_sd1^2 + (df$mix_mu1 - mix_mean)^2) +
            (1 - df$mix_w1)*(df$mix_sd2^2 + (df$mix_mu2 - mix_mean)^2)
mu0_base[is_m]     <- mix_mean[is_m]
sigma20_base[is_m] <- mix_var[is_m]

# Weakly-informative globals (fixed across runs)
theta_h <- rep(1, K)
a_h     <- matrix(1,   M, K)
b_h     <- matrix(1,   M, K)
u_h     <- matrix(50,  M, K)
v_h     <- matrix(0.3, M, K)
alpha_h <- matrix(5,   M, K)
beta_h  <- matrix(750, M, K)
hyper_h <- list(theta_h, a_h, b_h, u_h, v_h, alpha_h, beta_h)

make_view <- function(d_row, t_row, rho_i, tau_i, M){
  M_obs  <- which(d_row==1 & !is.na(t_row) & t_row>=rho_i & t_row<=tau_i)
  M_part <- which(d_row==1 & !is.na(t_row) & t_row< rho_i)
  all    <- seq_len(M); M_unobs <- setdiff(all, union(M_obs, M_part))
  list(M_obs=M_obs, M_part=M_part, M_unobs=M_unobs,
       d_obs=rep.int(1L, length(M_obs)),
       t_obs=if(length(M_obs)) t_row[M_obs] else numeric(0),
       d_part=rep.int(1L, length(M_part)))
}
make_view_subset <- function(d_row, t_row, rho_i, tau_i, M, cols){
  v <- make_view(d_row, t_row, rho_i, tau_i, M)
  v$M_obs  <- intersect(v$M_obs,  cols)
  v$M_part <- intersect(v$M_part, cols)
  v$M_unobs <- setdiff(seq_len(M), union(v$M_obs, v$M_part))
  v$d_obs  <- rep.int(1L, length(v$M_obs))
  v$t_obs  <- if(length(v$M_obs)) t_row[v$M_obs] else numeric(0)
  v$d_part <- rep.int(1L, length(v$M_part))
  v
}
hungarian_align <- function(pred, truth, K){
  K_used <- max(K, max(pred), max(truth))
  tab <- table(factor(pred, levels=1:K_used), factor(truth, levels=1:K_used))
  map <- as.integer(clue::solve_LSAP(tab, maximum = TRUE))
  map[pred]
}
safe_auc <- function(y, p){ if(length(unique(y))>1) as.numeric(pROC::auc(pROC::roc(y,p,quiet=TRUE))) else NA_real_ }

directory <- "src/ablation_5run"
if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)

# Core runner (fit once, evaluate, return top/shrink + prediction objects)
run_fit <- function(mu0, sigma20, tag){
  save_train     <- file.path(directory, sprintf("train_%s.rds", tag))
  save_test_full <- file.path(directory, sprintf("test_full_%s.rds", tag))

  need_train <- TRUE
  if (file.exists(save_train)) {
    posterior_train <- readRDS(save_train)
    if (!is.null(posterior_train$posterior.parameters$gap_sigma2_star)) need_train <- FALSE
  }
  if (need_train){
    set.seed(42 + as.integer(runif(1,1e6,1e7)))
    N_train <- train_data$N
    init_C <- t(rdirichlet(N_train, rep(1, K))); init_C <- t(init_C)
    init_D <- matrix(runif(N_train*M,0,1), N_train, M)
    init_p <- matrix(runif(N_train*M,0,10), N_train, M)
    init_q <- matrix(runif(N_train*M,1,2),  N_train, M)
    init_r <- matrix(runif(N_train*M,0.01,0.02), N_train, M)

    posterior_train <- VB_gaussian_update_d(
      t_obs=train_data$t, d=train_data$d, rho=train_data$rho, tau=train_data$tau,
      iota=train_data$iota, hyperparameters=hyper_h,
      initial_Cstar=init_C, initial_Dstar=init_D, initial_pstar=init_p,
      initial_qstar=init_q, initial_rstar=init_r,
      N=train_data$N, M=M, K=K, epsilon=epsilon,
      mu0=mu0, sigma20=pmax(sigma20,1e-8),
      sex=train_data$sex, birth_conds=train_data$birth_conds,
      male_conds=train_data$male_conds, female_conds=train_data$female_conds,
      cond_list=cond_list
    )
    saveRDS(posterior_train, save_train)
  }

  pp <- posterior_train$posterior.parameters
  theta_post <- pp$gamma_alpha / sum(pp$gamma_alpha)
  a_post <- pp$pi_a; b_post <- pp$pi_b
  u_post <- pp$mu_u; v_post <- pp$v_star
  alpha_post <- pp$mu_alpha; beta_post <- pp$mu_beta
  hyper_post <- list(theta_post, a_post, b_post, u_post, v_post, alpha_post, beta_post)

  # Shrinkage (train)
  gap_mu_tr  <- pp$gap_mu_star
  gap_var_tr <- pp$gap_sigma2_star
  ed_tr      <- pp$expected_d
  denom_w <- pmax(colSums(ed_tr, na.rm=TRUE), 1e-8)

  SR   <- colMeans(gap_var_tr, na.rm=TRUE) / sigma20
  PP   <- (colMeans(gap_mu_tr, na.rm=TRUE) - mu0) / sqrt(sigma20)
  SR_w <- (colSums(gap_var_tr*ed_tr, na.rm=TRUE) / denom_w) / sigma20
  PP_w <- ((colSums(gap_mu_tr*ed_tr, na.rm=TRUE)/denom_w) - mu0) / sqrt(sigma20)
  shrink_df <- data.frame(tag=tag, condition=cond_list,
                          family=delay_family, SR=SR, SR_w=SR_w, PP=PP, PP_w=PP_w)

  # Cluster recovery (test)
  N_test <- test_data$N
  phi_mat <- matrix(NA_real_, N_test, K)
  for (i in 1:N_test){
    v <- make_view(test_data$d[i,], test_data$t[i,], test_data$rho[i], test_data$tau[i], M)
    pred <- VB_gaussian_predictive_density_d(
      hyperparameters=hyper_post,
      M_obs=v$M_obs, M_part=v$M_part, M_unobs=v$M_unobs,
      d_obs=v$d_obs, t_obs=v$t_obs, d_part=v$d_part,
      rho=test_data$rho[i], tau=test_data$tau[i], M=M,
      mu0=mu0, sigma20=sigma20
    )
    phi_mat[i,] <- pred$phi
  }
  truth <- data_all$z[test_idx]
  aligned <- hungarian_align(max.col(phi_mat, ties.method="first"), truth, K)
  acc <- mean(aligned==truth); ari <- mclust::adjustedRandIndex(aligned, truth)
  nmi <- aricode::NMI(aligned, truth)

  # Forward prediction (after random cut age)
  set.seed(4242)
  cut_ages <- runif(N_test, 50, 90)
  cut_mat  <- matrix(cut_ages, nrow=N_test, ncol=M)
  phi_pre <- matrix(NA_real_, N_test, K)
  pred_list_pre <- vector("list", N_test)
  for (i in 1:N_test){
    v <- make_view(test_data$d[i,], test_data$t[i,], test_data$rho[i], cut_ages[i], M)
    pred <- VB_gaussian_predictive_density_d(
      hyperparameters=hyper_post,
      M_obs=v$M_obs, M_part=v$M_part, M_unobs=v$M_unobs,
      d_obs=v$d_obs, t_obs=v$t_obs, d_part=v$d_part,
      rho=test_data$rho[i], tau=cut_ages[i], M=M,
      mu0=mu0, sigma20=sigma20
    )
    phi_pre[i,] <- pred$phi; pred_list_pre[[i]] <- pred
  }
  P_after <- matrix(0, N_test, M); E_t_after <- matrix(NA_real_, N_test, M)
  for (i in 1:N_test){
    P_after[i,] <- probability_LTHC_by_T_d(pred_list_pre[[i]], hyper_post,
                                           T=test_data$tau[i], tau=cut_ages[i], M=M, mu0=mu0, sigma20=sigma20)
    E_t_after[i,] <- expected_LTHC_t_after_tau_d(pred_list_pre[[i]], hyper_post,
                                                tau=cut_ages[i], M=M, mu0=mu0, sigma20=sigma20)
  }

  d_true <- test_data$d; t_true <- test_data$t; tau_true <- test_data$tau
  is_pos <- (d_true==1) & (t_true>cut_mat) & (t_true<=tau_true)
  is_neg <- (d_true==0) | ((d_true==1)&(t_true<=cut_mat))
  is_unk <- (d_true==1) & (t_true>tau_true)
  eval_mask <- (is_pos | is_neg) & !is_unk

  y_true <- as.integer(is_pos[eval_mask])
  y_prob <- as.numeric(P_after[eval_mask])
  y_pred <- as.integer(y_prob >= 0.5)
  acc_after <- mean(y_pred==y_true)
  auroc_after <- safe_auc(y_true, y_prob)

  mae_mask <- is_pos & !is.na(E_t_after)
  mae_vals <- abs(E_t_after[mae_mask] - t_true[mae_mask])
  mae_after <- if (length(mae_vals)>0) mean(mae_vals) else NA_real_

  top <- data.frame(tag=tag, cluster_acc=acc, cluster_ari=ari, cluster_nmi=nmi,
                    after_acc=acc_after, after_auroc=auroc_after, after_mae=mae_after)

  saveRDS(list(top=top, shrink=shrink_df,
               hyper_post=hyper_post, mu0=mu0, sigma20=sigma20),
          save_test_full)
  list(top=top, shrink=shrink_df, hyper_post=hyper_post, mu0=mu0, sigma20=sigma20)
}

# RUN 1: Baseline (s=1)
mu0 <- mu0_base
sigma20 <- sigma20_base
res_base <- run_fit(mu0, sigma20, tag="baseline_s1")

# choose family & condition to focus on (most "needy" at baseline: highest median SR_w / highest SR_w)
sr_base <- res_base$shrink
fam_order <- sr_base %>% group_by(family) %>% summarise(median_SRw=median(SR_w, na.rm=TRUE)) %>%
  arrange(desc(median_SRw))
family_focus <- fam_order$family[1]

cond_order <- sr_base %>% arrange(desc(SR_w))
cond_focus <- cond_order$condition[1]
cond_focus_idx <- match(cond_focus, cond_list)

# scale sigma2 by family or condition
scale_sigma20_by_family <- function(sig2, fam, fam_name, scale){
  s2 <- sig2; s2[fam==fam_name] <- pmax(s2[fam==fam_name]*scale, 1e-8); s2
}
scale_sigma20_by_condition <- function(sig2, idx, scale){
  s2 <- sig2; s2[idx] <- pmax(s2[idx]*scale, 1e-8); s2
}

# RUN 2: Family tightened (0.25x)  
sigma20_f_tight <- scale_sigma20_by_family(sigma20_base, delay_family, family_focus, 0.25)
res_fam_tight <- run_fit(mu0_base, sigma20_f_tight, tag=sprintf("fam_%s_tight_s0p25", family_focus))

# RUN 3: Family loosened (4x)  
sigma20_f_loose <- scale_sigma20_by_family(sigma20_base, delay_family, family_focus, 4.0)
res_fam_loose <- run_fit(mu0_base, sigma20_f_loose, tag=sprintf("fam_%s_loose_s4", family_focus))

# RUN 4: Condition tightened (0.25x)  
sigma20_c_tight <- scale_sigma20_by_condition(sigma20_base, cond_focus_idx, 0.25)
res_cond_tight <- run_fit(mu0_base, sigma20_c_tight, tag=sprintf("cond_%03d_tight_s0p25", cond_focus_idx))

# RUN 5: Condition loosened (4x)  
sigma20_c_loose <- scale_sigma20_by_condition(sigma20_base, cond_focus_idx, 4.0)
res_cond_loose <- run_fit(mu0_base, sigma20_c_loose, tag=sprintf("cond_%03d_loose_s4", cond_focus_idx))

# PRIORITISATION (no fits): ARI with only chosen family/condition as evidence  
cluster_ari_with_cols <- function(cols_use, hyper_post, mu0, sigma20){
  N_test <- test_data$N
  phi_mat <- matrix(NA_real_, N_test, K)
  for (i in 1:N_test){
    v <- make_view_subset(test_data$d[i,], test_data$t[i,], test_data$rho[i], test_data$tau[i], M, cols_use)
    pred <- VB_gaussian_predictive_density_d(
      hyperparameters=hyper_post,
      M_obs=v$M_obs, M_part=v$M_part, M_unobs=v$M_unobs,
      d_obs=v$d_obs, t_obs=v$t_obs, d_part=v$d_part,
      rho=test_data$rho[i], tau=test_data$tau[i], M=M,
      mu0=mu0, sigma20=sigma20
    )
    phi_mat[i,] <- pred$phi
  }
  truth <- data_all$z[test_idx]
  aligned <- hungarian_align(max.col(phi_mat, ties.method="first"), truth, K)
  mclust::adjustedRandIndex(aligned, truth)
}

cols_family <- which(delay_family == family_focus)
ari_only_family <- cluster_ari_with_cols(cols_family, res_base$hyper_post, res_base$mu0, res_base$sigma20)

cond_idx <- cond_focus_idx
ari_only_condition <- cluster_ari_with_cols(cond_idx, res_base$hyper_post, res_base$mu0, res_base$sigma20)

w_m <- 1 / (1 + 10 * res_base$sigma20)  # implicit weights in current model

prioritisation_df <- data.frame(
  family_focus = family_focus,
  cond_focus   = cond_focus,
  cond_focus_idx = cond_focus_idx,
  ari_only_family = ari_only_family,
  ari_only_condition = ari_only_condition,
  w_cond_focus     = w_m[cond_focus_idx]
)

tops <- bind_rows(res_base$top, res_fam_tight$top, res_fam_loose$top,
                  res_cond_tight$top, res_cond_loose$top)
shrinks <- bind_rows(res_base$shrink, res_fam_tight$shrink, res_fam_loose$shrink,
                     res_cond_tight$shrink, res_cond_loose$shrink) %>%
           mutate(run_tag = tag)

write.csv(tops,    file.path(directory, "tops_5runs.csv"),    row.names=FALSE)
write.csv(shrinks, file.path(directory, "shrink_5runs.csv"),  row.names=FALSE)
write.csv(prioritisation_df, file.path(directory, "prioritisation_summary.csv"), row.names=FALSE)

cat("\n=== 5-RUN SUMMARY ===\n")
print(tops %>% mutate(across(where(is.numeric), ~round(.,4))))
cat("\nFamily focus (highest median SR_w): ", family_focus, "\n", sep="")
cat("Condition focus (highest SR_w): ", cond_focus, " [idx=", cond_focus_idx, "]\n", sep="")
cat("ARI with ONLY that family as evidence (baseline): ", round(ari_only_family,4), "\n", sep="")
cat("ARI with ONLY that condition as evidence (baseline): ", round(ari_only_condition,4), "\n", sep="")
