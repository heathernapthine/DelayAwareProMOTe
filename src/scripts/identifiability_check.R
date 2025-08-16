identifiability_swap_grid <- function(
  data_path  = "data/generated_promote_style_mixed_delays.rds",
  directory  = "src/ablationresultsnew",
  seed       = 42,
  train_prop = 0.8,
  save_csv   = TRUE
) {

  # Prior scales / tags (fixed)
  scales <- c(0.25, 0.5, 1, 2, 4)
  tags   <- sprintf("vscale_%s", gsub("\\.", "p", as.character(scales)))
  scale_from_tag <- function(tag) as.numeric(gsub("vscale_", "", gsub("p", ".", tag)))

  
  # Recreate the exact train split used in ablation study
  data_all <- readRDS(data_path)
  n_total  <- nrow(data_all$d)
  set.seed(seed)
  train_idx <- sample.int(n_total, size = floor(train_prop * n_total))

  d    <- data_all$d[train_idx, , drop = FALSE]
  tmat <- data_all$t[train_idx, , drop = FALSE]
  rho  <- data_all$rho[train_idx]
  tau  <- data_all$tau[train_idx]

  N <- nrow(d); M <- ncol(d)
  rho_mat <- matrix(rho, nrow = N, ncol = M, byrow = FALSE)
  tau_mat <- matrix(tau, nrow = N, ncol = M, byrow = FALSE)

  # Observed (non-censored) entries
  obs_mask <- (d == 1) & !is.na(tmat) & (tmat >= rho_mat) & (tmat <= tau_mat)
  if (!any(obs_mask)) stop("No observed (non-censored) diagnosis times under this split.")

  t_obs_vec <- as.numeric(tmat[obs_mask])
  vdt <- var(t_obs_vec, na.rm = TRUE)
  if (!is.finite(vdt) || vdt <= 0) stop("Degenerate diagnosis-time variance in observed set.")

  # Truncated-normal mean (lower bound 0) for delay
  egap_from_params <- function(mu, sigma2) {
    sigma2 <- pmax(sigma2, 1e-12)
    alpha  <- (-mu) / sqrt(sigma2)
    mu + sqrt(sigma2) * (dnorm(alpha) / pmax(1 - pnorm(alpha), 1e-12))
  }

  # Load each run and extract onset / E[gap] / sum on observed entries
  onset_list <- list(); egap_list <- list(); sum_list <- list()
  for (tg in tags) {
    pth <- file.path(directory, sprintf("train_%s.rds", tg))
    if (!file.exists(pth)) stop("Missing file: ", pth)
    pp <- readRDS(pth)$posterior.parameters

    t_onset <- pp$expected_t; t_onset[!obs_mask] <- NA_real_
    E_gap   <- egap_from_params(pp$gap_mu_star, pp$gap_sigma2_star); E_gap[!obs_mask] <- NA_real_
    S_sum   <- t_onset + E_gap

    onset_list[[tg]] <- as.numeric(t_onset[obs_mask])
    egap_list[[tg]]  <- as.numeric(E_gap[obs_mask])
    sum_list[[tg]]   <- as.numeric(S_sum[obs_mask])
  }

  # Matrices [n_obs x n_settings]
  onset_mat <- do.call(cbind, onset_list)
  egap_mat  <- do.call(cbind, egap_list)
  sum_mat   <- do.call(cbind, sum_list)
  colnames(onset_mat) <- colnames(egap_mat) <- colnames(sum_mat) <- tags

  
  # (1) Pairwise swap correlations
  all_rows <- list()
  cmb <- utils::combn(tags, 2, simplify = FALSE)
  for (pair in cmb) {
    t1 <- pair[1]; t2 <- pair[2]
    d_onset <- onset_mat[, t2] - onset_mat[, t1]
    d_gap   <- egap_mat[,  t2] - egap_mat[,  t1]
    swap_corr  <- suppressWarnings(cor(d_onset, d_gap, use = "complete.obs"))

    all_rows[[length(all_rows) + 1]] <- data.frame(
      tag1   = t1,
      tag2   = t2,
      scale1 = scale_from_tag(t1),
      scale2 = scale_from_tag(t2),
      swap_corr = swap_corr,
      stringsAsFactors = FALSE
    )
  }
  pairwise_df <- do.call(rbind, all_rows)
  swap_min <- min(pairwise_df$swap_corr, na.rm = TRUE)
  swap_max <- max(pairwise_df$swap_corr, na.rm = TRUE)

  
  # (2) Across-settings invariance score Inv%
  # For each (patient, condition) row (i.e., each observed diagnosis), compute the variance 
  # of the reconstructed diagnosis time across all prior settings.
  # Gives a vector showing how much the model's prediction changes across priors, 
  # per individual diagnosis.
  row_var_sum <- apply(sum_mat, 1, stats::var, na.rm = TRUE)  # Var_k(S_im^(k))
  inv_percent <- 100 * (1 - mean(row_var_sum, na.rm = TRUE) / vdt)         

  invariance_df <- data.frame(
    Inv_percent = inv_percent
  )

  if (save_csv) {
    utils::write.csv(pairwise_df[, c("tag1","tag2","scale1","scale2","swap_corr")],
                     file.path(directory, "ident_swaps_all_pairs.csv"), row.names = FALSE)
    utils::write.csv(invariance_df,
                     file.path(directory, "ident_invariance_summary.csv"), row.names = FALSE)
  }

  # SUmmary
  cat("\nACROSS-SETTINGS INVARIANCE (Inv%)\n")
  cat(sprintf("Inv%% = %.2f%%\n", inv_percent))

  cat("\nPAIRWISE SWAP CORRELATIONS ρ(ΔO, ΔD)\n")
  to_print <- pairwise_df[, c("scale1","scale2","swap_corr")]
  to_print$swap_corr <- round(to_print$swap_corr, 3)
  to_print <- to_print[order(to_print$scale1, to_print$scale2), ]
  print(to_print, row.names = FALSE)
  cat(sprintf("\nRange of swap correlations: [%.3f, %.3f]\n", swap_min, swap_max))

  invisible(list(
    pairwise_swap_corr = pairwise_df,
    invariance = invariance_df
  ))
}

# Run
identifiability_swap_grid()
