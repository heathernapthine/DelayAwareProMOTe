VB_gaussian_update_d <- function(
  d, t_obs, rho, tau, iota,
  hyperparameters,
  initial_Cstar, initial_Dstar, initial_pstar, initial_qstar, initial_rstar,
  N, M, K, epsilon, mu0, sigma20,
  sex = NULL, birth_conds = NULL, male_conds = NULL, female_conds = NULL, cond_list
) {

  # Unpack hyperparameters.
  theta_hyper <- hyperparameters[[1]]
  a_hyper     <- hyperparameters[[2]]
  b_hyper     <- hyperparameters[[3]]
  u_hyper     <- hyperparameters[[4]]
  v_hyper     <- hyperparameters[[5]]
  alpha_hyper <- hyperparameters[[6]]
  beta_hyper  <- hyperparameters[[7]]

  # Initialise variational parameters for globals.
  theta_star <- rep(Inf, K)
  a_star     <- matrix(Inf, M, K)
  b_star     <- matrix(Inf, M, K)
  u_star     <- matrix(Inf, M, K)
  v_star     <- matrix(Inf, M, K)
  alpha_star <- matrix(Inf, M, K)
  beta_star  <- matrix(Inf, M, K)

  # Initialise patient level latent variables.
  C_star <- initial_Cstar
  D_star <- initial_Dstar
  p_star <- initial_pstar
  q_star <- initial_qstar
  r_star <- initial_rstar

  # Initialise delay variational parameters.
  gap_mu_star     <- matrix(mu0, N, M, byrow = TRUE)
  gap_sigma2_star <- matrix(sigma20, N, M, byrow = TRUE)
  prev_expected_t_onset <- t_obs - matrix(mu0, N, M, byrow = TRUE)

  # Build censoring masks.
  left_censored  <- t_obs < rho & d == 1
  right_censored <- d == 0 & iota == 0
  not_censored   <- (t_obs >= rho & d == 1) | (d == 0 & iota == 1)

  # Override masks for birth only conditions.
  if (!is.null(birth_conds)) {
    left_censored[,  birth_conds] <- FALSE
    right_censored[, birth_conds] <- FALSE
    not_censored[,   birth_conds] <- TRUE
  }
  # Override masks for sex specific conditions.
  if (!is.null(sex)) {
    female <- which(sex == "Female")
    male   <- which(sex == "Male")
    if (!is.null(male_conds)) {
      left_censored[female,  male_conds] <- FALSE
      right_censored[female, male_conds] <- FALSE
      not_censored[female,   male_conds] <- TRUE
    }
    if (!is.null(female_conds)) {
      left_censored[male,  female_conds] <- FALSE
      right_censored[male, female_conds] <- FALSE
      not_censored[male,   female_conds] <- TRUE
    }
  }

  # Initialise loop controls.
  param_difference <- Inf
  n_steps <- 0

  # Initialise responsibilities and expected presence.
  expected_z <- exp(-C_star); expected_z <- expected_z / rowSums(expected_z)
  expected_d <- D_star / (1 + D_star)
  expected_d <- ((1 - (d == 0 & iota == 0)) * d) + ((d == 0 & iota == 0) * expected_d)

  # Allocate ELBO storage.
  elbos <- rep(NA_real_, 20000)

  # Precompute ELBO prior constants.
  elbo_gamma_prior_const <- lgamma(sum(theta_hyper)) - sum(lgamma(theta_hyper))
  elbo_pi_prior_const    <- sum(lgamma(a_hyper + b_hyper)) - sum(lgamma(a_hyper)) - sum(lgamma(b_hyper))
  elbo_musig_prior_const <- 0.5 * sum(log(v_hyper)) - 0.5 * sum(log(2 * pi)) +
                            sum(alpha_hyper * log(beta_hyper)) - sum(lgamma(alpha_hyper))

  while (param_difference > 0.1) {
    n_steps <- n_steps + 1

    # Update gap moments under truncated normal and derive onset time.
    alpha_tr <- (0 - gap_mu_star) / sqrt(gap_sigma2_star)
    phi_a    <- dnorm(alpha_tr)
    Phi_a    <- pnorm(alpha_tr)
    Z_trunc  <- pmax(1e-10, 1 - Phi_a)

    E_gap   <- gap_mu_star + sqrt(gap_sigma2_star) * (phi_a / Z_trunc)
    Var_gap <- gap_sigma2_star * (1 - (alpha_tr * phi_a / Z_trunc) - (phi_a / Z_trunc)^2)

    t_onset <- pmax(t_obs - E_gap, 0)

    # Compute expected diagnosis time and its square under censoring.
    safe_r_star <- pmax(r_star, 1e-8)

    logdnormtau   <- dnorm((tau - q_star / (2 * safe_r_star)) / (1 / sqrt(2 * safe_r_star)), log = TRUE)
    logpnormctau  <- pnorm((tau - q_star / (2 * safe_r_star)) / (1 / sqrt(2 * safe_r_star)), lower.tail = FALSE, log.p = TRUE)
    logdnormrho   <- dnorm((rho - q_star / (2 * safe_r_star)) / (1 / sqrt(2 * safe_r_star)), log = TRUE)
    logpnormrho   <- pnorm((rho - q_star / (2 * safe_r_star)) / (1 / sqrt(2 * safe_r_star)), log.p = TRUE)

    expected_t_right <- q_star / (2 * safe_r_star) + (1 / sqrt(2 * safe_r_star)) * exp(logdnormtau  - logpnormctau)
    expected_t_left  <- q_star / (2 * safe_r_star) - (1 / sqrt(2 * safe_r_star)) * exp(logdnormrho  - logpnormrho)
    expected_t <- (not_censored * t_onset) + (right_censored * expected_t_right) + (left_censored * expected_t_left)
    expected_t[is.nan(expected_t)] <- 0

    expected_dt <- expected_d * expected_t

    var_t_right  <- (1 / (2 * safe_r_star)) * (1 + ((tau - q_star / (2 * safe_r_star)) / (1 / sqrt(2 * safe_r_star))) * exp(logdnormtau - logpnormctau) - exp(logdnormtau - logpnormctau)^2)
    var_t_left   <- (1 / (2 * safe_r_star)) * (1 - ((rho - q_star / (2 * safe_r_star)) / (1 / sqrt(2 * safe_r_star))) * exp(logdnormrho - logpnormrho) - exp(logdnormrho - logpnormrho)^2)
    expected_t2_right <- var_t_right + expected_t_right^2
    expected_t2_left  <- var_t_left  + expected_t_left^2
    expected_t2 <- (not_censored * t_onset^2) + (right_censored * expected_t2_right) + (left_censored * expected_t2_left)
    expected_t2[is.nan(expected_t2)] <- 0

    expected_dt2 <- expected_d * expected_t2

    # Compute sufficient statistics for global updates.
    EdEz     <- t(expected_d)     %*% expected_z
    EdcEz    <- t(1 - expected_d) %*% expected_z
    EdtEz    <- t(expected_dt)    %*% expected_z
    Edt_sqEz <- t(expected_dt2)   %*% expected_z

    # Save parameters for convergence tracking.
    theta_star_old <- theta_star
    a_star_old     <- a_star
    b_star_old     <- b_star
    u_star_old     <- u_star
    v_star_old     <- v_star
    alpha_star_old <- alpha_star
    beta_star_old  <- beta_star
    gap_mu_star_old     <- gap_mu_star
    gap_sigma2_star_old <- gap_sigma2_star

    # Update global posteriors under conjugacy.
    theta_star <- theta_hyper + colSums(expected_z)
    a_star     <- a_hyper + EdEz
    b_star     <- b_hyper + EdcEz
    u_star     <- ((u_hyper * v_hyper) + EdtEz) / (v_hyper + EdEz)
    v_star     <- v_hyper + EdEz
    alpha_star <- alpha_hyper + (EdEz / 2)
    beta_star <- beta_hyper + 0.5 * (
      Edt_sqEz + v_hyper * (u_hyper^2) -
      ((v_hyper * u_hyper + EdtEz)^2) / (v_hyper + EdEz)
    )

    # Compute expectations under updated global posteriors.
    expected_log_gamma         <- digamma(theta_star) - digamma(sum(theta_star))
    expected_log_pi            <- digamma(a_star) - digamma(a_star + b_star)
    expected_log_pi_complement <- digamma(b_star) - digamma(a_star + b_star)
    expected_log_sigma_squared <- log(beta_star) - digamma(alpha_star)
     #Update expected_mu2/sigma2
    expected_mu2_sigma2 <- (1/v_star) + (u_star**2) * alpha_star/beta_star
    expected_mu_sigma2         <- u_star * alpha_star / beta_star
    expected_sigma2_inverse    <- alpha_star / beta_star

    # Update patient level latents given expectations.
    C_star <- -matrix(expected_log_gamma, N, K) -
      expected_d %*% expected_log_pi - (1 - expected_d) %*% expected_log_pi_complement +
      0.5 * expected_d %*% expected_log_sigma_squared +
      0.5 * expected_d %*% expected_mu2_sigma2 +
      0.5 * expected_dt2 %*% expected_sigma2_inverse -
      expected_dt %*% expected_mu_sigma2

    p_star <- expected_z %*% t(expected_log_pi - expected_log_pi_complement - 0.5 * expected_mu2_sigma2 - 0.5 * expected_log_sigma_squared - 0.5 * log(2 * pi))
    q_star <- expected_z %*% t(expected_mu_sigma2)
    r_star <- 0.5 * expected_z %*% t(expected_sigma2_inverse)
    D_star <- exp(p_star + q_star * expected_t - r_star * expected_t2)

    # Refresh responsibilities and expected presence.
    expected_z <- exp(-C_star); expected_z <- expected_z / rowSums(expected_z)
    expected_d <- D_star / (1 + D_star)
    expected_d <- ((1 - right_censored) * d) + (right_censored * expected_d)

    # Update delay posteriors by blending prior and observed gaps.
    for (j in 1:M) {
      for (i in 1:N) {
        if (!is.na(t_obs[i, j]) && d[i, j] == 1) {
          t_onset_ij <- expected_t[i, j]  # Current estimate of onset.
          observed_delay <- max(t_obs[i, j] - t_onset_ij, 0.01)
          prior_mean <- mu0[j]; prior_var <- sigma20[j]; prior_precision <- 1 / prior_var
          likelihood_precision <- 10.0
          post_precision <- likelihood_precision + prior_precision
          post_var  <- 1 / post_precision
          post_mean <- post_var * (likelihood_precision * observed_delay + prior_precision * prior_mean)
          gap_mu_star[i, j]     <- max(post_mean, 0.01)
          gap_sigma2_star[i, j] <- max(post_var,  0.01)
        }
      }
    }

    # Compute convergence measure across globals and gap parameters.
    if (n_steps > 1) {
      param_difference <-
        mean(abs(theta_star / sum(theta_star) - theta_star_old / sum(theta_star_old))) +
        mean(abs(a_star / (a_star + b_star) - a_star_old / (a_star_old + b_star_old))) +
        mean(abs(u_star - u_star_old)) +
        mean(abs(alpha_star - alpha_star_old)) +
        mean(abs(beta_star  - beta_star_old)) +
        mean(abs(gap_mu_star - gap_mu_star_old)) +
        mean(abs(gap_sigma2_star - gap_sigma2_star_old))
    }

    prev_expected_t_onset <- expected_t

    # Quick progress diagnostic.
    cat(n_steps, param_difference, "\n")

    #### =========================
    #### ELBO (moment-matched; aligned with updates)
    #### =========================
    log2pi <- log(2 * pi)

    # Priors over globals (same forms as before).
    elbo_gamma_prior <- elbo_gamma_prior_const + sum((theta_hyper - 1) * expected_log_gamma)
    elbo_pi_prior    <- elbo_pi_prior_const +
                        sum((a_hyper - 1) * expected_log_pi) +
                        sum((b_hyper - 1) * expected_log_pi_complement)
    elbo_musig_prior <- elbo_musig_prior_const -
                        sum((alpha_hyper + 3/2) * expected_log_sigma_squared) -
                        sum(beta_hyper * expected_sigma2_inverse) -
                        sum((v_hyper / 2) * expected_mu2_sigma2) +
                        sum(u_hyper * v_hyper * expected_mu_sigma2) -
                        sum((((u_hyper^2) * v_hyper) / 2) * expected_sigma2_inverse)

    # Likelihood under the same quadratic surrogate used in the updates.
    ones_KM <- matrix(1, K, M)
    elbo_likelihood <-
      sum(expected_z %*% expected_log_gamma) +
      sum(expected_d * (expected_z %*% t(expected_log_pi))) +
      sum((1 - expected_d) * (expected_z %*% t(expected_log_pi_complement))) -
      0.5 * log2pi * sum(expected_d * (expected_z %*% ones_KM)) -
      0.5 * sum(expected_d  * (expected_z %*% t(expected_log_sigma_squared))) -
      0.5 * sum(expected_d  * (expected_z %*% t(expected_mu2_sigma2))) +
      sum(expected_dt  * (expected_z %*% t(expected_mu_sigma2))) -
      0.5 * sum(expected_dt2 * (expected_z %*% t(expected_sigma2_inverse)))

    # Entropies for globals.
    elbo_qgamma <- lgamma(sum(theta_star)) - sum(lgamma(theta_star)) + sum((theta_star - 1) * expected_log_gamma)
    elbo_qpi    <- sum(lgamma(a_star + b_star)) - sum(lgamma(a_star)) - sum(lgamma(b_star)) +
                   sum((a_star - 1) * expected_log_pi) + sum((b_star - 1) * expected_log_pi_complement)
    elbo_qmusig <- 0.5 * sum(log(v_star)) - 0.5 * sum(log(2 * pi)) +
                   sum(alpha_star * log(beta_star)) - sum(lgamma(alpha_star)) -
                   sum((alpha_star + 3/2) * expected_log_sigma_squared) -
                   sum(beta_star * expected_sigma2_inverse) -
                   sum((v_star / 2) * expected_mu2_sigma2) +
                   sum(u_star * v_star * expected_mu_sigma2) -
                   sum((((u_star^2) * v_star) / 2) * expected_sigma2_inverse)

    # Entropy of responsibilities and presence variables.
    expected_z_safe <- pmax(expected_z, 1e-15)
    elbo_qz <- -sum(expected_z * log(expected_z_safe))
    expected_d_safe <- pmax(pmin(expected_d, 1 - 1e-15), 1e-15)
    elbo_qd <- 0
    if (sum(right_censored) > 0) {
      elbo_qd <- -sum(expected_d[right_censored] * log(expected_d_safe[right_censored]) +
                      (1 - expected_d[right_censored]) * log(1 - expected_d_safe[right_censored]))
    }

    # KL for delay posteriors vs priors (diagnosed entries only; fixed mask).
    diagnosed_mask <- (!is.na(t_obs)) & (d == 1)
    mu0_mat     <- matrix(mu0,     nrow = N, ncol = M, byrow = TRUE)
    sigma20_mat <- matrix(sigma20, nrow = N, ncol = M, byrow = TRUE)
    gap_sigma2_star_safe <- pmax(gap_sigma2_star, 1e-8)

    kl_qgap_pgap <- 0.5 * sum(
      log(sigma20_mat[diagnosed_mask]) - log(gap_sigma2_star_safe[diagnosed_mask]) - 1 +
      (gap_sigma2_star[diagnosed_mask] / sigma20_mat[diagnosed_mask]) +
      ((gap_mu_star[diagnosed_mask] - mu0_mat[diagnosed_mask])^2 / sigma20_mat[diagnosed_mask])
    )

    # Full ELBO.
    elbo <- elbo_gamma_prior +
            elbo_pi_prior +
            elbo_musig_prior +
            elbo_likelihood -
            elbo_qgamma -
            elbo_qpi -
            elbo_qmusig -
            elbo_qz -
            elbo_qd -
            kl_qgap_pgap

    elbos[n_steps] <- elbo
  }

  # Return posterior parameters and diagnostics.
  return(list(
    posterior.parameters = list(
      # Globals.
      theta_star = theta_star,
      a_star = a_star, b_star = b_star,
      u_star = u_star, v_star = v_star,
      alpha_star = alpha_star, beta_star = beta_star,

      # Rename to match notation in main script.
      gamma_alpha = theta_star,
      pi_a = a_star, pi_b = b_star,
      mu_u = u_star, mu_alpha = alpha_star, mu_beta = beta_star,

      # Patient level.
      C_star = C_star, p_star = p_star, q_star = q_star, r_star = r_star, D_star = D_star,

      # Diagnostics.
      gap_mu_star = gap_mu_star,
      gap_sigma2_star = gap_sigma2_star,
      expected_t = expected_t, expected_d = expected_d
    ),
    n_steps = n_steps,
    final_step_size = param_difference,
    elbo = elbos[1:n_steps],
    cond_list = cond_list
  ))
}
