VB_gaussian_update <- function(
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


  while (param_difference > 0.01 && n_steps < 100) {
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
    EdEz     <- t(expected_d)   %*% expected_z
    EdcEz    <- t(1 - expected_d) %*% expected_z
    EdtEz    <- t(expected_dt)  %*% expected_z
    Edt_sqEz <- t(expected_dt2) %*% expected_z

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
    expected_log_gamma        <- digamma(theta_star) - digamma(sum(theta_star))
    expected_log_pi           <- digamma(a_star) - digamma(a_star + b_star)
    expected_log_pi_complement<- digamma(b_star) - digamma(a_star + b_star)
    expected_log_sigma_squared<- log(beta_star) - digamma(alpha_star)
    expected_mu2_sigma2 <- ((1 / v_star) + (u_star^2)) * (alpha_star / beta_star)
    expected_mu_sigma2        <- u_star * alpha_star / beta_star
    expected_sigma2_inverse   <- alpha_star / beta_star

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

    # Print quick progress diagnostics.
    cat(n_steps, param_difference, "\n")

    }


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
      expected_t = expected_t, expected_d = expected_d
    ),
    n_steps = n_steps,
    final_step_size = param_difference,
    cond_list = cond_list
  ))
}
