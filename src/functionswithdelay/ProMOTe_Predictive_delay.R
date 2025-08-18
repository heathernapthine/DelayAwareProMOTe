VB_gaussian_predictive_density_d <- function(
  hyperparameters, 
  M_obs, M_part, M_unobs, 
  d_obs, t_obs, d_part, 
  rho, tau, M,
  mu0, sigma20   # New for delay model: delay prior (length M)
) {

 # Purpose: Posterior cluster responsibilities and predictive pieces from partial observations
#   under delay-aware model (Y = onset + gap with truncated-normal delay).
# Inputs: hyperparameters (theta,a,b,u,v,alpha,beta);
#   M_obs, M_part, M_unobs; d_obs, t_obs, d_part; rho, tau, M; mu0, sigma20.
# Outputs: list(phi = posterior cluster probs,
#               eta = P(event after tau | cluster, condition),
#               varpi = P(condition present | cluster)).


  # Unpack posteriors (M x K except theta)
  theta <- hyperparameters[[1]]
  a     <- hyperparameters[[2]]
  b     <- hyperparameters[[3]]
  u     <- hyperparameters[[4]]
  v     <- hyperparameters[[5]]
  alpha <- hyperparameters[[6]]
  beta  <- hyperparameters[[7]]

 
  theta <- as.numeric(theta)
  theta <- pmax(theta, 1e-12); theta <- theta / sum(theta)

  # Presence probability per condition/cluster
  varpi <- a / (a + b)  

  # Delay-aware location/scale for observed time Y = onset + delay
  # onset predictive variance
  sig_onset <- (beta * (v + 1)) / (alpha * v)

  # broadcast *truncated* delay moments across clusters
  sigma0    <- sqrt(pmax(sigma20, 1e-12))
  alpha_del <- -mu0 / sigma0
  # stable lambda(alpha) = φ(alpha) / (1 - Φ(alpha))
  lambda <- exp(dnorm(alpha_del, log = TRUE) -
                pnorm(alpha_del, lower.tail = FALSE, log.p = TRUE))
  delta  <- lambda * (lambda - alpha_del)

  mu_gap_trunc   <- mu0 + sigma0 * lambda
  sig_gap2_trunc <- pmax(sigma20 * (1 - delta), 1e-12)

  mu_gap  <- matrix(mu_gap_trunc,   nrow = M, ncol = ncol(u))
  sig_gap2<- matrix(sig_gap2_trunc, nrow = M, ncol = ncol(u))

  muY   <- u + mu_gap
  sigY2 <- sig_onset + sig_gap2
  sdY   <- sqrt(pmax(sigY2, 1e-12)) # guard

  # CDFs for window edges under Y ~ t_{2*alpha}(muY, sdY)
  Ftau <- matrix(pt((tau - muY) / sdY, df = 2 * alpha), nrow = M)
  Frho <- matrix(pt((rho - muY) / sdY, df = 2 * alpha), nrow = M)

  # A: contribution of unobserved conditions
  # A = P(no event by tau | cluster) * varpi + P(absent) = varpi*(1-Ftau) + (1-varpi)
  A <- (varpi * (1 - Ftau)) + (1 - varpi)

  # B: left-censored presences (event before rho): B = varpi * Frho
  B <- varpi * Frho

  # eta: P(event after tau | cluster, condition)
  eta <- (varpi * (1 - Ftau)) / A

  # Build psi (product over conditions without explicit times)
  # Start from cluster prior theta, then multiply A over M_unobs and B over M_part
  if (length(M_unobs) > 0) {
    Au <- A[M_unobs, , drop = FALSE]
    psi <- if (nrow(Au) == 1) theta * Au else theta * apply(Au, 2, prod)
  } else {
    psi <- theta
  }

  if (length(M_part) > 0) {
    Bp <- B[M_part, , drop = FALSE]
    psi <- if (nrow(Bp) == 1) psi * Bp else psi * apply(Bp, 2, prod)
  }

  # Fully observed conditions with times in [rho, tau]
  # For d_obs=1 we use the t-density of Y at t_obs; for d_obs=0 we use (1 - varpi)
  C <- 1
  if (length(M_obs) > 0) {
    Mobs <- M_obs
    dob  <- as.numeric(d_obs)
    tobs <- as.numeric(t_obs)

    # density for location–scale t: dt((x-mu)/sd, df)/sd
    dens <- dt((tobs - muY[Mobs, , drop = FALSE]) / sdY[Mobs, , drop = FALSE],
               df = 2 * alpha[Mobs, , drop = FALSE]) / sdY[Mobs, , drop = FALSE]

    term_present <- (varpi[Mobs, , drop = FALSE] * dens) ^ dob
    term_absent  <- (1 - varpi[Mobs, , drop = FALSE]) ^ (1 - dob)

    C_mat <- term_present * term_absent
    C <- if (nrow(C_mat) == 1) C_mat else apply(C_mat, 2, prod)
  }

  # Posterior cluster probabilities
  phi <- psi * C
  s <- sum(phi)
  if (!is.finite(s) || s <= 0) {
    phi <- rep(1 / length(theta), length(theta)) 
  } else {
    phi <- phi / s
  }

  list(phi = phi, eta = eta, varpi = varpi)
}
