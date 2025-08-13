expected_LTHC_t_after_tau_d <- function(parameters, hyperparameters, tau, M,
                                      mu0, sigma20) {
  # parameters: list(phi, eta, varpi) from VB_gaussian_predictive_density()
  # hyperparameters: list(theta, a, b, u, v, alpha, beta)
  # tau: left-truncation point (e.g., current/cut age)
  # mu0, sigma20: delay prior (length M), same ones used in training

  phi   <- parameters[[1]]   
  eta   <- parameters[[2]]   

  u     <- hyperparameters[[4]]  #(onset location)
  v     <- hyperparameters[[5]]  
  alpha <- hyperparameters[[6]]  
  beta  <- hyperparameters[[7]]  

  K <- ncol(u)

#   # Observed-time Y = onset + delay
#   # onset variance (Student-t scale^2) + delay variance
#   sig_onset <- (beta * (v + 1)) / (alpha * v)                 
#   mu_gap    <- matrix(mu0,     nrow = M, ncol = K)
#   sig_gap2  <- matrix(sigma20, nrow = M, ncol = K)

#   muY <- u + mu_gap                                            
#   sdY <- sqrt(pmax(sig_onset + sig_gap2, 1e-12))               
#   df  <- 2 * alpha   

  # Observed-time Y = onset + delay
  # onset variance (Student-t scale^2) + *truncated* delay variance
  sig_onset <- (beta * (v + 1)) / (alpha * v)
  sigma0    <- sqrt(pmax(sigma20, 1e-12))
  alpha_del <- -mu0 / sigma0
  lambda <- exp(dnorm(alpha_del, log = TRUE) -
                pnorm(alpha_del, lower.tail = FALSE, log.p = TRUE))
  delta  <- lambda * (lambda - alpha_del)
  mu_gap_trunc   <- mu0 + sigma0 * lambda
  sig_gap2_trunc <- pmax(sigma20 * (1 - delta), 1e-12)

  mu_gap  <- matrix(mu_gap_trunc,   nrow = M, ncol = K)
  sig_gap2<- matrix(sig_gap2_trunc, nrow = M, ncol = K)

  muY <- u + mu_gap
  sdY <- sqrt(pmax(sig_onset + sig_gap2, 1e-12))
  df  <- 2 * alpha


  # E[Y | Y > tau] per (m,k) using left-truncated Student-t expectation
  Et_mk <- expected_lst_lefttrunc(df, muY, sdY, tau)           

  # Mix across clusters GIVEN event after tau:
  # weights w_{m,k} proportional to phi_k * eta_{m,k}
  w <- eta * matrix(phi, nrow = M, ncol = K, byrow = TRUE)     
  sw <- rowSums(w)
  sw[sw <= 1e-12] <- 1 # avoid divide-by-zero
  w_norm <- w / sw

  Et <- rowSums(w_norm * Et_mk)                                

  return(Et)
}
