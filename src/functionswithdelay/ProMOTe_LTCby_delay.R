probability_LTHC_by_T_d <- function(parameters, hyperparameters, T, tau, M, 
                                  mu0, sigma20) {

# Purpose: Predict P(diagnosed in (tau, T]] for each condition under delay-aware model.
# Inputs: parameters = list(phi, eta); hyperparameters (u,v,alpha,beta);
#   T, tau, M; delay priors mu0, sigma20.
# Outputs: prob_m (length M) = probability of diagnosis by T given survival past tau,
#   marginalised over clusters using phi.

  phi   <- parameters[[1]]          
  eta   <- parameters[[2]]            

  u     <- hyperparameters[[4]]       
  v     <- hyperparameters[[5]]       
  alpha <- hyperparameters[[6]]       
  beta  <- hyperparameters[[7]]       

  K <- ncol(u)

  # Predictive for observed time Y = onset + delay (with *truncated* delay moments)
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

  muY  <- u + mu_gap
  sdY  <- sqrt(pmax(sig_onset + sig_gap2, 1e-12))

  # CDFs under location-scale Student-t with df = 2*alpha
  Ftau <- matrix(pt((tau - muY) / sdY, df = 2 * alpha), nrow = M)
  FT   <- matrix(pt((T   - muY) / sdY, df = 2 * alpha), nrow = M)

  # Conditional CDF for Y in (tau, T] given Y > tau
  den <- pmax(1 - Ftau, 1e-12)
  GT  <- (FT - Ftau) / den
  GT  <- pmin(pmax(GT, 0), 1) # clamp into [0,1]

  # If eta==0 (no chance of presence by tau under the cluster), GT irrelevant
  GT[eta == 0] <- 1

  # Mix over clusters with phi
  prob_mk <- GT * eta                                         
  prob_m  <- prob_mk %*% matrix(phi, ncol = 1)               

  return(as.numeric(prob_m))
}
