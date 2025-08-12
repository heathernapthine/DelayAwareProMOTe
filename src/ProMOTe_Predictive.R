VB_gaussian_predictive_density <- function(hyperparameters, M_obs, M_part, M_unobs, d_obs, t_obs, d_part, rho, tau, M) {
  # #Purpose: 
  # Estimate the parameters of the posterior predictive distribution for a new individual
  #
  # #Inputs:
  # hyperparameters - a list of hyperparameters of the posterior
  # M_obs - indices of the individual's fully observed conditions
  # M_part - indices of the individual's partially observed conditions
  # M_unobs - indices of the individual's unobserved conditions
  # d_obs - a vector containing the absence/presence of the fully observed conditions as 0s and 1s
  # t_obs - a vector containing the onset times of the fully observed conditions
  # d_part - a vector containing the absence/presence of the partially observed conditions as 0s and 1s
  # rho - a length N vector containing observation start dates
  # tau - a length N vector containing observation end dates
  # M - the number of conditions in the data
  #
  # #Outputs
  # phi - a vector of cluster probabilities
  # eta - a matrix of condition probabilities given cluster
  # varpi - a vector of useful intermediate condition probabilities conditional on clusters
  
  theta <- hyperparameters[[1]]
  a <- hyperparameters[[2]]
  b <- hyperparameters[[3]]
  u <- hyperparameters[[4]]
  v <- hyperparameters[[5]]
  alpha <- hyperparameters[[6]]
  beta <- hyperparameters[[7]]
  
  n_M_obs <- length(M_obs)
  n_M_part <- length(M_part)
  n_M_unobs <- M - n_M_obs - n_M_part
  
  varpi <- a/(a + b)
  
  sig <- (beta*(v+1))/(alpha*v)
  
  Ftau <- matrix(plst(tau, 2*alpha, u, sqrt(sig)), nrow = M)
  Frho <- matrix(plst(rho, 2*alpha, u, sqrt(sig)), nrow = M)
  
  A <- (varpi * (1 - Ftau)) + 1 - varpi
  B <- varpi * Frho
  
  eta <- (varpi * (1 - Ftau))/A
  
  if (length(A[M_unobs,]) == length(theta)) {
    psi <- theta*A[M_unobs,]
  } else {
    psi <- theta*apply(A[M_unobs,], 2, prod)
  }
  
  if (length(B[M_part,]) == length(theta)) {
    psi <- psi*B[M_part,]
  } else {
    psi <- psi*apply(B[M_part,], 2, prod)
  }
  
  C <- (varpi[M_obs,] * (dlst(t_obs, 2*alpha[M_obs,], u[M_obs,], sig[M_obs,])))**d_obs
  
  C <- C * (1-varpi[M_obs,])**(1-d_obs)
  
  if (length(C) == length(psi)){
    phi <- psi * C
  } else {
    phi <- psi * apply(C, 2, prod)
  }
  
  phi <- phi/sum(phi)
  
  return(list(phi = phi, eta = eta, varpi = varpi))
  
}