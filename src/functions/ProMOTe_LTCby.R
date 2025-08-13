probability_LTHC_by_T <- function(parameters, hyperparameters, T, tau, M){
  # #Purpose: 
  # Estimate the probability of conditions occurring in an individual by time T
  #
  # #Inputs:
  # parameters - a list of posterior predictive parameters (the output of VB_gaussian_predictive_density)
  # hyperparameters - a list of hyperparameters of the posterior
  # T - the time by which the LTCs should occur
  # tau - the individual's age at end of observation
  # M - the number of conditions
  #
  # #Outputs
  # prob - a vector of probabilities
  
  phi <- parameters[[1]]
  eta <- parameters[[2]]
  
  u <- hyperparameters[[4]]
  v <- hyperparameters[[5]]
  alpha <- hyperparameters[[6]]
  beta <- hyperparameters[[7]]
  
  sig <- (beta*(v+1))/(alpha*v)
  
  Ftau <- matrix(plst(tau, 2*alpha, u, sqrt(sig)), nrow = M)
  
  FT <- matrix(plst(T, 2*alpha, u, sqrt(sig)), nrow = M)
  
  GT <- (FT-Ftau)/(1-Ftau)
  
  GT[eta == 0] <- 1
  
  prob <- GT * eta
  
  prob <- prob %*% matrix(phi, ncol = 1)
  
  return(c(prob))
  
}