expected_LTHC_t_after_tau <- function(parameters, hyperparameters, tau, M){
  # #Purpose: 
  # Estimate the expected time of occurrence of conditions in an individual after tau
  #
  # #Inputs:
  # parameters - a list of posterior predictive parameters (the output of VB_gaussian_predictive_density)
  # hyperparameters - a list of hyperparameters of the posterior
  # tau - the individual's age at end of observation
  # M - the number of conditions
  #
  # #Outputs:
  # Et - Expected onset times
  
  phi <- parameters[[1]]
  eta <- parameters[[2]]
  
  u <- hyperparameters[[4]]
  v <- hyperparameters[[5]]
  alpha <- hyperparameters[[6]]
  beta <- hyperparameters[[7]]
  
  df <- 2*alpha
  mu <- u
  sigma <- sqrt((beta*(v+1))/(alpha*v))
  
  Et <- expected_lst_lefttrunc(df, mu, sigma, tau)
  
  Et <- c(phi %*% t(Et))
  
  return(Et)
  
}
