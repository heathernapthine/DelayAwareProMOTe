expected_lst_lefttrunc <- function(df, mu, sigma, tau){
  # #Purpose: 
  # A helper function that calculates the expected value of a left truncated location scale t distribution
  #
  # #Inputs:
  # df - degrees of freedom of the t distribution
  # mu - location parameter of the t distribution
  # sigma - scale parameter of the t distribution
  # tau - the left truncation point
  #
  # #Outputs
  # u - the expected value
  
  taudot <- (tau - mu)/sigma
  
  udot <- ((df/(df-1)) * exp(lgamma((df+1)/2) - lgamma(df/2))/sqrt(df*pi) * ((df+(taudot**2))/df)**((1-df)/2))/pt(taudot, df, lower.tail = FALSE)
  
  u <- (udot * sigma) + mu
  
  return(u)
}
