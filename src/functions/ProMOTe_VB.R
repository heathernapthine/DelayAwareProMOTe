VB_gaussian_update <- function(d, t, rho, tau, iota, hyperparameters, initial_Cstar, initial_Dstar, initial_pstar, initial_qstar, initial_rstar, N, M, K, epsilon, sex = NULL, birth_conds = NULL, male_conds = NULL, female_conds = NULL, cond_list) {
  # #Purpose: 
  # Perform Variational Bayes (VB) updates for the Gaussian latent class model with fixed K and censored data
  #
  # #Inputs:
  # d - a NxM matrix containing data about presence of conditions
  # t - a NxM matrix containing data about onset times of conditions
  # rho - a length N vector containing study start dates
  # tau - a length N vector containing study end dates
  # iota - a length N vector indicating if individuals are alive/deceased at time tau
  # hyperparameters - a list of hyperparameters of the prior
  # initial_Cstar - a NxK matrix containing an initial value to initialise the latent variable z
  # initial_Dstar - a NxM matrix containing an initial value to initialise the latent variable d
  # initial_pstar - a NxM matrix containing an initial value to initialise the latent variable d
  # initial_qstar - a NxM matrix containing an initial value to initialise the latent variable t
  # initial_rstar - a NxM matrix containing an initial value to initialise the latent variable t
  # N - the number of individuals in the data
  # K - the number of clusters
  # M - the number of conditions in the data
  # epsilon - a number used to determine the stopping condition
  # sex - the sex of individuals in the data
  # birth_conds - column indices for d indicating conditions which only occur at birth
  # male_conds - column indices for d indicating conditions which only occur in males
  # female_conds - column indices for d indicating conditions which only occur in females
  # cond_list - a string vector containing the names of the conditions
  #
  # #Outputs
  # theta_star - a vector of parameters for the VB posterior of gamma
  # a_star - a MxK matrix of parameters for the VB posterior of pi
  # b_star - a MxK matrix of parameters for the VB posterior of pi
  # u_star - a MxK matrix of parameters for the VB posterior of mu and sigma^2
  # v_star - a MxK matrix of parameters for the VB posterior of mu and sigma^2
  # alpha_star - a MxK matrix of parameters for the VB posterior of mu and sigma^2
  # beta_star - a MxK matrix of parameters for the VB posterior of mu and sigma^2
  # C_star - a NxK matrix of parameters for the VB posterior of z
  # p_star - a NxM matrix of parameters for the VB posterior of d
  # q_star - a NxM matrix of parameters for the VB posterior of t and d
  # r_star - a NxM matrix of parameters for the VB posterior of t and d
  # D_star - a NxM matrix of parameters for the VB posterior of d
  # n_steps - the number of iterations required to achieve the stopping condition
  
  # read parameters
  theta_hyper <- hyperparameters[[1]]
  a_hyper <- hyperparameters[[2]]
  b_hyper <- hyperparameters[[3]]
  u_hyper <- hyperparameters[[4]]
  v_hyper <- hyperparameters[[5]]
  alpha_hyper <- hyperparameters[[6]]
  beta_hyper <- hyperparameters[[7]]
  
  #initialise variables
  theta_star <- rep(Inf, K)
  theta_star_naught <- Inf
  a_star <- matrix(Inf, M, K)
  b_star <- matrix(Inf, M, K)
  u_star <- matrix(Inf, M, K)
  v_star <- matrix(Inf, M, K)
  alpha_star <- matrix(Inf, M, K)
  beta_star <- matrix(Inf, M, K)
  C_star <- initial_Cstar
  D_star <- initial_Dstar
  p_star <- initial_pstar
  q_star <- initial_qstar
  r_star <- initial_rstar
  
  #detect censoring
  left_censored <- t < rho & d==1
  right_censored <- d==0 & iota==0
  not_censored <- (t >= rho & d==1) | (d==0 & iota==1)
  
  #hardcode lifetime conditions
  if (!is.null(birth_conds)){
    left_censored[,birth_conds] <- FALSE
    right_censored[,birth_conds] <- FALSE
    not_censored[,birth_conds] <- TRUE
  }
  
  if(!is.null(sex)) {
    female <- which(sex == "Female")
    if(!is.null(male_conds)){
      left_censored[female, male_conds] <- FALSE
      right_censored[female, male_conds] <- FALSE
      not_censored[female, male_conds] <- TRUE
    }
    if(!is.null(female_conds)){
      male <- which(sex == "Male")
      left_censored[male, female_conds] <- FALSE
      right_censored[male, female_conds] <- FALSE
      not_censored[male, female_conds] <- TRUE
    }
  }
  
  #initialise stopping condition difference
  param_difference <- Inf
  elbo_difference <- Inf
  
  #initialise working variables
  expected_log_gamma <- NA
  expected_log_pi <- NA
  expected_log_pi_complement <- NA
  expected_mu <- NA
  expected_mu_squared <- NA
  expected_log_sigma_squared <- NA
  expected_sigma2_inverse <- NA
  
  n_steps <- 0
  
  #Update expected_z
  expected_z <- exp(-C_star)
  expected_z_row_sums <- rowSums(expected_z)
  expected_z <- expected_z/expected_z_row_sums
  
  #Update expected_d
  expected_d <- D_star
  expected_d <- expected_d/(1+expected_d)
  expected_d <- ((1-right_censored) * d) + (right_censored * expected_d)
  
  logdnormtau <- dnorm((tau - q_star/(2*r_star))/(1/sqrt(2*r_star)), 0, 1, log = TRUE)
  logpnormctau <- pnorm((tau - q_star/(2*r_star))/(1/sqrt(2*r_star)), 0, 1, lower.tail = FALSE, log.p = TRUE)
  logdnormrho <- dnorm((rho - q_star/(2*r_star))/(1/sqrt(2*r_star)), 0, 1, log = TRUE)
  logpnormrho <- pnorm((rho - q_star/(2*r_star))/(1/sqrt(2*r_star)), 0, 1, log.p = TRUE)
  
  #Update expected_t
  expected_t_right <- q_star/(2*r_star) + (1/sqrt(2*r_star))*exp(logdnormtau - logpnormctau)
  expected_t_left <- q_star/(2*r_star) - (1/sqrt(2*r_star))*exp(logdnormrho - logpnormrho)
  expected_t <- (not_censored * t) + (right_censored * expected_t_right) + (left_censored*expected_t_left)
  expected_t[is.nan(expected_t)] <- 0
  
  #Update expected_dt
  expected_dt <- expected_d * expected_t
  
  #Update expected_t^2
  var_t_right <- (1/(2*r_star)) * (1+((tau - q_star/(2*r_star))/(1/sqrt(2*r_star)) * exp(logdnormtau - logpnormctau)) - exp(logdnormtau - logpnormctau)**2)
  var_t_left <- (1/(2*r_star)) * (1-((rho - q_star/(2*r_star))/(1/sqrt(2*r_star)) * exp(logdnormrho - logpnormrho)) - exp(logdnormrho - logpnormrho)**2)
  expected_t2_right <- var_t_right + (expected_t_right)**2
  expected_t2_left <- var_t_left + (expected_t_left)**2
  expected_t2 <- (not_censored * t**2) + (right_censored*expected_t2_right) + (left_censored*expected_t2_left)
  expected_t2[is.nan(expected_t2)] <- 0
  
  #Update expected_dt^2
  expected_dt2 <- expected_d * expected_t2
  
  #Update working variables
  EdEz <- t(expected_d)%*%expected_z
  EdcEz <- t(1 - expected_d)%*%expected_z
  EdtEz <- t(expected_dt)%*%expected_z
  Edt <- expected_dt
  Edt_sqEz <- t(expected_dt2)%*%expected_z
  
  elbo_gamma_prior_const <- lgamma(sum(theta_hyper)) - sum(lgamma(theta_hyper))
  elbo_pi_prior_const <- sum(lgamma(a_hyper+b_hyper)) - sum(lgamma(a_hyper)) - sum(lgamma(b_hyper))
  elbo_musig_prior_const <- 0.5 * sum(log(v_hyper)) - 0.5*sum(log(2*pi)) + sum(alpha_hyper*log(beta_hyper)) - sum(lgamma(alpha_hyper))
  elbo_qt_const <- -0.5*log(2*pi)
  elbos <- rep(NA, 20000)
  
  
  cat("n_steps")
  cat(" ")
  cat("param_difference")
  cat(" ")
  cat("elbo_difference")
  cat(" ")
  cat("elbo")
  cat("\n")
  
  #perform updates
  while(param_difference > epsilon && n_steps < 70) {
    
    n_steps <- n_steps + 1
    
    #Update theta_star
    theta_star_old <- theta_star
    theta_star <- theta_hyper + colSums(expected_z)
    theta_star_naught <- sum(theta_star)
    
    #Update a_star
    a_star_old <- a_star
    a_star <- a_hyper + EdEz
    
    #Update b_star
    b_star_old <- b_star
    b_star <- b_hyper + EdcEz
    
    #Update u_star
    u_star_old <- u_star
    u_star <- ((u_hyper * v_hyper) + EdtEz)/(v_hyper + EdEz)
    
    #Update v_star
    v_star_old <- v_star
    v_star <- v_hyper + EdEz
    
    #Update alpha_star
    alpha_star_old <- alpha_star
    alpha_star <- alpha_hyper + (EdEz/2)
    
    #Update beta_star
    beta_star_old <- beta_star
    beta_star <- beta_hyper + Edt_sqEz/2
    beta_star <- beta_star - (EdtEz**2)/(2*EdEz)
    beta_star <- beta_star + ((v_hyper*(u_hyper**2))*(EdEz))/(2*(v_hyper + EdEz))
    beta_star <- beta_star + ((v_hyper*u_hyper)*(EdtEz))/(v_hyper + EdEz)
    beta_star <- beta_star + (v_hyper*((EdtEz)**2))/(2*(v_hyper + EdEz)*(EdEz))
    
    #Update expected_log_gamma
    expected_log_gamma <- digamma(theta_star) - digamma(theta_star_naught)
    
    #Update expected_log_pi
    expected_log_pi <- digamma(a_star) - digamma(a_star+b_star)
    
    #Update expected_log_pi_complement
    expected_log_pi_complement <- digamma(b_star) - digamma(a_star+b_star)
    
    #Update expected_log_sigma_squared
    expected_log_sigma_squared <- log(beta_star) - digamma(alpha_star)
    
    #Update expected_mu2/sigma2
    expected_mu2_sigma2 <- (1/v_star) + (u_star**2) * alpha_star/beta_star
    
    #Update expected_mu/sigma2
    expected_mu_sigma2 <- u_star*alpha_star/beta_star
    
    #Update expected_sigma2_inverse
    expected_sigma2_inverse <- alpha_star/beta_star
    
    #Update C_star
    C_star_old <- C_star
    C_star <- -matrix(expected_log_gamma, N, K) - expected_d%*%expected_log_pi
    C_star <- C_star - (1-expected_d)%*%expected_log_pi_complement
    C_star <- C_star + 0.5 * expected_d%*%expected_log_sigma_squared
    C_star <- C_star + 0.5 * expected_d%*%expected_mu2_sigma2
    C_star <- C_star + 0.5 * expected_dt2%*%expected_sigma2_inverse
    C_star <- C_star - expected_dt%*%expected_mu_sigma2
    
    #Update p_star #
    p_star_old <- p_star
    p_star <- expected_log_pi - expected_log_pi_complement - 0.5*expected_mu2_sigma2 - 0.5*(log(beta_star) - digamma(alpha_star)) - 0.5*log(2*pi)
    p_star <- expected_z %*% t(p_star)
    
    #Update q_star
    q_star_old <- q_star
    q_star <- expected_z %*% t(expected_mu_sigma2)
    
    #Update r_star
    r_star_old <- r_star
    r_star <- 0.5*expected_z %*% t(expected_sigma2_inverse)
    
    #Update D_star
    D_star_old <- D_star
    
    D_star <- exp(p_star+q_star*expected_t-r_star*expected_t2)
    
    #Update expected_z
    expected_z <- exp(-C_star)#/(1-exp(-C_star))
    expected_z_row_sums <- rowSums(expected_z)
    expected_z <- expected_z/expected_z_row_sums
    
    #Update expected_d
    expected_d <- D_star
    expected_d <- expected_d/(1+expected_d)
    expected_d <- ((1-right_censored) * d) + (right_censored * expected_d)
    
    logdnormtau <- dnorm((tau - q_star/(2*r_star))/(1/sqrt(2*r_star)), 0, 1, log = TRUE)
    logpnormctau <- pnorm((tau - q_star/(2*r_star))/(1/sqrt(2*r_star)), 0, 1, lower.tail = FALSE, log.p = TRUE)
    logdnormrho <- dnorm((rho - q_star/(2*r_star))/(1/sqrt(2*r_star)), 0, 1, log = TRUE)
    logpnormrho <- pnorm((rho - q_star/(2*r_star))/(1/sqrt(2*r_star)), 0, 1, log.p = TRUE)
    
    #Update expected_t
    expected_t_right <- q_star/(2*r_star) + (1/sqrt(2*r_star))*exp(logdnormtau - logpnormctau)
    expected_t_left <- q_star/(2*r_star) - (1/sqrt(2*r_star))*exp(logdnormrho - logpnormrho)
    expected_t <- (not_censored * t) + (right_censored * expected_t_right) + (left_censored*expected_t_left)
    expected_t[is.nan(expected_t)] <- 0
    
    #Update expected_dt
    expected_dt <- expected_d * expected_t
    
    #Update expected_t^2
    var_t_right <- (1/(2*r_star)) * (1+((tau - q_star/(2*r_star))/(1/sqrt(2*r_star)) * exp(logdnormtau - logpnormctau)) - exp(logdnormtau - logpnormctau)**2)
    var_t_left <- (1/(2*r_star)) * (1-((rho - q_star/(2*r_star))/(1/sqrt(2*r_star)) * exp(logdnormrho - logpnormrho)) - exp(logdnormrho - logpnormrho)**2)
    expected_t2_right <- var_t_right + (expected_t_right)**2
    expected_t2_left <- var_t_left + (expected_t_left)**2
    expected_t2 <- (not_censored * t**2) + (right_censored*expected_t2_right) + (left_censored*expected_t2_left)
    expected_t2[is.nan(expected_t2)] <- 0
    
    #Update expected_dt^2
    expected_dt2 <- expected_d * expected_t2
    
    #Update working variables
    EdEz <- t(expected_d)%*%expected_z
    EdcEz <- t(1 - expected_d)%*%expected_z
    EdtEz <- t(expected_dt)%*%expected_z
    Edt <- expected_dt
    Edt_sqEz <- t(expected_dt2)%*%expected_z
    
    if(n_steps>1) elbo_gamma_prior_old <- elbo_gamma_prior
    elbo_gamma_prior <- elbo_gamma_prior_const + sum((theta_hyper-1)*expected_log_gamma)
    if(n_steps>1) elbo_pi_prior_old <- elbo_pi_prior
    elbo_pi_prior <- elbo_pi_prior_const + sum((a_hyper-1)*expected_log_pi) + sum((b_hyper-1)*expected_log_pi_complement)
    if(n_steps>1) elbo_musig_prior_old <- elbo_musig_prior
    elbo_musig_prior <- elbo_musig_prior_const - sum((alpha_hyper+3/2)*expected_log_sigma_squared) - sum(beta_hyper*expected_sigma2_inverse) - sum((v_hyper/2) * expected_mu2_sigma2) + sum(u_hyper*v_hyper*expected_mu_sigma2) - sum((((u_hyper**2)*v_hyper)/2) * expected_sigma2_inverse)
    if(n_steps>1) elbo_likelihood_old <- elbo_likelihood

    log2pi <- log(2 * pi)
    elbo_likelihood <- sum(expected_z %*% expected_log_gamma) +
      sum(expected_d * (expected_z %*% t(expected_log_pi))) +
      sum((1 - expected_d) * (expected_z %*% t(expected_log_pi_complement))) -
      0.5 * log2pi * sum(expected_d * (expected_z %*% matrix(1, K, M))) -
      0.5 * sum(expected_d * (expected_z %*% t(expected_log_sigma_squared))) -
      0.5 * sum(expected_d * (expected_z %*% t(expected_mu2_sigma2))) +
      sum(expected_dt * (expected_z %*% t(expected_mu_sigma2))) -
      0.5 * sum(expected_dt2 * (expected_z %*% t(expected_sigma2_inverse)))

    if(n_steps>1) elbo_qgamma_old <- elbo_qgamma
    elbo_qgamma <- lgamma(sum(theta_star)) - sum(lgamma(theta_star)) + sum((theta_star-1)*expected_log_gamma)
    if(n_steps>1) elbo_qpi_old <- elbo_qpi
    elbo_qpi <- sum(lgamma(a_star+b_star)) - sum(lgamma(a_star)) - sum(lgamma(b_star)) + sum((a_star-1)*expected_log_pi) + sum((b_star-1)*expected_log_pi_complement)
    if(n_steps>1) elbo_qmusig_old <- elbo_qmusig
    elbo_qmusig <- 0.5*sum(log(v_star)) - 0.5*sum(log(2*pi)) + sum(alpha_star*log(beta_star)) - sum(lgamma(alpha_star)) - sum((alpha_star+3/2)*expected_log_sigma_squared) - sum(beta_star*expected_sigma2_inverse) - sum((v_star/2) * expected_mu2_sigma2) + sum(u_star*v_star*expected_mu_sigma2) - sum((((u_star**2)*v_star)/2) * expected_sigma2_inverse)
    if(n_steps>1) elbo_qz_old <- elbo_qz
    elbo_qz <- sum(expected_z * log(expected_z))
    if(n_steps>1) elbo_qd_old <- elbo_qd
    elbo_qd <- sum((expected_d * log(expected_d))[right_censored]) + sum(((1-expected_d)*log(1-expected_d))[right_censored])
    if(n_steps>1) elbo_qt_old <- elbo_qt
    elbo_qt <- elbo_qt_const - 0.5*log(1/(2*r_star)) - r_star*expected_t2 + q_star*expected_t - (q_star**2)/(4*r_star)
    elbo_qt_left <- left_censored*(elbo_qt - pnorm(rho, q_star/(2*r_star), sqrt(1/(2*r_star)), log.p = TRUE))
    elbo_qt_right <- right_censored*(elbo_qt  - pnorm(tau, q_star/(2*r_star), sqrt(1/(2*r_star)), log.p = TRUE, lower.tail = FALSE))
    elbo_qt <- sum(elbo_qt_left) + sum(elbo_qt_right)
    
    elbo <- elbo_gamma_prior + elbo_pi_prior + elbo_musig_prior + elbo_likelihood - elbo_qgamma - elbo_qpi - elbo_qmusig - elbo_qz - elbo_qd - elbo_qt
    elbos[n_steps] <- elbo
    
    #Update stopping condition difference
    if(n_steps>1) param_difference <- mean(abs(theta_star/sum(theta_star) - theta_star_old/sum(theta_star_old))) + mean(abs(a_star/(a_star+b_star) - a_star_old/(a_star_old + b_star_old))) + mean(abs(2*alpha_star - 2*alpha_star_old)) + mean(abs(u_star - u_star_old)) + mean(abs((beta_star*(v_star+1))/(alpha_star*v_star) - (beta_star_old*(v_star_old+1))/(alpha_star_old * v_star_old)))
    #print(param_difference)
    if(n_steps>1) elbo_difference <- elbo - elbos[n_steps-1]
    cat(n_steps)
    cat(" ")
    cat(param_difference)
    cat(" ")
    cat(elbo_difference)
    cat(" ")
    cat(elbo)
    cat("\n")
  }
  elbos <- elbos[1:n_steps]
  #Return outputs
  return(list(posterior.parameters = list(theta_star = theta_star, a_star = a_star, b_star = b_star, u_star = u_star, v_star = v_star, alpha_star = alpha_star, beta_star = beta_star, C_star = C_star, p_star = p_star, q_star = q_star, r_star = r_star, D_star = D_star, expected_t = expected_t), n_steps = n_steps, final_step_size = param_difference, elbo = elbos[1:n_steps], cond_list = cond_list))
  
}

