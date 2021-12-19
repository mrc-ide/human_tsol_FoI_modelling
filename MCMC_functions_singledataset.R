#=============================================================================================================================#
#                           MCMC function (single dataset)                                                                    #
#=============================================================================================================================#



#============================================================#
#       Simple FoI model (single dataset)                    #
#============================================================#

# Prior Function
prior <- function(par) {
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  lambda = lambda
  
  lambda_prior = dunif(lambda, min = 0.000001, max = 12, log = T)
  
  se = se
  se_prior = dbeta(se, alpha_se, beta_se, log = T)
  
  sp = sp
  sp_prior = dbeta(sp, alpha_sp, beta_sp, log = T)
  
  return(sum(c(lambda_prior, se_prior, sp_prior)))
  
}


# Likelihood Function
loglike_simple <- function(data, par){
  predicted_seroprevalence = predicted_prev_func(age = data$age, par)
  sum(dbinom(data$pos_pigs, data$n_pigs, predicted_seroprevalence, log=T))
}



#loglike_simple(data,c(0.1, 0.9, 0.9))

# Posterior Function
posterior_function_simple <- function(data, par){
  loglike_simple(data, par) + prior(par)
}

#posterior_function_simple(data, c(0.1, 0.945, 0.99))


proposal_function_simple <- function(par, cov) {
  #lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  ## draw propopsals all from a multivariate normal
  repeat {
    proposed <- mvrnorm(1, par, cov)
    if(all(proposed[2]>0 & all(proposed[2]<1))& all(proposed[3]>0 & all(proposed[3]<1)) & all(proposed[1]>0 & all(proposed[1]<12))){break}
  }
  
  
  return(proposed)
  
}  


# Run MCMC Function 
MCMC_simple_model <- function(inits,  number_of_iterations, cov) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol=length(inits)) # ncol is number of parameters
  MCMC_output[1,] <- inits
  Acceptances <- vector(length = number_of_iterations + 1)
  LogLikelihood_storage <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    proposed_parameter_value <- proposal_function_simple(MCMC_output[i,], cov)  #new proposed paramater value(s) with a given s.d. (step size)
    
    current_likelihood <- loglike_simple(data, MCMC_output[i,]) # likelihood 
    
    current_posterior <- posterior_function_simple(data, MCMC_output[i,]) # current posterior likelihood from MCMC
    
    proposed_posterior <- posterior_function_simple(data, proposed_parameter_value) # proposed posterior likelihood with new proposed par value
    
    likelihood_ratio = exp(proposed_posterior - current_posterior);
    
    if(i %% (number_of_iterations/20) == 0){
      message(round(i/number_of_iterations*100), ' % completed')
    }
    
    if(runif(1) < likelihood_ratio) {
      
      MCMC_output[i + 1,] <- proposed_parameter_value
      Acceptances[i] <- 1
      LogLikelihood_storage[i + 1] <- current_likelihood
      Logposterior_storage[i + 1] <- proposed_posterior
      
    } else{
      
      MCMC_output[i + 1,] <- MCMC_output[i,]
      Acceptances[i] <- 0
      LogLikelihood_storage[i + 1] <- current_likelihood
      Logposterior_storage[i + 1] <- current_posterior
      
    }
    
  } # likelihood ratio comparison step (exponentiated because on log scale) - accept (1) if posterior proposed improved LL on posterior current (runif- generates random variable between 0 -1)
  
  list <- list()
  list[["MCMC_Output"]] <- MCMC_output
  list[["Acceptances"]] <- Acceptances
  list[["Likelihood_Output"]] <- LogLikelihood_storage
  list[["Posterior_Output"]] <- Logposterior_storage
  return(list)
  
}