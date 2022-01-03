#=============================================================================================================================#
#                           MCMC function (single dataset)                                                                    #
#=============================================================================================================================#

#============================================================#
#       Simple FoI model (single dataset)                    #
#============================================================#

# # prior distributions for each parameter #
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


# likelihood calculation given each parameter set & data
loglike_simple <- function(data, par){
  predicted_seroprevalence = predicted_prev_func(age = data$age, par)
  sum(dbinom(data$pos, data$n, predicted_seroprevalence, log=T))
}


# Posterior Function
posterior_function_simple <- function(data, par){
  loglike_simple(data, par) + prior(par)
}

# proposal function for MCMC
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
MCMC_simple_model <- function(inits,  number_of_iterations, cov, fitting) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol=length(inits)) # ncol is number of parameters
  MCMC_output[1,] <- inits
  Acceptances <- vector(length = number_of_iterations + 1)
  LogLikelihood_storage <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    if(fitting == "single dataset"){
    proposed_parameter_value <- proposal_function_simple(MCMC_output[i,], cov)  #new proposed paramater value(s) with a given s.d. (step size)
    }
    
    if(fitting == "multiple datasets"){
      proposed_parameter_value <- proposal_function_simple_multidata(MCMC_output[i,], cov)  #new proposed paramater value(s) with a given s.d. (step size)
    }
    
    if(fitting == "single dataset"){
    current_likelihood <- loglike_simple(data, MCMC_output[i,]) # likelihood 
    }
    
    if(fitting == "multiple datasets"){
      current_likelihood <- loglike_simple_multidata(data, MCMC_output[i,]) # likelihood 
    }
    
    if(fitting == "single dataset"){
    current_posterior <- posterior_function_simple(data, MCMC_output[i,]) # current posterior likelihood from MCMC
    }
    
    if(fitting == "multiple datasets"){
      current_posterior <- posterior_function_simple_multidata(data, MCMC_output[i,]) # current posterior likelihood from MCMC
    }
    
    if(fitting == "single dataset"){
    proposed_posterior <- posterior_function_simple(data, proposed_parameter_value) # proposed posterior likelihood with new proposed par value
    }
    
    if(fitting == "multiple datasets"){
      proposed_posterior <- posterior_function_simple_multidata(data, proposed_parameter_value) # proposed posterior likelihood with new proposed par value
    }
    
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

#===========================================================#
#       Reversible FoI model (single dataset)               #
#===========================================================#

# prior distributions for each parameter #
# specify lambda median from simple model to inform the lognromal prior
prior_function_reversible <- function(par, simple_lambda_median) {
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  lambda = lambda
  lambda_prior = dlnorm(lambda, log(simple_lambda_median), 1, log=TRUE)
  
  se = se
  se_prior = dbeta(se, alpha_se, beta_se, log = T)
  
  sp = sp
  sp_prior = dbeta(sp, alpha_sp, beta_sp, log = T)
  
  rho = rho
  rho_prior = dunif(rho, min = 0.000001, max = 12, log= T)
  
  return(sum(c(lambda_prior, se_prior, sp_prior, rho_prior)))
  
}

# when calculating postrior after MCMC fitting (i.e. replace simple median lambda with new fitted lambda)
prior_function_reversible2 <- function(par) {
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  lambda = lambda
  lambda_prior = dlnorm(lambda, log(lambda), 1, log=TRUE)
  
  se = se
  se_prior = dbeta(se, alpha_se, beta_se, log = T)
  
  sp = sp
  sp_prior = dbeta(sp, alpha_sp, beta_sp, log = T)
  
  rho = rho
  rho_prior = dunif(rho, min = 0.000001, max = 12, log= T)
  
  return(sum(c(lambda_prior, se_prior, sp_prior, rho_prior)))
  
}

# likelihood calculation given each parameter set & data
log_lik_func_reversible <- function(data, par){
  
  predicted_seroprevalence = predicted_prev_reversible_func(data$age, par)
  
  sum(dbinom(data$pos, data$n, predicted_seroprevalence, log=T))
}

# posterior calculation (for MCMC fitting)
posterior_function_reversible<- function(data, par, simple_lambda_median){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  log_lik_func_reversible(data, par) + prior_function_reversible(par, simple_lambda_median)
}

# posterior calculation (post mCMC fitting when have fitted lambda)
posterior_function_reversible2 <- function(data, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  log_lik_func_reversible(data, par) + prior_function_reversible2(par)
}

# proposal function for MCMC
proposal_function_reversible <- function(par, cov) {
  
  ## draw propopsals all from a multivariate normal
  repeat {
    proposed <- mvrnorm(1, par, cov)
    if(all(proposed[2]>0 & all(proposed[2]<1)) & all(proposed[3]>0 & all(proposed[3]<1)) & all(proposed[1]>0 & all(proposed[1]<12)) & all(proposed[4]>0 & all(proposed[4]<12))){break}
  }
  
  return(proposed)
  
}  

MCMC_reversible_model <- function(inits,  number_of_iterations, cov, simple_lambda_median) {
  
  # Storage for Output
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol=length(inits)) # ncol is number of parameters
  MCMC_output[1,] <- inits
  Acceptances <- vector(length = number_of_iterations + 1)
  LogLikelihood_storage <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    proposed_parameter_value <- proposal_function_reversible(MCMC_output[i,], cov)  #new proposed paramater value(s) with a given s.d. (step size)
    
    current_likelihood <- log_lik_func_reversible(data, MCMC_output[i,]) # likelihood 
    
    current_posterior <- posterior_function_reversible(data, MCMC_output[i,], simple_lambda_median) # current posterior likelihood from MCMC
    
    proposed_posterior <- posterior_function_reversible(data, proposed_parameter_value, simple_lambda_median) # proposed posterior likelihood with new proposed par value
    
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

#============================================================#
#       Simple FoI model (multiple dataset)                  #
#============================================================#

prior_simple_multidata <- function(par) {
  # Sp and Se are the first to parameters inthe vector of parameters
  sp <- par[1]
  se <- par[2]
  
  # Site-specific lamdas are the remianing parameters in the vector of parameters
  lambda_par <- par[3:length(par)]
  
  lambda_prior = sum(dunif(lambda_par, min = 0.00001, max = 12, log = T))
  
  # diagnostic priors
  se_prior = dbeta(se,alpha_se, beta_se, log = T) # 88.89*0.185, 11.11*0.335
  sp_prior = dbeta(sp,alpha_sp, beta_sp, log = T)
  
  return(sum(c(lambda_prior, se_prior, sp_prior)))
}

# likelihood calculation given each parameter set & data
loglike_simple_multidata <- function(data, par){
  predicted_seroprevalence = predicted_prev_func_multidataset(age = data$age, par)
  sum(dbinom(data$pos, data$n, predicted_seroprevalence, log=T))
}

# Posterior Function (multidata)
posterior_function_simple_multidata <- function(data, par){
  loglike_simple_multidata(data, par) + prior_simple_multidata(par)
}

# proposal func for multiple data
proposal_function_simple_multidata <- function(par, cov) {
  
  ## draw propopsals all from a multivariate normal & re-draw if any <=0 or > 12 (for upper limit of FoI)
  repeat {
    proposed <- mvrnorm(1, par, cov)
    if(all(proposed>0 & all(proposed[1:length(par)]<12))){break} 
  }
  
  return(proposed)
  
}  
