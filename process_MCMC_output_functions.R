#=============================================================================================================================#
#                                      Process MCMC outputs functions                                                         #
#=============================================================================================================================#

# function to remove burnin
Burn <- function(chains, burnin){
  chains[-(1:burnin),]
}

# function to sub-sample
Downsample <- function(chains, sample){
  chains[seq(1, nrow(chains), sample),]
}

# function calling burnin and sub-sampling
Process_chains <- function(run1, run2, burnin, sample){
  C1 <- Burn(run1, burnin)
  C2<-Burn(run2, burnin)
  
  C1 <- Downsample(C1, sample)
  C2 <- Downsample(C2, sample)
  
  
  return(list(C1, C2))
}

# function to calculate autocorrelation for each parameter (simple model)
determine_autocorrelation_func1 <- function(processed_chain, number_datasets){
  
  if(number_datasets == 1){
  # Repplot auto/-correlation post processing #
  autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
  autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
  autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
  
  
  a <- signif(autocor.par1,3)
  b <- signif(autocor.par2,3)
  c <- signif(autocor.par3,3)
  
  return(list(a, b, c))
  }
  
  if(number_datasets == 5){
    # Repplot auto/-correlation post processing #
    autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
    autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
    autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
    autocor.par4 <- cor(processed_chain[[1]][,4][-1], processed_chain[[1]][,4][-length(processed_chain[[1]][,4])])
    autocor.par5 <- cor(processed_chain[[1]][,5][-1], processed_chain[[1]][,5][-length(processed_chain[[1]][,5])])
    autocor.par6 <- cor(processed_chain[[1]][,6][-1], processed_chain[[1]][,6][-length(processed_chain[[1]][,6])])
    autocor.par7 <- cor(processed_chain[[1]][,7][-1], processed_chain[[1]][,7][-length(processed_chain[[1]][,7])])
    
    a <- signif(autocor.par1,3)
    b <- signif(autocor.par2,3)
    c <- signif(autocor.par3,3)
    d <- signif(autocor.par4,3)
    e <- signif(autocor.par5,3)
    f <- signif(autocor.par6,3)
    g <- signif(autocor.par7,3)
    
    return(list(a, b, c, d, e, f, g))
  }
  
}

# function to calculate autocorrelation for each parameter (reversible model)
determine_autocorrelation_func2 <- function(processed_chain, number_datasets){
  
  if(number_datasets == 1){
  # Repplot auto/-correlation post processing #
  autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
  autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
  autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
  autocor.par4 <- cor(processed_chain[[1]][,4][-1], processed_chain[[1]][,4][-length(processed_chain[[1]][,4])])
  
  a <- signif(autocor.par1,3)
  b <- signif(autocor.par2,3)
  c <- signif(autocor.par3,3)
  d <- signif(autocor.par3,4)
  
  return(list(a, b, c, d))
  }
}


# function to obtain median and credible intervals for each parameter

obtain_parameter_values_func <- function(processed_chains, model, number_datasets){
  
  if(model == "simple" && number_datasets == 1){
  lambda_simple<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025,0.5,0.975))
  se_simple<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025,0.5,0.975))
  sp_simple<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025,0.5,0.975))

  lambda.median<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.5))
  se.median<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.5))
  sp.median<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.5))

  lambda.credible<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025, 0.975))
  se.credible<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025, 0.975))
  sp.credible<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025, 0.975))  
  
  return(list(lambda.median, se.median, sp.median,
            lambda_simple, se_simple, sp_simple,
            lambda.credible, se.credible, sp.credible))
  }
  
  if(model == "reversible" && number_datasets == 1){
    rho_reversible <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025,0.5,0.975))
    rho.median <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.5))
    rho.credible <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025, 0.975))
    
    return(list(lambda.median, se.median, sp.median, rho.median,
                lambda_simple, se_simple, sp_simple, rho_reversible,
                lambda.credible, se.credible, sp.credible, rho.credible))
  }
  
  if(model == "simple" && number_datasets == 5){
  
  sp_simple<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025,0.5,0.975))
  se_simple<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025,0.5,0.975))
  lambda1_simple<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025,0.5,0.975))
  lambda2_simple<-quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025,0.5,0.975))
  lambda3_simple<-quantile(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), c(0.025,0.5,0.975))
  lambda4_simple<-quantile(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), c(0.025,0.5,0.975))
  lambda5_simple<-quantile(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), c(0.025,0.5,0.975))
  
  sp.median <-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.5))
  se.median <-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.5))
  lambda1.median <-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.5))
  lambda2.median <-quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.5))
  lambda3.median <-quantile(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), c(0.5))
  lambda4.median <-quantile(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), c(0.5))
  lambda5.median <-quantile(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), c(0.5))
  
  sp.credible <-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025,0.975))
  se.credible <-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025,0.975))
  lambda1.credible <-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025,0.975))
  lambda2.credible <-quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025,0.975))
  lambda3.credible <-quantile(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), c(0.025,0.975))
  lambda4.credible <-quantile(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), c(0.025,0.975))
  lambda5.credible <-quantile(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), c(0.025,0.975))
  
  return(list(sp.median, se.median, lambda1.median, lambda2.median, lambda3.median, lambda4.median, lambda5.median,
              sp_simple, se_simple, lambda1_simple, lambda2_simple, lambda3_simple, lambda4_simple, lambda5_simple,
              sp.credible, se.credible, lambda1.credible, lambda2.credible, lambda3.credible, lambda4.credible, lambda5.credible))
  
  }
  
  if(model == "reversible" && number_datasets == 5){
    
    rho_reversible <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025,0.5,0.975))
    rho.median <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.5))
    rho.credible <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025, 0.975))
    
  }
  
}

# calculate model fit : deviance information criterion (DIC) statistic #

calculate_DIC_func <- function(chain, burnin, subsample, parameters, model, number_datasets){
  
  #remove burnin from loglik and logposterior chains
  loglike.chain.no.burnin_model1 <- chain$Likelihood[-(1:burnin)]
  logpost.chain.no.burnin_model1 <- chain$Posterior_Output[-(1:burnin)]
  
  ## thinning 
  loglike.chain.sub_model1 <- loglike.chain.no.burnin_model1[seq(1,length(loglike.chain.no.burnin_model1),subsample)]
  logpost.chain.sub_model1 <- logpost.chain.no.burnin_model1[seq(1,length(logpost.chain.no.burnin_model1),subsample)]
  
  ## computing goodness of fit (mean deviance) - deviance for given theta --> want to calculate mean of these 
  D_bar_model1 <- mean(-2 * loglike.chain.sub_model1) # mean deviance
  
  
  if(model == "simple" && number_datasets == 1){
  # Posterior log likelihood
  modal_posterior_likelihood <- posterior_function_simple(data=data, c(parameters[[1]], 
                                                                       parameters[[2]], 
                                                                       parameters[[3]]))
  }
  
  if(model == "simple" && number_datasets == 5){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_simple_multidata(data=data, c(parameters[[1]], parameters[[2]], 
                                                                         parameters[[3]], parameters[[4]], parameters[[5]],
                                                                         parameters[[6]], parameters[[7]]))
  }
  
  if(model == "reversible" && number_datasets == 1){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_reversible2(data=data, c(parameters[[1]], 
                                                                         parameters[[2]], 
                                                                         parameters[[3]],
                                                                         parameters[[4]]))
  }
  
  modal_posterior_deviance <- -2 * modal_posterior_likelihood
  
  # calculating the DIC 
  DIC <- 2 * D_bar_model1 - modal_posterior_deviance
  
  return(list(D_bar_model1, modal_posterior_likelihood, modal_posterior_deviance, DIC))
  
  }
