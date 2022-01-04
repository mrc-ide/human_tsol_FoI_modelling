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
  C2<- Burn(run2, burnin)
  
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
  
  if(number_datasets == 5){
    # Repplot auto/-correlation post processing #
    autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
    autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
    autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
    autocor.par4 <- cor(processed_chain[[1]][,4][-1], processed_chain[[1]][,4][-length(processed_chain[[1]][,4])])
    autocor.par5 <- cor(processed_chain[[1]][,5][-1], processed_chain[[1]][,5][-length(processed_chain[[1]][,5])])
    autocor.par6 <- cor(processed_chain[[1]][,6][-1], processed_chain[[1]][,6][-length(processed_chain[[1]][,6])])
    autocor.par7 <- cor(processed_chain[[1]][,7][-1], processed_chain[[1]][,7][-length(processed_chain[[1]][,7])])
    autocor.par8 <- cor(processed_chain[[1]][,8][-1], processed_chain[[1]][,8][-length(processed_chain[[1]][,8])])
    autocor.par9 <- cor(processed_chain[[1]][,9][-1], processed_chain[[1]][,9][-length(processed_chain[[1]][,9])])
    autocor.par10 <- cor(processed_chain[[1]][,10][-1], processed_chain[[1]][,10][-length(processed_chain[[1]][,10])])
    autocor.par11 <- cor(processed_chain[[1]][,11][-1], processed_chain[[1]][,11][-length(processed_chain[[1]][,11])])
    autocor.par12 <- cor(processed_chain[[1]][,12][-1], processed_chain[[1]][,12][-length(processed_chain[[1]][,12])])

    
    a <- signif(autocor.par1,3)
    b <- signif(autocor.par2,3)
    c <- signif(autocor.par3,3)
    d <- signif(autocor.par4,3)
    e <- signif(autocor.par5,3)
    f <- signif(autocor.par6,3)
    g <- signif(autocor.par7,3)
    h <- signif(autocor.par8,3)
    i <- signif(autocor.par9,3)
    j <- signif(autocor.par10,3)
    k <- signif(autocor.par11,3)
    l <- signif(autocor.par12,3)
    
    return(list(a, b, c, d, e, f, g, h, i, j, k, l))
  }
}


# function to obtain median and credible intervals for each parameter

obtain_parameter_values_func <- function(processed_chains, model, number_datasets){
  
  if(number_datasets == 1){
  lambda_simple<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025,0.5,0.975))
  se_simple<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025,0.5,0.975))
  sp_simple<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025,0.5,0.975))

  lambda.median<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.5))
  se.median<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.5))
  sp.median<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.5))

  lambda.credible<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025, 0.975))
  se.credible<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025, 0.975))
  sp.credible<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025, 0.975)) 
  }
  
  if(model == "simple" && number_datasets == 1){
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
  
  if(number_datasets == 5){
  
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
  }
  
  if(model == "simple" ){
  return(list(sp.median, se.median, lambda1.median, lambda2.median, lambda3.median, lambda4.median, lambda5.median,
              sp_simple, se_simple, lambda1_simple, lambda2_simple, lambda3_simple, lambda4_simple, lambda5_simple,
              sp.credible, se.credible, lambda1.credible, lambda2.credible, lambda3.credible, lambda4.credible, lambda5.credible))
  }
  
  if(model == "reversible" && number_datasets == 5){
    
    rho1_reversible <- quantile(c(processed_chains[[1]][,8], processed_chains[[2]][,8]), c(0.025,0.5,0.975))
    rho2_reversible <- quantile(c(processed_chains[[1]][,9], processed_chains[[2]][,9]), c(0.025,0.5,0.975))
    rho3_reversible <- quantile(c(processed_chains[[1]][,10], processed_chains[[2]][,10]), c(0.025,0.5,0.975))
    rho4_reversible <- quantile(c(processed_chains[[1]][,11], processed_chains[[2]][,11]), c(0.025,0.5,0.975))
    rho5_reversible <- quantile(c(processed_chains[[1]][,12], processed_chains[[2]][,12]), c(0.025,0.5,0.975))

    rho1.median <- quantile(c(processed_chains[[1]][,8], processed_chains[[2]][,8]), c(0.5))
    rho2.median <- quantile(c(processed_chains[[1]][,9], processed_chains[[2]][,9]), c(0.5))
    rho3.median <- quantile(c(processed_chains[[1]][,10], processed_chains[[2]][,10]), c(0.5))
    rho4.median <- quantile(c(processed_chains[[1]][,11], processed_chains[[2]][,11]), c(0.5))
    rho5.median <- quantile(c(processed_chains[[1]][,12], processed_chains[[2]][,12]), c(0.5))
    
    rho1.credible <- quantile(c(processed_chains[[1]][,8], processed_chains[[2]][,8]), c(0.025, 0.975))
    rho2.credible <- quantile(c(processed_chains[[1]][,9], processed_chains[[2]][,9]), c(0.025, 0.975))
    rho3.credible <- quantile(c(processed_chains[[1]][,10], processed_chains[[2]][,10]), c(0.025, 0.975))
    rho4.credible <- quantile(c(processed_chains[[1]][,11], processed_chains[[2]][,11]), c(0.025, 0.975))
    rho5.credible <- quantile(c(processed_chains[[1]][,12], processed_chains[[2]][,12]), c(0.025, 0.975))
    
    return(list(sp.median, se.median, lambda1.median, lambda2.median, lambda3.median, lambda4.median, lambda5.median, rho1.median,
                rho2.median, rho3.median, rho4.median, rho5.median,
                sp_simple, se_simple, lambda1_simple, lambda2_simple, lambda3_simple, lambda4_simple, lambda5_simple, rho1_reversible,
                rho2_reversible, rho3_reversible, rho4_reversible, rho5_reversible,
                sp.credible, se.credible, lambda1.credible, lambda2.credible, lambda3.credible, lambda4.credible, lambda5.credible,
                rho1.credible, rho2.credible, rho3.credible, rho4.credible, rho5.credible))
    
  }
}

# calculate model fit : deviance information criterion (DIC) statistic #
# for simple model
calculate_DIC_func1 <- function(chain, burnin, subsample, parameters, number_datasets){
  
  #remove burnin from loglik and logposterior chains
  loglike.chain.no.burnin_model1 <- chain$Likelihood[-(1:burnin)]
  logpost.chain.no.burnin_model1 <- chain$Posterior_Output[-(1:burnin)]
  
  ## thinning 
  loglike.chain.sub_model1 <- loglike.chain.no.burnin_model1[seq(1,length(loglike.chain.no.burnin_model1),subsample)]
  logpost.chain.sub_model1 <- logpost.chain.no.burnin_model1[seq(1,length(logpost.chain.no.burnin_model1),subsample)]
  
  ## computing goodness of fit (mean deviance) - deviance for given theta --> want to calculate mean of these 
  D_bar_model1 <- mean(-2 * loglike.chain.sub_model1) # mean deviance
  
  
  if(number_datasets == 1){
  # Posterior log likelihood
  modal_posterior_likelihood <- posterior_function_simple(data=data, c(parameters[[1]], 
                                                                       parameters[[2]], 
                                                                       parameters[[3]]))
  }
  
  if(number_datasets == 5){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_simple_multidata(data=data, c(parameters[[1]], parameters[[2]], 
                                                                         parameters[[3]], parameters[[4]], parameters[[5]],
                                                                         parameters[[6]], parameters[[7]]))
  }
  
 modal_posterior_deviance <- -2 * modal_posterior_likelihood
  
  # calculating the DIC 
  DIC <- 2 * D_bar_model1 - modal_posterior_deviance
  
  return(list(D_bar_model1, modal_posterior_likelihood, modal_posterior_deviance, DIC))
  
}

# for reversible model
calculate_DIC_func2 <- function(chain, burnin, subsample, parameters, number_datasets, simple_lambda_median){
  
  #remove burnin from loglik and logposterior chains
  loglike.chain.no.burnin_model1 <- chain$Likelihood[-(1:burnin)]
  logpost.chain.no.burnin_model1 <- chain$Posterior_Output[-(1:burnin)]
  
  ## thinning 
  loglike.chain.sub_model1 <- loglike.chain.no.burnin_model1[seq(1,length(loglike.chain.no.burnin_model1),subsample)]
  logpost.chain.sub_model1 <- logpost.chain.no.burnin_model1[seq(1,length(logpost.chain.no.burnin_model1),subsample)]
  
  ## computing goodness of fit (mean deviance) - deviance for given theta --> want to calculate mean of these 
  D_bar_model1 <- mean(-2 * loglike.chain.sub_model1) # mean deviance
  
  
 if(number_datasets == 1){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_reversible(data=data, c(parameters[[1]], 
                                                                             parameters[[2]], 
                                                                             parameters[[3]],
                                                                             parameters[[4]]), simple_lambda_median)
  }
  
  if(number_datasets == 5){
    # Posterior log likelihood
    modal_posterior_likelihood <- posterior_function_reversible_multidata(data=data, c(parameters[[1]], parameters[[2]], parameters[[3]], 
                                                                                       parameters[[4]], parameters[[5]], parameters[[6]], 
                                                                                       parameters[[7]], parameters[[8]], parameters[[9]], 
                                                                                       parameters[[10]], parameters[[11]], parameters[[12]]),
                                                                          simple_lambda_median)
  }
  
  
  modal_posterior_deviance <- -2 * modal_posterior_likelihood
  
  # calculating the DIC 
  DIC <- 2 * D_bar_model1 - modal_posterior_deviance
  
  return(list(D_bar_model1, modal_posterior_likelihood, modal_posterior_deviance, DIC))
  
}
