#=============================================================================================================================#
#                                      Process MCMC outputs functions                                                         #
#=============================================================================================================================#

#===================================================#
#  Simple FoI moodel: single dataset                #

# function to plot chains
plot_chains <- function(run1, run2){
  par(mfrow=c(ncol(run1),1))
  
  for(i in 1:ncol(run1)){
    plot(run1[,i], t='l', col='deeppink',
         ylim=c(min(c(run1[,i], run2[,i])),max(c(run1[,i], run2[,i]))),
         xlab='', ylab=paste('Parameter', i, sep=' '))
    lines(run2[,i], col='dodgerblue')
  }
  
}

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

# function to calculate autocorrelation for each parameter (simple model - single dataset)
determine_autocorrelation_func1 <- function(processed_chain){
  
  # Repplot auto/-correlation post processing #
  autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
  autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
  autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
  
  
  a <- signif(autocor.par1,3)
  b <- signif(autocor.par2,3)
  c <- signif(autocor.par3,3)
  
  return(list(a, b, c))
}

# function to calculate autocorrelation for each parameter (reversible model - single dataset)
determine_autocorrelation_func2 <- function(processed_chain){
  
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


# function to obtain median and credible intervals for each parameter

obtain_parameter_values_func <- function(processed_chains, model){
  
  lambda_simple<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025,0.5,0.975))
  se_simple<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025,0.5,0.975))
  sp_simple<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025,0.5,0.975))

  lambda.median<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.5))
  se.median<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.5))
  sp.median<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.5))

  lambda.credible<-quantile(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), c(0.025, 0.975))
  se.credible<-quantile(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), c(0.025, 0.975))
  sp.credible<-quantile(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), c(0.025, 0.975))
  
  if(model == "reversible"){
    rho_reversible <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025,0.5,0.975))
    rho.median <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.5))
    rho.credible <- quantile(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), c(0.025, 0.975))
    
  }
  
  if(model == "simple"){
  return(list(lambda.median, se.median, sp.median,
            lambda_simple, se_simple, sp_simple,
            lambda.credible, se.credible, sp.credible))
         }
  
  if(model == "reversible"){
    return(list(lambda.median, se.median, sp.median, rho.median,
                lambda_simple, se_simple, sp_simple, rho_reversible,
                lambda.credible, se.credible, sp.credible, rho.credible))
           }
  }

# calculate model fit : deviance information criterion (DIC) statistic #

calculate_DIC_func <- function(chain, burnin, subsample, parameters, model){
  
  #remove burnin from loglik and logposterior chains
  loglike.chain.no.burnin_model1 <- chain$Likelihood[-(1:burnin)]
  logpost.chain.no.burnin_model1 <- chain$Posterior_Output[-(1:burnin)]
  
  ## thinning 
  loglike.chain.sub_model1 <- loglike.chain.no.burnin_model1[seq(1,length(loglike.chain.no.burnin_model1),subsample)]
  logpost.chain.sub_model1 <- logpost.chain.no.burnin_model1[seq(1,length(logpost.chain.no.burnin_model1),subsample)]
  
  ## computing goodness of fit (mean deviance) - deviance for given theta --> want to calculate mean of these 
  D_bar_model1 <- mean(-2 * loglike.chain.sub_model1) # mean deviance
  
  
  if(model == "simple"){
  # Posterior log likelihood
  modal_posterior_likelihood <- posterior_function_simple(data=data, c(parameters[[1]], 
                                                                       parameters[[2]], 
                                                                       parameters[[3]]))
  }
  
  if(model == "reversible"){
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
