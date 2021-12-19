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

# function to calculate autocorrelation for each parameter
determine_autocorrelation_func <- function(processed_chain){
  
  # Repplot auto/-correlation post processing #
  autocor.par1 <- cor(processed_chain[[1]][,1][-1], processed_chain[[1]][,1][-length(processed_chain[[1]][,1])])
  autocor.par2 <- cor(processed_chain[[1]][,2][-1], processed_chain[[1]][,2][-length(processed_chain[[1]][,2])])
  autocor.par3 <- cor(processed_chain[[1]][,3][-1], processed_chain[[1]][,3][-length(processed_chain[[1]][,3])])
  
  
  a <- signif(autocor.par1,3)
  b <- signif(autocor.par2,3)
  c <- signif(autocor.par3,3)
  
  return(list(a, b, c))
}

# function to obtain median and credible intervals for each parameter

obtain_parameter_values_func <- function(processed_chains){
  
  lambda_simple<-quantile(c(PC_simple[[1]][,1], PC_simple[[2]][,1]), c(0.025,0.5,0.975))
  se_simple<-quantile(c(PC_simple[[1]][,2], PC_simple[[2]][,2]), c(0.025,0.5,0.975))
  sp_simple<-quantile(c(PC_simple[[1]][,3], PC_simple[[2]][,3]), c(0.025,0.5,0.975))

  lambda.median<-quantile(c(PC_simple[[1]][,1], PC_simple[[2]][,1]), c(0.5))
  se.median<-quantile(c(PC_simple[[1]][,2], PC_simple[[2]][,2]), c(0.5))
  sp.median<-quantile(c(PC_simple[[1]][,3], PC_simple[[2]][,3]), c(0.5))

  lambda.credible<-quantile(c(PC_simple[[1]][,1], PC_simple[[2]][,1]), c(0.025, 0.975))
  se.credible<-quantile(c(PC_simple[[1]][,2], PC_simple[[2]][,2]), c(0.025, 0.975))
  sp.credible<-quantile(c(PC_simple[[1]][,3], PC_simple[[2]][,3]), c(0.025, 0.975))

return(list(lambda_simple, se_simple, sp_simple,
            lambda.median, se.median, sp.median,
            lambda.credible, se.credible, sp.credible))
}

# calculate model fit : deviance information criterion (DIC) statistic #

calculate_DIC_func <- function(chain, burnin, subsample){
  
  #remove burnin from loglik and logposterior chains
  loglike.chain.no.burnin_model1 <- chain$Likelihood[-(1:burnin)]
  logpost.chain.no.burnin_model1 <- chain$Posterior_Output[-(1:burnin)]
  
  ## thinning 
  loglike.chain.sub_model1 <- loglike.chain.no.burnin_model1[seq(1,length(loglike.chain.no.burnin_model1),subsample)]
  logpost.chain.sub_model1 <- logpost.chain.no.burnin_model1[seq(1,length(logpost.chain.no.burnin_model1),subsample)]
  
  ## computing goodness of fit (mean deviance) - deviance for given theta --> want to calculate mean of these 
  D_bar_model1 <- mean(-2 * loglike.chain.sub_model1) # mean deviance
  
  # Posterior log likelihood
  modal_posterior_likelihood <- posterior_function_simple(data=data, c(simple_model_parameters[[1]], 
                                                                       simple_model_parameters[[2]], 
                                                                       simple_model_parameters[[3]]))
  modal_posterior_deviance <- -2 * modal_posterior_likelihood
  
  # calculating the DIC 
  DIC <- 2 * D_bar_model1 - modal_posterior_deviance
  
  return(list(D_bar_model1, modal_posterior_likelihood, modal_posterior_deviance, DIC))
  
  }
