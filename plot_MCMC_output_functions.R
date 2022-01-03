#=============================================================================================================================#
#                                       Plot MCMC outputs                                                                     #
#=============================================================================================================================#

#======================================#
# simple model - single dataset        #

#### chains plot #####

chains_plot_func <- function(inits1, chain1, chain2, model) {
  
  if(model == "simple"){
  par(mfrow = (c(1, length(inits1))))
  for (i in 1:length(inits1)) {
    if (i == 1) {
      ylab = "lambda"
    } else if (i == 2) {
      ylab = "se"
    } else {
      ylab = "sp"
    }
    plot(
      simple_out_chain1$MCMC_Output[, i],
      type = "l",
      ylab = ylab,
      xlab = "iter"
    )
    lines(simple_out_chain2$MCMC_Output[, i], col = "red")
    
  }}
  
  if(model == "reversible"){
    par(mfrow=(c(1,length(inits1))))
    for (i in 1:length(inits1)) {
      if (i==1) {
        ylab="lambda"
      } else if (i==2) {
        ylab="se"
      } else if (i==3) {
        ylab="sp"
      } else {
        ylab="rho"
      }
      plot(reversible_out_chain1$MCMC_Output[,i], type = "l", ylab=ylab, xlab="iter", col = "black")
      lines(reversible_out_chain2$MCMC_Output[,i], col="red")
    }}
  
}

#### posterior histogram plots #####

histogram_plot_func <- function(inits1, burnin, niter, chain1, chain2, model) {
    
  if(model == "simple"){
  par(mfrow = (c(1, length(inits1))))
    for (i in 1:length(inits1)) {
      if (i == 1) {
        ylab = "lambda"
      } else if (i == 2) {
        ylab = "se"
      } else {
        ylab = "sp"
      }
      
      hist(
        c(simple_out_chain1$MCMC_Output[burnin:niter, i], simple_out_chain2$MCMC_Output[burnin:niter, i]),
        xlab = ylab,
        main = ""
      )
      
    }}
  
  if(model == "reversible"){
    par(mfrow=(c(1,length(inits1))))
    for (i in 1:length(inits1)) {
      if (i==1) {
        ylab="lambda"
      } else if (i==2) {
        ylab="se"
      } else if (i==3) {
        ylab="sp"
      } else {
        ylab="rho"
      }
      
      hist(c(reversible_out_chain1$MCMC_Output[burnin:niter,i],reversible_out_chain2$MCMC_Output[burnin:niter,i]), 
           xlab = ylab, main="")
      
    }}
  
  }

##### plot log likelihood #####

loglikchains_plot_func <- function(chain1, chain2){
  par(mfrow=c(1,1))
  plot(chain1$Likelihood_Output, t='l', ylab='Loglikelihood', xlab='iteration')
  lines(chain2$Likelihood_Output, col='red')
}

#### plot posterior distributions ###

plot_posterior_distrib_func <- function(processed_chains, model){
  
  if(model == "simple"){
  par(mfrow=c(1,3))
  
  hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='Lambda')  # Parameter 1 - lambda 
  hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity')      # Parameter 2 - se
  hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='specificity')      # Parameter 3 - sp
  }
  
  if(model == "reversible"){
    par(mfrow=c(1,4))
    hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='Lambda')  # Parameter 1 - lambda 
    hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity')      # Parameter 2 - se
    hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='specificity')      # Parameter 3 - sp
    hist(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), breaks=30, xlab='rho')      # Parameter 4 - rho
  }
  
  
}


