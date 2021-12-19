#=============================================================================================================================#
#                                       Plot MCMC outputs                                                                     #
#=============================================================================================================================#



#### chains plot #####

chains_plot_func1 <- function(inits1, chain1, chain2) {
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
    
  }
  
}

#### posterior histogram plots #####

histogram_plot_func1 <- function(inits1, burnin, niter, chain1, chain2) {
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
      
    }
  }

##### plot log likelihood #####

loglikchains_plot_func1 <- function(chain1, chain2){
  par(mfrow=c(1,1))
  plot(simple_out_chain1$Likelihood_Output, t='l', ylab='simple LL')
  lines(simple_out_chain2$Likelihood_Output, col='red')
}

#### plot posterior distributions ###

plot_posterior_distrib_func1 <- function(processed_chains, number_pars_to_plot){
  
  par(mfrow=c(1,number_pars_to_plot))
  
  hist(c(PC_simple[[1]][,1], PC_simple[[2]][,1]), breaks=30, xlab='Lambda')  # Parameter 1 - lambda 
  hist(c(PC_simple[[1]][,2], PC_simple[[2]][,2]), breaks=30, xlab='sensitivity')      # Parameter 2 - rho
  hist(c(PC_simple[[1]][,3], PC_simple[[2]][,3]), breaks=30, xlab='specificity')      # Parameter 2 - rho
}

