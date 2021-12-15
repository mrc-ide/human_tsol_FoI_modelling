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

histogram_plot_func1 <-
  function(inits1, burnin, niter, chain1, chain2) {
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
        c(
          simple_out_chain1$MCMC_Output[burnin:niter, i],
          simple_out_chain2$MCMC_Output[burnin:niter, i]
        ),
        xlab = ylab,
        main = ""
      )
      
    }
  }
