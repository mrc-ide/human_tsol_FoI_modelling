#=============================================================================================================================#
#                                       Plot MCMC outputs                                                                     #
#=============================================================================================================================#

#======================================#
#       Single datasets                #

#### chains plot (single dataset) - not processed #####

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


#### posterior histogram plots - not processed #####

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

#### function to plot chains - post-processing ####
plot_chains <- function(run1, run2){
  par(mfrow=c(ncol(run1),1))
  
  for(i in 1:ncol(run1)){
    plot(run1[,i], t='l', col='deeppink',
         ylim=c(min(c(run1[,i], run2[,i])),max(c(run1[,i], run2[,i]))),
         xlab='', ylab=paste('Parameter', i, sep=' '))
    lines(run2[,i], col='dodgerblue')
  }
  
}

plot_chains_multidatasets <- function(run1, run2, inits, number_datasets){
  par(mfrow=c(2,length(inits)/2))
  
  if(number_datasets == 5){
    for(i in 1:ncol(run1)){
      if (i == 1) {
        ylab = "sp"
      } else if (i == 2) {
        ylab = "se"
      }  else if (i ==3) {
        ylab="lambda (site 1)"
      } else if (i==4) {
        ylab="lambda (site 2)"
      } else if (i==5) {
        ylab="lambda (site 3)"
      } else if (i==6) {
        ylab="lambda (site 4)"
      } else {
        ylab="lambda (site 5)"
      }
      plot(run1[,i], t='l', col='deeppink',
           ylim=c(min(c(run1[,i], run2[,i])),max(c(run1[,i], run2[,i]))),
           xlab='', ylab = ylab)
      lines(run2[,i], col='dodgerblue')
    }}
  
}

#### plot posterior distributions (after processing) ###

plot_posterior_distrib_func <- function(processed_chains, model, number_datasets){
  
  if(model == "simple" && number_datasets == 1){
  par(mfrow=c(1,3))
  
  hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='Lambda', main="")  # Parameter 1 - lambda 
  hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity', main="")      # Parameter 2 - se
  hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='specificity', main="")      # Parameter 3 - sp
  }
  
  if(model == "simple" && number_datasets == 5){
    par(mfrow=c(2,ncol(PC_simple[[1]])/2))
    
    hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='specifcity', main="")  # Parameter 1  
    hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity', main="")      # Parameter 2 
    hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='lambda (site 1)', main="")      # Parameter 3 
    hist(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), breaks=30, xlab='lambda (site 2)', main="")  
    hist(c(processed_chains[[1]][,5], processed_chains[[2]][,5]), breaks=30, xlab='lambda (site 3)', main="")     
    hist(c(processed_chains[[1]][,6], processed_chains[[2]][,6]), breaks=30, xlab='lambda (site 4)', main="")      
    hist(c(processed_chains[[1]][,7], processed_chains[[2]][,7]), breaks=30, xlab='lambda (site 5)', main="")      
  }
  
  if(model == "reversible"){
    par(mfrow=c(1,4))
    hist(c(processed_chains[[1]][,1], processed_chains[[2]][,1]), breaks=30, xlab='Lambda', main="")  # Parameter 1 - lambda 
    hist(c(processed_chains[[1]][,2], processed_chains[[2]][,2]), breaks=30, xlab='sensitivity', main="")      # Parameter 2 - se
    hist(c(processed_chains[[1]][,3], processed_chains[[2]][,3]), breaks=30, xlab='specificity', main="")      # Parameter 3 - sp
    hist(c(processed_chains[[1]][,4], processed_chains[[2]][,4]), breaks=30, xlab='rho', main="")      # Parameter 4 - rho
  }
  
  
}

#================================================#
#           Multiple datasets                   #

#### chains plot (multiple datasets) #####

chains_plot_multidatasets_func <- function(inits1, chain1, chain2, model, number_datasets) {
  
  if(model == "simple" && number_datasets == 5){
    par(mfrow = (c(2, length(inits1)/2)))
    for (i in 1:length(inits1)) {
      if (i == 1) {
        ylab = "sp"
      } else if (i == 2) {
        ylab = "se"
      }  else if (i ==3) {
        ylab="lambda (site 1)"
    } else if (i==4) {
      ylab="lambda (site 2)"
    } else if (i==5) {
      ylab="lambda (site 3)"
    } else if (i==6) {
      ylab="lambda (site 4)"
    } else {
      ylab="lambda (site 5)"
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
    par(mfrow=(c(2,length(inits1)/2)))
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

histogram_plot_multidatasets_func <- function(inits1, burnin, niter, chain1, chain2, model, number_datasets) {
  
  if(model == "simple" && number_datasets == 5){
    par(mfrow = (c(2, length(inits1)/2)))
    for (i in 1:length(inits1)) {
      if (i == 1) {
        ylab = "sp"
      } else if (i == 2) {
        ylab = "se"
      }  else if (i ==3) {
        ylab="lambda (site 1)"
      } else if (i==4) {
        ylab="lambda (site 2)"
      } else if (i==5) {
        ylab="lambda (site 3)"
      } else if (i==6) {
        ylab="lambda (site 4)"
      } else {
        ylab="lambda (site 5)"
      }
      hist(
        c(simple_out_chain1$MCMC_Output[burnin:niter, i], simple_out_chain2$MCMC_Output[burnin:niter, i]),
        xlab = ylab,
        main = ""
      )
      
    }}
  
  if(model == "reversible"){
    par(mfrow=(c(2,length(inits1)/2)))
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

