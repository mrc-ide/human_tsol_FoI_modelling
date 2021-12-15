#=============================================================================================================================#
#                                      Process MCMC outputs functions                                                         #
#=============================================================================================================================#

#===================================================#
#  Simple FoI moodel: single dataset                #


plot_chains<-function(run1, run2){
  par(mfrow=c(ncol(run1),1))
  
  for(i in 1:ncol(run1)){
    plot(run1[,i], t='l', col='deeppink',
         ylim=c(min(c(run1[,i], run2[,i])),max(c(run1[,i], run2[,i]))),
         xlab='', ylab=paste('Parameter', i, sep=' '))
    lines(run2[,i], col='dodgerblue')
  }
  
}

Burn<-function(chains, burnin){
  chains[-(1:burnin),]
}

Downsample<-function(chains, sample){
  chains[seq(1, nrow(chains), sample),]
}

Process_chains<-function(run1, run2, burnin, sample){
  C1<-Burn(run1, burnin)
  C2<-Burn(run2, burnin)
  
  
  C1<-Downsample(C1, sample)
  C2<-Downsample(C2, sample)
  
  
  return(list(C1, C2))
}
