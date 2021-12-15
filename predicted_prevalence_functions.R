#=============================================================================================================================#
#                                   Predicted prevalence functions                                                            #
#=============================================================================================================================#

#===============================#
# Simple Model (single dataset) #

predicted_prev_func <- function(data, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  tp <-  1 - exp(-lambda * (data$age))   # true prevalence
  op <- (1-sp) + (se+sp-1)*tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

#predicted_prev_func(data, c(0.1,0.9,0.9)) 