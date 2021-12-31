#=============================================================================================================================#
#                                   Predicted prevalence functions                                                            #
#=============================================================================================================================#

#===============================#
# Simple Model (single dataset) #

# for MCMC 
predicted_prev_func <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  tp <-  1 - exp(-lambda * (data$age))   # true prevalence
  op <- (1-sp) + (se+sp-1)*tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

# for plotting predicted prevalence

predicted_prev_func2 <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]
  
  tp <-  1 - exp(-lambda * (age))   # true prevalence
  op <- (1-sp) + (se+sp-1)*tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}


calculate_predicted_prevalence_function <- function (max_age_toplot, data, pars, processed_chains) {
  
  # specify posterior median parameters to enable predicted prevalence calculation
  lambda.median <- pars[[4]]
  se.median <- pars[[5]]
  sp.median <- pars[[6]]
  
  # set up dataframe and calculate (median) predicted prevalence 
  age_dum <- seq(from=0, to = max_age_toplot, by=0.005)  ## If not already performed this step for simple catalytic model
  fitted_curve_df <- as.data.frame(age_dum) ## make sequence of numbers (mean ages) for predicted variable
  names(fitted_curve_df)[names(fitted_curve_df)=="age_dum"] <- "age"
  
  predicted_median_curve <- full_join(fitted_curve_df, data) 
  predicted_median_curve$predicted <- sapply(1:nrow(predicted_median_curve), 
                                             function(i) predicted_prev_func2(age = predicted_median_curve$age[i], 
                                                                                               c(lambda.median, se.median, 
                                                                                                 sp.median)))
  
  # create uncertainty (credible interval) of model run resulting from posterior #
  subsampled_model_outputs <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
  
  for (i in 1:length(processed_chains[[1]][,1])){
    
    single_model_output <- predicted_prev_func2(seq(0, max_age_toplot, 0.005),c(processed_chains[[1]][i,1],
                                                                       processed_chains[[1]][i,2],
                                                                       processed_chains[[1]][i,3]))
    subsampled_model_outputs[i, ] <- single_model_output
    
    }

  lower_credible_interval_processed <- apply(subsampled_model_outputs, MARGIN = 2, quantile, prob = 0.025)
  upper_credible_interval_processed <- apply(subsampled_model_outputs, MARGIN = 2, quantile, prob = 0.975)
  
  Lower_obs <- as.data.frame(lower_credible_interval_processed)
  Upper_obs <- as.data.frame(upper_credible_interval_processed)
  
  predicted_CrI <- cbind(Lower_obs, Upper_obs)
  predicted_CrI <- as.data.frame(predicted_CrI)
  predicted_CrI$age <- fitted_curve_df$age
  
  # plot predicted prevalence curve & uncertainty band vs data
  p <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,1.0)
  
  p2 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,0.5)
  
  p3 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,0.25)
  
  return(list(p, p2, p3, predicted_median_curve, predicted_CrI))

}


#===================================#
# reversible model (single dataset) #

predicted_prev_reversible_func <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  tp <-  (lambda/(lambda + rho)) * (1 - exp(-(lambda + rho) *(data$age)))   # true prevalence
  op <- (1-sp) + (se+sp-1) * tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

