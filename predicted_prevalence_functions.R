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

#=========================================#
# Simple Model (multiple dataset fitting) #

predicted_prev_func_multidataset <- function(age, par){
  # Sp and Se are the first to parameters in the vector of parameters
  sp <- par[1]
  se <- par[2]
  # Site-specific lamdas are the remianing parameters in the vector of parameters
  lambda_par <- par[3:length(par)]
  # Repeat each site specific lambda for each age group in each dataset
  lambda_all <- lambda_par[data$dataset]
  # Prediction
  tp <-  1 - exp(-lambda_all * data$age)  # true prevalence
  op <- (1-sp) + (se+sp-1)*tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

predicted_prev_func2_multidataset <- function(age, par){
  # Sp and Se are the first to parameters in the vector of parameters
  sp <- par[1]
  se <- par[2]
  # Site-specific lamdas are the remianing parameters in the vector of parameters
  lambda <- par[3]

  # Prediction
  tp <-  1 - exp(-lambda * age)  # true prevalence
  op <- (1-sp) + (se+sp-1)*tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}


#===================================#
# reversible model (single dataset) #

predicted_prev_reversible_func <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  tp <-  (lambda/(lambda + rho)) * (1 - exp(-(lambda + rho) *(data$age)))   # true prevalence
  op <- (1-sp) + (se+sp-1) * tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

predicted_prev_reversible_func2 <- function(age, par){
  
  lambda <- par[1]; se<- par[2]; sp <- par[3]; rho <- par[4]
  
  tp <-  (lambda/(lambda + rho)) * (1 - exp(-(lambda + rho) *(age)))   # true prevalence
  op <- (1-sp) + (se+sp-1) * tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}

#=====================================#
# reversible model (multiple dataset) #

# Model
predicted_prev_reversible_func_multidataset <- function(age, par){
  # Sp and Se are the first to parameters in the vector of parameters
  sp <- par[1]
  se <- par[2]
  # Site-specific lamdas are parameters 3: number of datasets in the vector of parameters
  lambda_par <- par[3:(length(unique(data$dataset))+2)]
  # Repeat each site specific lambda for each age group in each dataset
  lambda_all <- lambda_par[data$dataset]
  
  # Site-specific rhos are parameters 
  rho_par <- par[(2+length(unique(data$dataset))+1):(7+length(unique(data$dataset)))]
  # Repeat each site specific lambda for each age group in each dataset
  rho_all <- rho_par[data$dataset]
  
  # Prediction
  tp <-  (lambda_all/(lambda_all + rho_all)) *(1 - exp(-(lambda_all + rho_all) * (data$age)))  # true prevalence
  op <- (1-sp) + (se+sp-1)*tp      # observed prevalence Diggle et al 2011 Epi Res Int
  op
}


#================================================================#
# produce predicted prevalence curves (fitted to single dataset) #

calculate_predicted_prevalence_function <- function (max_age_toplot, data, pars, processed_chains, model) {
  
  if(model == "simple"){
  # specify posterior median parameters to enable predicted prevalence calculation
  lambda.median <- pars[[1]]
  se.median <- pars[[2]]
  sp.median <- pars[[3]]
  
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
  
  }
  
  if(model == "reversible"){
    # specify posterior median parameters to enable predicted prevalence calculation
    lambda.median <- pars[[1]]
    se.median <- pars[[2]]
    sp.median <- pars[[3]]
    rho.median <- pars[[4]]
    
    # set up dataframe and calculate (median) predicted prevalence 
    age_dum <- seq(from=0, to = max_age_toplot, by=0.005)  ## If not already performed this step for simple catalytic model
    fitted_curve_df <- as.data.frame(age_dum) ## make sequence of numbers (mean ages) for predicted variable
    names(fitted_curve_df)[names(fitted_curve_df)=="age_dum"] <- "age"
    
    predicted_median_curve <- full_join(fitted_curve_df, data) 
    predicted_median_curve$predicted <- sapply(1:nrow(predicted_median_curve), 
                                               function(i) predicted_prev_reversible_func2(
                                                 age = predicted_median_curve$age[i], 
                                                 c(lambda.median, se.median, sp.median, rho.median)))
    
    # create uncertainty (credible interval) of model run resulting from posterior #
    subsampled_model_outputs <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
    
    for (i in 1:length(processed_chains[[1]][,1])){
      
      single_model_output <- predicted_prev_reversible_func2(seq(0, max_age_toplot, 0.005),c(processed_chains[[1]][i,1],
                                                                                  processed_chains[[1]][i,2],
                                                                                  processed_chains[[1]][i,3],
                                                                                  processed_chains[[1]][i,4]))
      subsampled_model_outputs[i, ] <- single_model_output
      
    }
    
    lower_credible_interval_processed <- apply(subsampled_model_outputs, MARGIN = 2, quantile, prob = 0.025)
    upper_credible_interval_processed <- apply(subsampled_model_outputs, MARGIN = 2, quantile, prob = 0.975)
    
    Lower_obs <- as.data.frame(lower_credible_interval_processed)
    Upper_obs <- as.data.frame(upper_credible_interval_processed)
    
    predicted_CrI <- cbind(Lower_obs, Upper_obs)
    predicted_CrI <- as.data.frame(predicted_CrI)
    predicted_CrI$age <- fitted_curve_df$age
    
  }
  
  # plot predicted prevalence curve & uncertainty band vs data
  p <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,1.0)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  p2 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,0.5)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  p3 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower_credible_interval_processed,
                                       ymax=upper_credible_interval_processed), fill="purple", alpha=0.1)+
    ylim(0,0.25)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  return(list(p, p2, p3, predicted_median_curve, predicted_CrI))

}

#===============================================================================#
# produce predicted prevalence curves (simple model - fitted to single dataset) #


calculate_predicted_prevalence_multipledatasets_function1 <- function (max_age_toplot, data, pars, processed_chains, number_datasets) {
  
  if(number_datasets == 5){
    
    # specify posterior median parameters to enable predicted prevalence calculation
    sp.median <- pars[[1]]
    se.median <- pars[[2]]
    lambda1.median <- pars[[3]]
    lambda2.median <- pars[[4]]
    lambda3.median <- pars[[5]]
    lambda4.median <- pars[[6]]
    lambda5.median <- pars[[7]]
    
    # set up dataframes for predicted prevalence calculations (n = 5 datasets)
    n <- length(unique(data$dataset))
    eval(parse(text = paste0("data", seq(1:n), " <- ", split(data, data$dataset))))
    eval(parse(text = paste0("data", seq(1:n), " <-  as.data.frame(data", seq(1:number_datasets), ")"))) # change seq to 1: number of individual datasets
    
    age_dum <- seq(from=0, to=max_age_toplot, by=0.005)  
  
    fitted_curve_df <- as.data.frame(age_dum) # make sequence of numbers (mean ages) for predicted variable
     
    names(fitted_curve_df)[names(fitted_curve_df)=="age_dum"] <- "age"
     
    predicted_median_curve1 <- full_join(fitted_curve_df, data1) 
    predicted_median_curve2 <- full_join(fitted_curve_df, data2) 
    predicted_median_curve3 <- full_join(fitted_curve_df, data3)
    predicted_median_curve4 <- full_join(fitted_curve_df, data4)
    predicted_median_curve5 <- full_join(fitted_curve_df, data5)
    
    # calculate predicted prevalence for each dataset
    predicted_median_curve1$predicted <- sapply(1:nrow(predicted_median_curve1), function(i) predicted_prev_func2_multidataset(age = predicted_median_curve1$age[i], 
                                                                                                                              c(sp.median, se.median, lambda1.median)))
    predicted_median_curve2$predicted <- sapply(1:nrow(predicted_median_curve2), function(i) predicted_prev_func2_multidataset(age = predicted_median_curve2$age[i], 
                                                                                                                              c(sp.median, se.median, lambda2.median)))
    predicted_median_curve3$predicted <- sapply(1:nrow(predicted_median_curve3), function(i) predicted_prev_func2_multidataset(age = predicted_median_curve3$age[i], 
                                                                                                                              c(sp.median, se.median, lambda3.median)))
    predicted_median_curve4$predicted <- sapply(1:nrow(predicted_median_curve4), function(i) predicted_prev_func2_multidataset(age = predicted_median_curve4$age[i], 
                                                                                                                              c(sp.median, se.median, lambda4.median)))
    predicted_median_curve5$predicted <- sapply(1:nrow(predicted_median_curve5), function(i) predicted_prev_func2_multidataset(age = predicted_median_curve5$age[i], 
                                                                                                                              c(sp.median, se.median, lambda5.median)))
    # make master dataframe combining all datasets
    predicted_median_curve1$dataset_name <- rep(as.factor("data1"))
    predicted_median_curve2$dataset_name <- rep(as.factor("data2"))
    predicted_median_curve3$dataset_name <- rep(as.factor("data3"))
    predicted_median_curve4$dataset_name <- rep(as.factor("data4"))
    predicted_median_curve5$dataset_name <- rep(as.factor("data5"))
    
    
    predicted_median_curve <- rbind(predicted_median_curve1, predicted_median_curve2, predicted_median_curve3, predicted_median_curve4, predicted_median_curve5)
    
    # chains for each parameter to sample uncertainty from
    sp.chain <- c(PC_simple[[1]][,1], PC_simple[[2]][,1])
    se.chain <- c(PC_simple[[1]][,2], PC_simple[[2]][,2])
    lambda1.chain <- c(PC_simple[[1]][,3], PC_simple[[2]][,3])
    lambda2.chain <- c(PC_simple[[1]][,4], PC_simple[[2]][,4])
    lambda3.chain <- c(PC_simple[[1]][,5], PC_simple[[2]][,5])
    lambda4.chain <- c(PC_simple[[1]][,6], PC_simple[[2]][,6])
    lambda5.chain <- c(PC_simple[[1]][,7], PC_simple[[2]][,7])
    
    # calculate uncertainty (credible interval) of model run resulting from posterior #
    # dataset 1 #
    subsampled_model_outputs_dat1 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
    
    for (i in 1:length(processed_chains[[1]][,1])){
      
      single_model_output <- predicted_prev_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], 
                                                                                               sp.chain[i], 
                                                                                               lambda1.chain[i]))
      subsampled_model_outputs_dat1[i, ] <- single_model_output
      
    }
    
    lower_credible_interval_processed <- apply(subsampled_model_outputs_dat1, MARGIN = 2, quantile, prob = 0.025)
    upper_credible_interval_processed <- apply(subsampled_model_outputs_dat1, MARGIN = 2, quantile, prob = 0.975)
    
    Lower_obs <- as.data.frame(lower_credible_interval_processed)
    Upper_obs <- as.data.frame(upper_credible_interval_processed)
    
    predicted_CrI_dat1 <- cbind(Lower_obs, Upper_obs)
    predicted_CrI_dat1 <- as.data.frame(predicted_CrI_dat1)
    predicted_CrI_dat1$age <- fitted_curve_df$age
    predicted_CrI_dat1$dataset_name <- rep(as.factor("data1"))
    names(predicted_CrI_dat1) <- c("lower", "upper", "age", "dataset_name")
    
    # dataset 2 #
    subsampled_model_outputs_dat2 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
    
    for (i in 1:length(processed_chains[[1]][,1])){
      
      single_model_output <- predicted_prev_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], 
                                                                                               sp.chain[i], 
                                                                                               lambda2.chain[i]))
      subsampled_model_outputs_dat2[i, ] <- single_model_output
      
    }
    
    lower_credible_interval_processed <- apply(subsampled_model_outputs_dat2, MARGIN = 2, quantile, prob = 0.025)
    upper_credible_interval_processed <- apply(subsampled_model_outputs_dat2, MARGIN = 2, quantile, prob = 0.975)
    
    Lower_obs <- as.data.frame(lower_credible_interval_processed)
    Upper_obs <- as.data.frame(upper_credible_interval_processed)
    
    predicted_CrI_dat2 <- cbind(Lower_obs, Upper_obs)
    predicted_CrI_dat2 <- as.data.frame(predicted_CrI_dat2)
    predicted_CrI_dat2$age <- fitted_curve_df$age
    predicted_CrI_dat2$dataset_name <- rep(as.factor("data2"))
    names(predicted_CrI_dat2) <- c("lower", "upper", "age", "dataset_name")
    
    # dataset 3 #
    subsampled_model_outputs_dat3 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
    
    for (i in 1:length(processed_chains[[1]][,1])){
      
      single_model_output <- predicted_prev_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], 
                                                                                               sp.chain[i], 
                                                                                               lambda3.chain[i]))
      subsampled_model_outputs_dat3[i, ] <- single_model_output
      
    }
    
    lower_credible_interval_processed <- apply(subsampled_model_outputs_dat3, MARGIN = 2, quantile, prob = 0.025)
    upper_credible_interval_processed <- apply(subsampled_model_outputs_dat3, MARGIN = 2, quantile, prob = 0.975)
    
    Lower_obs <- as.data.frame(lower_credible_interval_processed)
    Upper_obs <- as.data.frame(upper_credible_interval_processed)
    
    predicted_CrI_dat3 <- cbind(Lower_obs, Upper_obs)
    predicted_CrI_dat3 <- as.data.frame(predicted_CrI_dat3)
    predicted_CrI_dat3$age <- fitted_curve_df$age
    predicted_CrI_dat3$dataset_name <- rep(as.factor("data3"))
    names(predicted_CrI_dat3) <- c("lower", "upper", "age", "dataset_name")
    
    # dataset 4 #
    subsampled_model_outputs_dat4 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
    
    for (i in 1:length(processed_chains[[1]][,1])){
      
      single_model_output <- predicted_prev_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], 
                                                                                               sp.chain[i], 
                                                                                               lambda4.chain[i]))
      subsampled_model_outputs_dat4[i, ] <- single_model_output
      
    }
    
    lower_credible_interval_processed <- apply(subsampled_model_outputs_dat4, MARGIN = 2, quantile, prob = 0.025)
    upper_credible_interval_processed <- apply(subsampled_model_outputs_dat4, MARGIN = 2, quantile, prob = 0.975)
    
    Lower_obs <- as.data.frame(lower_credible_interval_processed)
    Upper_obs <- as.data.frame(upper_credible_interval_processed)
    
    predicted_CrI_dat4 <- cbind(Lower_obs, Upper_obs)
    predicted_CrI_dat4 <- as.data.frame(predicted_CrI_dat4)
    predicted_CrI_dat4$age <- fitted_curve_df$age
    predicted_CrI_dat4$dataset_name <- rep(as.factor("data4"))
    names(predicted_CrI_dat4) <- c("lower", "upper", "age", "dataset_name")
    
    # dataset 5 #
    subsampled_model_outputs_dat5 <- matrix(NA, nrow = length(processed_chains[[1]][,1]), ncol = length(seq(0, 90, 0.005)))
    
    for (i in 1:length(processed_chains[[1]][,1])){
      
      single_model_output <- predicted_prev_func2_multidataset(seq(0, max_age_toplot, 0.005),c(sp.chain[i], 
                                                                                               sp.chain[i], 
                                                                                               lambda5.chain[i]))
      subsampled_model_outputs_dat5[i, ] <- single_model_output
      
    }
    
    lower_credible_interval_processed <- apply(subsampled_model_outputs_dat5, MARGIN = 2, quantile, prob = 0.025)
    upper_credible_interval_processed <- apply(subsampled_model_outputs_dat5, MARGIN = 2, quantile, prob = 0.975)
    
    Lower_obs <- as.data.frame(lower_credible_interval_processed)
    Upper_obs <- as.data.frame(upper_credible_interval_processed)
    
    predicted_CrI_dat5 <- cbind(Lower_obs, Upper_obs)
    predicted_CrI_dat5 <- as.data.frame(predicted_CrI_dat5)
    predicted_CrI_dat5$age <- fitted_curve_df$age
    predicted_CrI_dat5$dataset_name <- rep(as.factor("data5"))
    names(predicted_CrI_dat5) <- c("lower", "upper", "age", "dataset_name")
    
    # combined uncertainty dataframes across datasets
    predicted_CrI <- rbind(predicted_CrI_dat1, predicted_CrI_dat2, predicted_CrI_dat3, predicted_CrI_dat4, predicted_CrI_dat5)
    
    
  }
  
 # plot predicted prevalence curve & uncertainty band vs data
  p <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower, ymax=upper), fill="purple", alpha=0.1)+
    ylim(0,1.0)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    facet_wrap(~dataset_name, scales = "free")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  p2 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower, ymax=upper), fill="purple", alpha=0.1)+
    ylim(0,0.5)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    facet_wrap(~dataset_name, scales = "free")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  p3 <- ggplot() +   
    theme(legend.position = 'none') +   
    geom_point(data=predicted_median_curve, aes(x=age, y=prev))+
    geom_errorbar(data=predicted_median_curve,aes(x=age, y=prev, ymin=lower, ymax=upper), width=0.8)+
    geom_line(data=predicted_median_curve,aes(x=age, y=predicted), size= 1.1, colour='purple')+
    geom_ribbon(data=predicted_CrI,aes(x=age, ymin=lower, ymax=upper), fill="purple", alpha=0.1)+
    ylim(0,0.25)+
    labs(x="Age (months) of human host", y="(Sero)prevalence (0 -1)")+
    facet_wrap(~dataset_name, scales = "free")+
    theme_bw() +
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text=element_text(size=11, face= "bold"),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          plot.title = element_text(size = 18, hjust=0.5),
          axis.title.x = element_text(size = 18, face= "bold"),
          axis.title.y = element_text(size = 16, angle = 90, face= "bold"),
          legend.position = c(0.83, 0.15),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=16))
  
  
  return(list(p, p2, p3, predicted_median_curve, predicted_CrI))
  
}

