#==============================================================================================================================#
#                                       MASTER SCRIPT  - multiple dataset fitting                                              #
#==============================================================================================================================#

rm(list = ls())

#=============================================#
#       Load data files                       #
#=============================================#

#=======================#                                                                                             
# Initiatie sub-scripts #                                                                                             
source('libraries.R')
source('plot_data_modelfit_functions.R')
source('diagnostic_parameter_functions.R')
source('predicted_prevalence_functions.R')
source('MCMC_functions_singledataset.R')
source('plot_MCMC_output_functions.R')
source('process_MCMC_output_functions.R')

#======================#
#    load data         #
#======================#

test_data_multipledataset <- read.csv("~/human_tsol_FoI_modelling/data/test_data_multipledataset.csv")

data <- test_data_multipledataset

data <- data[,c(which(colnames(data)=="ID"),which(colnames(data)!="ID"))]

names(data) <- c("dataset", "age", "pos", "n","prev","lower","upper","ref")

data$ref <- as.factor(data$ref)

plot_ageprev_func2(data) # plot age-prevalence data

#===============================================#
#  optimise parameters for diagnostic priors    #

# sensitivity #
sensitivity_parameters <- estimate_alpha_beta_par_diagnostic(input_alpha = 97, input_beta = 3, 
                                                             target_l = 0.95, target_u = 1.00)

alpha_se <- sensitivity_parameters[[1]]
alpha_se 
beta_se <- sensitivity_parameters[[2]]
beta_se

# specifcity #
specificity_parameters <- estimate_alpha_beta_par_diagnostic(input_alpha = 97, input_beta = 3, 
                                                             target_l = 0.94, target_u = 0.99)

alpha_sp <- specificity_parameters[[1]]
alpha_sp
beta_sp <- specificity_parameters[[2]]
beta_sp

#===========================================================================================================================#
#==================================================#
#  run MCMC (multiple datasets; simple FoI model)  #

inits1 <- c(0.95, 0.955, 0.001, 0.1,0.01,0.01,0.1)   # Initial parameter values to initiate chains
inits2 <- c(0.99, 0.995, 0.01, 0.01, 0.1, 0.1,0.01) # e.g. Ab−EITB, rT24H (se: 0.96 (0.93-0.99), sp: 0.98 (0.96-1), Noh et al. 2014) 
sd <- 0.0003 # set standard deviation of proposal distribution; aim for 0.25 acceptance
cov <- diag(sd^2, 2+length(unique(data$dataset))) # covariance
niter <- 100000 # number of iterations
burnin <- 50000 # burnin (chains to discard before convergence)

# run MCMC (chain 1)
set.seed(123) # for reproducibility
simple_out_chain1 <- MCMC_simple_model(inits1, niter, cov, fitting = "multiple datasets")  # initiate the MCMC

# run MCMC (chain 1)
set.seed(123)
simple_out_chain2 <- MCMC_simple_model(inits2, niter, cov, fitting = "multiple datasets")  # initiate the MCMC

# whats the acceptance ratio (aiming for 0.25)
sum(simple_out_chain1$Acceptances)/niter
sum(simple_out_chain2$Acceptances)/niter

# plot MCMC outputs #
chains_plot <- chains_plot_multidatasets_func(inits1 = inits1, chain1 = simple_out_chain1, chain2 = simple_out_chain2,
                                model = "simple", number_datasets = 5)

histograms_plot <- histogram_plot_multidatasets_func(inits1 = inits1, burnin = burnin, niter = niter, 
                                       chain1 = simple_out_chain1, chain2 = simple_out_chain2, model = "simple",
                                       number_datasets = 5)

#=============================================# 
#         process MCMC outputs                # 

chains1_output <- simple_out_chain1$MCMC_Output
chains2_output <- simple_out_chain2$MCMC_Output

# remove burnin and proceed with reducing autocorrelation (thinning by sub-sampling)
PC_simple <- Process_chains(chains1_output, chains2_output, burnin = burnin, sample = 100) # set burnin to 0 if already

# View the process chains (from the autocorrelation plots, sampling every 20th value seems appropriate)
plot_chains_multidatasets(PC_simple[[1]], PC_simple[[2]], inits = inits1, number_datasets = 5)

# check autocorrelation of chains for each parameter (to inform sub-sampling)
check_autocorr <- determine_autocorrelation_func1(processed_chain = PC_simple, number_datasets = 5)

check_autocorr[1] # autocorrelation significance parameter 1 
check_autocorr[2] # autocorrelation significance parameter 2 
check_autocorr[3] # autocorrelation significance parameter 3 
check_autocorr[4] # autocorrelation significance parameter 4 
check_autocorr[5] # autocorrelation significance parameter 5 
check_autocorr[6] # autocorrelation significance parameter 6 
check_autocorr[7] # autocorrelation significance parameter 7 

# plot loglikelihood
loglikchains_plot_func(chain1 = simple_out_chain1, chain2 = simple_out_chain2)

#==============================================================================#
# Obtain parameter values (median & credible) & plot posterior distributions   #

simple_model_parameters <- obtain_parameter_values_func(processed_chains = PC_simple, model = "simple", number_datasets = 5)
simple_model_parameters

plot_posterior_distrib_func(processed_chains = PC_simple, model = "simple", number_datasets = 5)

#============================================================#
# calculate Deviance Information Criterion (DIC) - model fit #

DIC_result <- calculate_DIC_func(chain = simple_out_chain1, burnin = burnin, subsample = 50, 
                                 parameters = simple_model_parameters, model = "simple", number_datasets = 5)
DIC_result  # 1) D bar model1, 2) modal posterior likelihood, 3) modal posterior deviance, 4) DIC

#====================================================================================================================#
# calculate (with posterior parameter estimates) predicted (sero)prevalence and unceetainty intervals (for plotting) #

predicted_prev_output <- calculate_predicted_prevalence_multipledatasets_function1(max_age_toplot = 90, data = data, 
                                                                 pars = simple_model_parameters,
                                                                 processed_chains = PC_simple, number_datasets = 5)

predicted_prev_output[[1]] # predicted prevalence plot (ylim 0-100% prev)
predicted_prev_output[[2]] # predicted prevalence plot (ylim 0-50% prev)
predicted_prev_output[[3]] # predicted prevalence plot (ylim 0-25% prev)

#===========================================================================================================================#
#======================================================#
#  run MCMC (multiple datasets; reversible FoI model)  #

inits1 <- c(0.99, 0.95, 0.00003, 0.001,0.005,0.0001,0.0001,0.1,0.1,0.1,0.1,0.1)  # Initial parameter values to initiate chains
inits2 <- c(0.96, 0.99, 0.001, 0.04, 0.2, 0.04,0.04,0.01,0.01,0.01,0.01,0.01) # e.g. Ab−EITB, rT24H (se: 0.96 (0.93-0.99), sp: 0.98 (0.96-1), Noh et al. 2014) 
sd <- 0.0006 # set standard deviation of proposal distribution; aim for 0.25 acceptance
cov <- diag(sd^2, 2+(2*length(unique(data$dataset)))) # covariance
niter <- 1000000 # number of iterations
burnin <- 500000 # burnin (chains to discard before convergence)

# run MCMC (chain 1)
set.seed(123) # for reproducibility
reversible_out_chain1 <- MCMC_reversible_model(inits1, niter, cov, simple_lambda_median = c(0.00007756672,
                                                                                            0.003205126,
                                                                                            0.008291567,
                                                                                            0.003112397,
                                                                                            0.003481667), 
                                               fitting = "multiple datasets")  # initiate the MCMC

# run MCMC (chain 1)
set.seed(123)
reversible_out_chain2 <- MCMC_reversible_model(inits2, niter, cov, simple_lambda_median = c(0.00007756672,
                                                                                            0.003205126,
                                                                                            0.008291567,
                                                                                            0.003112397,
                                                                                            0.003481667),
                                               fitting = "multiple datasets")  # initiate the MCMC

# whats the acceptance ratio (aiming for 0.25)
sum(reversible_out_chain1$Acceptances)/niter
sum(reversible_out_chain2$Acceptances)/niter

# plot MCMC outputs #
chains_plot <- chains_plot_func(inits1 = inits1, chain1 = reversible_out_chain1, chain2 = reversible_out_chain2, 
                                model = "reversible")

histograms_plot <- histogram_plot_func(inits1 = inits1, burnin = burnin, niter = niter, 
                                       chain1 = reversible_out_chain1, chain2 = reversible_out_chain2,
                                       model = "reversible")

#=============================================# 
#         process MCMC outputs                # 

chains1_output <- reversible_out_chain1$MCMC_Output
chains2_output <- reversible_out_chain2$MCMC_Output

# remove burnin and proceed with reducing autocorrelation (thinning by sub-sampling)
PC_reversible <- Process_chains(chains1_output, chains2_output, burnin = burnin, sample = 50) # set burnin to 0 if already

# View the process chains (from the autocorrelation plots, sampling every 20th value seems appropriate)
plot_chains(PC_reversible[[1]], PC_reversible[[2]])

# check autocorrelation of chains for each parameter (to inform sub-sampling)
check_autocorr <- determine_autocorrelation_func2(processed_chain = PC_reversible, number_datasets = 1)

check_autocorr[1] # autocorrelation significance parameter 1 (e.g. lambda)
check_autocorr[2] # autocorrelation significance parameter 2 (e.g. se)
check_autocorr[3] # autocorrelation significance parameter 1 (e.g. sp)
check_autocorr[4] # autocorrelation significance parameter 1 (e.g. sp)

# plot loglikelihood
loglikchains_plot_func(chain1 = reversible_out_chain1, chain2 = reversible_out_chain2)

#==============================================================================#
# Obtain parameter values (median & credible) & plot posterior distributions   #

reversible_model_parameters <- obtain_parameter_values_func(processed_chains = PC_reversible, model = "reversible", number_datasets = 5)
reversible_model_parameters

plot_posterior_distrib_func(processed_chains = PC_reversible, model = "reversible", number_datasets = 1)

#============================================================#
# calculate Deviance Information Criterion (DIC) - model fit #

DIC_result <- calculate_DIC_func(chain = reversible_out_chain1, burnin = burnin, subsample = 50, 
                                 parameters = reversible_model_parameters, model = "reversible", number_datasets = 1)
DIC_result  # 1) D bar model1, 2) modal posterior likelihood, 3) modal posterior deviance, 4) DIC

#====================================================================================================================#
# calculate (with posterior parameter estimates) predicted (sero)prevalence and unceetainty intervals (for plotting) #

predicted_prev_output <- calculate_predicted_prevalence_function(max_age_toplot = 90, data = data, 
                                                                 pars = reversible_model_parameters,
                                                                 processed_chains = PC_reversible, model = "reversible")

predicted_prev_output[[1]] # predicted prevalence plot (ylim 0-100% prev)
predicted_prev_output[[2]] # predicted prevalence plot (ylim 0-50% prev)
predicted_prev_output[[3]] # predicted prevalence plot (ylim 0-25% prev)


