#==============================================================================================================================#
#                                       MASTER SCRIPT                                                                          #
#==============================================================================================================================#

rm(list = ls())

#=============================================#
#       Load data files                       #
#=============================================#

#=======================#                                                                                             
# Initiatie sub-scripts #                                                                                             
source('libraries.R')
source('diagnostic_parameter_functions.R')
source('predicted_prevalence_functions.R')
source('MCMC_functions_singledataset.R')
source('plot_MCMC_output_functions.R')
source('process_MCMC_output_functions.R')

#======================#
#    load data         #
#======================#

test_data_singledataset <- read.csv("~/human_tsol_FoI_modelling/data/test_data_singledataset.csv")

data <- test_data_singledataset

names(data) <- c("age", "pos_pigs", "n_pigs","prev","lower","upper")

#===============================================#
#  optimise parameters for diagnostic priors    #

# sensitivity #
sensitivity_parameters <- estimate_alpha_beta_par_diagnostic(input_alpha = 96, input_beta = 4, 
                                                             target_l = 0.93, target_u = 0.99)

alpha_se <- sensitivity_parameters[[1]]
alpha_se 
beta_se <- sensitivity_parameters[[2]]
beta_se

# specifcity #
specificity_parameters <- estimate_alpha_beta_par_diagnostic(input_alpha = 98, input_beta = 2, 
                                                             target_l = 0.96, target_u = 1)

alpha_sp <- specificity_parameters[[1]]
alpha_sp
beta_sp <- specificity_parameters[[2]]
beta_sp

#===============================================#
#  run MCMC (single dataset; simple FoI model)  #

inits1 <- c(0.1, 0.93, 0.96)   # Initial parameter values to initiate chains
inits2 <- c(0.0001, 0.99, 0.999) # e.g. Ab−EITB, rT24H (se: 0.96 (0.93-0.99), sp: 0.98 (0.96-1), Noh et al. 2014) 
sd <- 0.001 # set standard deviation of proposal distribution; aim for 0.25 acceptance
cov <- diag(sd^2, 3) # covariance
niter <- 100000 # number of iterations
burnin <- 50000 # burnin (chains to discard before convergence)

# run MCMC (chain 1)
set.seed(123) # for reproducibility
simple_out_chain1 <- MCMC_simple_model(inits1, niter, cov)  # initiate the MCMC

# run MCMC (chain 1)
set.seed(123)
simple_out_chain2 <- MCMC_simple_model(inits2, niter, cov)  # initiate the MCMC

# whats the acceptance ratio (aiming for 0.25)
sum(simple_out_chain1$Acceptances)/niter
sum(simple_out_chain2$Acceptances)/niter

# plot MCMC outputs #
chains_plot <- chains_plot_func1(inits1 = inits1, chain1 = simple_out_chain1, chain2 = simple_out_chain2)

histograms_plot <- histogram_plot_func1(inits1 = inits1, burnin = burnin, niter = niter, 
                                        chain1 = simple_out_chain1, chain2 = simple_out_chain2)

#=============================================# 
#         process MCMC outputs                # 

chains1_output <- simple_out_chain1$MCMC_Output
chains2_output <- simple_out_chain2$MCMC_Output

# remove burnin and proceed with reducing autocorrelation (thinning by sub-sampling)
PC_simple <-Process_chains(chains1_output, chains2_output, burnin=50000, sample=50) # set burnin to 0 if already

# View the process chains (from the autocorrelation plots, sampling every 20th value seems appropriate)
plot_chains(PC_simple[[1]], PC_simple[[2]])

check_autocorr <- determine_autocorrelation_func(processed_chain = PC_simple)

check_autocorr[1] # autocorrelation significance parameter 1 (e.g. lambda)
check_autocorr[2] # autocorrelation significance parameter 2 (e.g. se)
check_autocorr[3] # autocorrelation significance parameter 1 (e.g. sp)

