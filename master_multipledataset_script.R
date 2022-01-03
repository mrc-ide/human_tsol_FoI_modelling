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

names(data) <- c("dataset", "age", "pos", "n_","prev","lower","upper","ref")

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
