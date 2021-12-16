#=============================================================================================================================#
#                                Estimating diagnostic (prior) parameters                                                     #
#=============================================================================================================================#

estimate_alpha_beta_par_diagnostic <-
  function(input_alpha, input_beta, target_l, target_u) {
    
    # Define targets we want to fit to
    #test_alpha <- 84.5 # original
    #test_beta <- 15.5 # original
    
    # mean (of beta prior distribution given by alpha and beta shape parameter as:)
    target_mean <- input_alpha / (input_alpha + input_beta)
    
    # Lower CI
    #target_l <- qbeta(0.025, input_alpha, input_beta)
    # Upper CI
    #target_u <- qbeta(0.975, input_alpha, input_beta)
    
    # Actual upper and low targets to calibrate on for below #
    target_l <- target_l
    target_u <- target_u
    
    # Find parameters
    # Vector of possible alphas
    alpha <- seq(0, 100, 0.1)
    # Vector of corresponding betas, such that alpha / (alpha + beta) = m
    beta <- (alpha - (alpha * target_mean)) / target_mean
    # Lower CI | alpha, beta
    l <- qbeta(0.025, alpha, beta)
    # Upper CI | alpha, beta
    u <- qbeta(0.975, alpha, beta)
    # Best fit
    best <- which.min(abs(l - target_l) + abs(u - target_u))
    
    # Recover input parameters
    alpha[best]
    beta[best]
    
    actual_alpha <- alpha[best] # updated, best fit
    actual_beta <- beta[best] # updated, best fit
    
    actual_mean <- actual_alpha / (actual_alpha + actual_beta) # mean
    actual_mean
    actual_lower <- qbeta(0.025, actual_alpha, actual_beta)
    actual_lower
    actual_upper <- qbeta(0.975, actual_alpha, actual_beta)
    actual_upper
    
    # 0.6592385 (l) and  0.9939759 (u) # given by 11 for alpha and 1.2612613 for beta
    # p = seq(0, 1, length = 100)
    # #par(mfrow = c(1, 1)) #
    # 
    # a <- plot(
    #   p,
    #   dbeta(p, actual_alpha, actual_beta),
    #   ylab = "density",
    #   type = "l",
    #   col = 4
    # ) #
    
    

    return(list(actual_alpha, actual_beta))
  }
