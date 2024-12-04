#################################
# Same log-likelihood fn as parametric models
#################################

log_lik_fn_mixture <- function(params, data, fn) {
  x <- data$Age
  Dx <- data$Deaths
  Ex <- data$Exposure
  mux <- fn(params, data)
  # add small constant to avoid log(0)
  mux <- pmax(mux, .Machine$double.eps)
  
  # log lik assuming deaths are poisson distributed
  log_lik <- sum( (Dx*log(mux)) - (Ex*mux) )
  
  return(log_lik) #* 20^-6) # scale down according to HMD procedure 
}

#################################
# Models for localized flexibility
#################################

# First model type: mixture of gammas, as recommended by Wrigley-Field (2014)
# Actually trying mixture of gamma-Gompertz here but component functions can 
# easily be adjusted if needed

mixture_fn <- function(params, data) {
  w1 <- params[1]
  w2 <- 1-w1
  theta1 <- params[2:4] # parameters for first g-g model
  theta2 <- params[5:7] # parameters for second g-g model
  
  x <- data$Age
  f1 <- (theta1[1] * exp(theta1[2] * (x - 80))) / (1 + theta1[3] * exp(theta1[2] * (x - 80)))  # Component 1
  f2 <- (theta2[1] * exp(theta2[2] * (x - 80))) / (1 + theta2[3] * exp(theta2[2] * (x - 80)))  # Component 2
  
  mixture_rate <- w1 * f1 + w2 * f2
  
  return(mixture_rate)
}

mixture_get_mx <- function(params) {
  w1 <- params[1]
  w2 <- 1-w1
  theta1 <- params[2:4] # parameters for first g-g model
  theta2 <- params[5:7] # parameters for second g-g model
  
  x <- 80:115
  f1 <- (theta1[1] * exp(theta1[2] * (x - 80))) / (1 + theta1[3] * exp(theta1[2] * (x - 80)))  # Component 1
  f2 <- (theta2[1] * exp(theta2[2] * (x - 80))) / (1 + theta2[3] * exp(theta2[2] * (x - 80)))  # Component 2
  
  mixture_rate <- w1 * f1 + w2 * f2
  Mx <- cbind(x, mixture_rate)
  
  return(Mx)
}

# initial params will use random samples rather than grid search to conserve memory
mixture_get_init_params <- function(log_lik_fn, data, fn, 
                                    n_samples = 1000,  # Number of random samples
                                    w_range = c(0.1, 0.9), 
                                    a_range = c(.Machine$double.eps, 1),
                                    b_range = c(.Machine$double.eps, 1),
                                    d_range = c(.Machine$double.eps, 1)) {
  # Generate random samples for parameters
  random_samples <- data.frame(
    w = runif(n_samples, w_range[1], w_range[2]),
    a1 = runif(n_samples, a_range[1], a_range[2]),
    b1 = runif(n_samples, b_range[1], b_range[2]),
    d1 = runif(n_samples, d_range[1], d_range[2]),
    a2 = runif(n_samples, a_range[1], a_range[2]),
    b2 = runif(n_samples, b_range[1], b_range[2]),
    d2 = runif(n_samples, d_range[1], d_range[2])
  )
  
  # Evaluate log-likelihood for each sample
  random_samples$log_lik <- apply(random_samples, 1, function(row) {
    params <- c(row['w'], row['a1'], row['b1'], row['d1'], row['a2'], row['b2'], row['d2'])
    return(log_lik_fn(params, data, fn))
  })
  
  # Find the combination with the highest log-likelihood
  best_params <- random_samples[which.max(random_samples$log_lik), ]
  initial_params <- unlist(best_params)
  return(initial_params)
}


fit_mixture <- function(init_param_fn, log_lik_fn, fn, mx_fn, data) {
  initial_params <- mixture_get_init_params(log_lik_fn = log_lik_fn, data = data, fn = fn)
  optim_result <- optim(par = initial_params, 
                        fn = function(params) log_lik_fn(params, data, fn),
                        method = "L-BFGS-B", 
                        lower = c(0, rep(.Machine$double.eps, 6)),  # Lower bounds
                        upper = c(1, rep(5, 6)),                   # Upper bounds
                        control = list(fnscale=-1))
  
  fitted_mx <- mx_fn(optim_result$par)
  opt_params <- optim_result$par
  log_lik_opt <- log_lik_fn(opt_params, data, fn)
  num_params <- length(opt_params)
  
  n <- nrow(data)
  aic <- 2 * num_params - 2 * log_lik_opt
  bic <- log(n) * num_params - 2 * log_lik_opt
  
  results <- list(fitted_mx = fitted_mx, AIC = aic, BIC = bic)
  return(results)
}


apply_mixture_model <- function(data) {
  result <- fit_mixture(
    init_param_fn = mixture_get_init_params, 
    log_lik_fn = log_lik_fn_mixture, 
    fn = mixture_fn,
    mx_fn = mixture_get_mx, 
    data = data
  )
  return(result)
}

fit_all_mixture <- function(data) {
  data |>
    group_by(Cohort, Country, Sex) |>
    do({
      model_result <- apply_mixture_model(.)
      tibble(
        fitted_rates = list(model_result$fitted_mx),
        AIC = model_result$AIC,
        BIC = model_result$BIC
      )
    }) |>
    ungroup()
}

