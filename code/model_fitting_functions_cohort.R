##########################
# Generic functions 
##########################

log_lik_fn <- function(params, data, fn) {
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

# ##########################
# # Global variables
# ##########################
# 
# low_age = 80
# fra_high_age = 115
# scand_high_age = 108

##########################
# Model specific functions 
##########################

#### Kannisto ####
get_initial_params <- function(log_lik_fn, data, fn,
                               a_values=seq(0.001, 1, length.out = 20), 
                               b_values=seq(0.001, 1, length.out = 20)) {
  # Perform grid search
  grid_search_results <- expand.grid(a = a_values, b = b_values)
  
  # Evaluate log-likelihood for each combination of a and b
  grid_search_results$log_lik <- apply(grid_search_results, 1, function(row) {
    params <- c(row['a'], row['b'])
    return(log_lik_fn(params, data, fn))
  })
  
  # Find the combination with the highest log-likelihood 
  best_params <- grid_search_results[which.max(grid_search_results$log_lik), ]
  
  # Extract the best initial values for a and b
  initial_params <- c(a = best_params$a, b = best_params$b)
  
  return(initial_params)
}

kannisto_fn <- function(params, data) {
  a <- params[1]
  b <- params[2]
  x <- data$Age
  
  ( a*exp(b*( (x+0.5) - 80) ) ) / ( 1 + a*exp(b*( (x+0.5) - 80)) )
}


kannisto_get_grad <- function(params, data, fn, epsilon = 1e-8) {
  a <- params[1]
  b <- params[2]
  
  x <- data$Age
  Dx <- data$Deaths
  Ex <- data$Exposure
  
  Q = b * (x + 0.5 - 80)
  P = exp(Q)
  N = a * P
  mux = N / (1 + N)
  M = log(mux)
  G = Dx * M - Ex * mux
  L = sum(G)
  
  dLdG = rep(1, length(G))
  dGdmux = Dx / mux - Ex
  dmuxdN = 1 / (1+N)^2
  dNda = P
  dNdP = a
  dPdQ = exp(Q)
  dQdb = x + 0.5 - 80
  
  grad_a = sum( (dLdG * dGdmux * dmuxdN * dNda) )
  grad_b = sum( (dLdG * dGdmux * dmuxdN * dNdP * dPdQ * dQdb) )
  
  grad <- c(grad_a, grad_b) #* 10^-6 
  
  if (any(!is.finite(grad))) {
       print(params)
       print(grad)
       stop("Non-finite gradient value encountered")
     }

  return(grad)
}

kannisto_get_mx <- function(a,b) {
  
  x <- 80:115
  Mx <- ( a*exp(b*( (x+0.5) - 80) ) ) / ( 1 + a*exp(b*( (x+0.5) - 80)) )
  Mx <- cbind(x, Mx)
  
  return(Mx)
}

fit_kannisto <- function(init_param_fn, log_lik_fn, grad_fn, fn, data) {
  initial_params <- get_initial_params(log_lik_fn = log_lik_fn, data = data, fn = fn)
  optim_result <- optim(par = initial_params, 
                        fn = function(params) log_lik_fn(params, data, fn),
                        gr = function(params) grad_fn(params, data, fn), 
                        method = "L-BFGS-B", 
                        lower = c(.Machine$double.eps, .Machine$double.eps), 
                        upper = c(5, 5),
                        control = list(fnscale=-1))
  fitted_mx <- kannisto_get_mx(optim_result$par[1], optim_result$par[2])
  opt_params <- optim_result$par
  log_lik_opt <- log_lik_fn(opt_params, data, fn)
  num_params <- length(opt_params)
  
  n <- nrow(data)
  aic <- 2 * num_params - 2 * log_lik_opt
  bic <- log(n) * num_params - 2 * log_lik_opt
  results <- list(fitted_mx = fitted_mx, AIC = aic, BIC = bic)
  
  return(results)
}


apply_kannisto_model <- function(data) {
  result <- fit_kannisto(
    init_param_fn = get_initial_params, 
    log_lik_fn = log_lik_fn, 
    grad_fn = kannisto_get_grad, 
    fn = kannisto_fn, 
    data = data
  )
  
  return(result)
}

fit_all_kannisto <- function(data) {
  data |>
    group_by(Cohort, Country, Sex) |>
    do({
      model_result <- apply_kannisto_model(.)
      tibble(
        fitted_rates = list(model_result$fitted_mx),
        AIC = model_result$AIC,
        BIC = model_result$BIC
      )
    }) |>
    ungroup()
}

#### Beard ####
beard_get_init_params <- function(log_lik_fn, data, fn,
                                  a_values=seq(.Machine$double.eps, 1, length.out = 20), 
                                  b_values=seq(.Machine$double.eps, 1, length.out = 20),
                                  d_values=seq(.Machine$double.eps, 1, length.out = 20)) {
  # Perform grid search
  grid_search_results <- expand.grid(a = a_values, b = b_values, d = d_values)
  
  # Evaluate log-likelihood for each combination of a and b
  grid_search_results$log_lik <- apply(grid_search_results, 1, function(row) {
    params <- c(row['a'], row['b'], row['d'])
    return(log_lik_fn(params, data, fn))
  })
  
  # Find the combination with the highest log-likelihood (least negative)
  best_params <- grid_search_results[which.max(grid_search_results$log_lik), ]
  
  # Extract the best initial values for a and b
  initial_params <- c(a = best_params$a, b = best_params$b, d = best_params$d)
  
  return(initial_params)
}
  
beard_fn <- function(params, data) {
  a <- params[1]
  b <- params[2]
  d <- params[3]
  x <- data$Age
  
  ( a*exp(b*( (x+0.5) - 80) ) ) / ( 1 + d*exp(b*( (x+0.5) - 80)) )
}

beard_get_grad <- function(params, data, fn, epsilon = 1e-8) {
  a <- params[1]
  b <- params[2]
  d <- params[3]
  
  x <- data$Age
  Dx <- data$Deaths
  Ex <- data$Exposure
  
  Q = b * (x + 0.5 - 80)
  P = exp(Q)
  N = a * P
  R = d * P
  mux = N / (1 + R)
  M = log(mux)
  G = Dx * M - Ex * mux
  L = sum(G)
  
  dLdG = rep(1, length(G))
  dGdmux = Dx / mux - Ex
  dmuxdN = 1 / (1+N)^2
  dmuxdR = -N / (1+R)^2
  dRdd = P
  dNda = P
  dNdP = a
  dPdQ = exp(Q)
  dQdb = x + 0.5 - 80
  
  grad_a = sum( (dLdG * dGdmux * dmuxdN * dNda) )
  grad_b = sum( (dLdG * dGdmux * dmuxdN * dNdP * dPdQ * dQdb) )
  grad_d = sum( (dLdG * dGdmux * dmuxdR * dRdd) )
  
  grad <- c(grad_a, grad_b, grad_d) #* 10^-6 
  
  if (any(!is.finite(grad))) {
    print(params)
    print(grad)
    stop("Non-finite gradient value encountered")
  }
  
  return(grad)
}

beard_get_mx <- function(a,b,d) {
  
  x <- 80:115
  Mx <- ( a*exp(b*( (x+0.5) - 80) ) ) / ( 1 + d*exp(b*( (x+0.5) - 80)) )
  Mx <- cbind(x, Mx)
  
  return(Mx)
}

fit_beard <- function(init_param_fn, log_lik_fn, grad_fn, fn, data) {
  initial_params <- beard_get_init_params(log_lik_fn = log_lik_fn, data = data, fn = fn)
  optim_result <- optim(par = initial_params, 
                        fn = function(params) log_lik_fn(params, data, fn),
                        gr = function(params) grad_fn(params, data, fn), 
                        method = "L-BFGS-B", 
                        lower = c(.Machine$double.eps, .Machine$double.eps), 
                        upper = c(5, 5),
                        control = list(fnscale=-1))
  # high_age <- case_when(
  #   data$Country == "fra" ~ fra_high_age,
  #   data$Country != "fra" ~ scand_high_age
  # )
  # high_age <- high_age[1]
  
  fitted_mx <- beard_get_mx(optim_result$par[1], optim_result$par[2], optim_result$par[3])
  opt_params <- optim_result$par
  log_lik_opt <- log_lik_fn(opt_params, data, fn)
  num_params <- length(opt_params)
  
  n <- nrow(data)
  aic <- 2 * num_params - 2 * log_lik_opt
  bic <- log(n) * num_params - 2 * log_lik_opt
  # data <- data |>
  #   mutate(observed_rates = Deaths/Exposure)
  # observed <- data$observed_rates
  # fitted <- fitted_mx[, 2]
  # sse <- sum((observed - fitted)^2)
  
  results <- list(fitted_mx = fitted_mx, AIC = aic, BIC = bic)
  
  return(results)
}

apply_beard_model <- function(data) {
  #fit_beard(
  result <- fit_beard(
    init_param_fn = beard_get_init_params, 
    log_lik_fn = log_lik_fn, 
    grad_fn = beard_get_grad, 
    fn = beard_fn, 
    data = data
  )
  
  return(result)
}

# fit_all_beard <- function(data) {
#   data |>
#     group_by(Cohort, Country, Sex) |>
#     do(fitted_rates = apply_beard_model(.)) |>
#     ungroup()
# }

fit_all_beard <- function(data) {
  data |>
    group_by(Cohort, Country, Sex) |>
    do({
      model_result <- apply_beard_model(.)
      tibble(
        fitted_rates = list(model_result$fitted_mx),
        AIC = model_result$AIC,
        BIC = model_result$BIC
      )
    }) |>
    ungroup()
}




#### Gompertz ####
gompertz_get_init_params <- function(log_lik_fn, data, fn,
                                  a_values=seq(0.001, 1, length.out = 20), 
                                  b_values=seq(0.001, 1, length.out = 20)) {
  # Perform grid search
  grid_search_results <- expand.grid(a = a_values, b = b_values)
  
  # Evaluate log-likelihood for each combination of a and b
  grid_search_results$log_lik <- apply(grid_search_results, 1, function(row) {
    params <- c(row['a'], row['b'])
    return(log_lik_fn(params, data, fn))
  })
  
  # Find the combination with the highest log-likelihood (least negative)
  best_params <- grid_search_results[which.max(grid_search_results$log_lik), ]
  
  # Extract the best initial values for a and b
  initial_params <- c(a = best_params$a, b = best_params$b)
  
  return(initial_params)
}

gompertz_fn <- function(params, data) {
  a <- params[1]
  b <- params[2]
  x <- data$Age
  
  ( a * exp(b*( (x+0.5) - 80) ) )
}

gompertz_get_grad <- function(params, data, fn, epsilon = 1e-8) {
  a <- params[1]
  b <- params[2]
  
  x <- data$Age
  Dx <- data$Deaths
  Ex <- data$Exposure
  
  Q = b * (x + 0.5 - 80)
  P = exp(Q)
  N = a * P
  mux = N
  M = log(mux)
  G = Dx * M - Ex * mux
  L = sum(G)
  
  dLdG = rep(1, length(G))
  dGdmux = Dx / mux - Ex
  dmuxdN = 1
  dNda = P
  dNdP = a
  dPdQ = exp(Q)
  dQdb = x + 0.5 - 80
  
  grad_a = sum( (dLdG * dGdmux * dmuxdN * dNda) )
  grad_b = sum( (dLdG * dGdmux * dmuxdN * dNdP * dPdQ * dQdb) )
  
  grad <- c(grad_a, grad_b) #* 10^-6 
  
  if (any(!is.finite(grad))) {
    print(params)
    print(grad)
    stop("Non-finite gradient value encountered")
  }
  
  return(grad)
}

gompertz_get_mx <- function(a,b) {
  
  x <- 80:115
  Mx <- a*exp(b*( (x+0.5) - 80) ) 
  Mx <- cbind(x, Mx)
  
  return(Mx)
}

fit_gompertz <- function(init_param_fn, log_lik_fn, grad_fn, fn, data) {
  initial_params <- gompertz_get_init_params(log_lik_fn = log_lik_fn, data = data, fn = fn)
  optim_result <- optim(par = initial_params, 
                        fn = function(params) log_lik_fn(params, data, fn),
                        gr = function(params) grad_fn(params, data, fn), 
                        method = "L-BFGS-B", 
                        lower = c(.Machine$double.eps, .Machine$double.eps), 
                        upper = c(5, 5),
                        control = list(fnscale=-1))
  # high_age <- case_when(
  #   data$Country == "fra" ~ fra_high_age,
  #   data$Country != "fra" ~ scand_high_age
  # )
  # high_age <- high_age[1]
  opt_params <- optim_result$par
  fitted_mx <- gompertz_get_mx(optim_result$par[1], optim_result$par[2])
  
  log_lik_opt <- log_lik_fn(opt_params, data, fn)
  num_params <- length(opt_params)
  
  n <- nrow(data)
  aic <- 2 * num_params - 2 * log_lik_opt
  bic <- log(n) * num_params - 2 * log_lik_opt
  # data <- data |>
  #   mutate(observed_rates = Deaths/Exposure)
  # observed <- data$observed_rates
  # fitted <- fitted_mx[, 2]
  # sse <- sum((observed - fitted)^2)
  
  results <- list(fitted_mx = fitted_mx, AIC = aic, BIC = bic)
  
  return(results)
}

apply_gompertz_model <- function(data) {
  #fit_beard(
  result <- fit_gompertz(
    init_param_fn = gompertz_get_init_params, 
    log_lik_fn = log_lik_fn, 
    grad_fn = gompertz_get_grad, 
    fn = gompertz_fn, 
    data = data
  )
  
  return(result)
}

# fit_all_beard <- function(data) {
#   data |>
#     group_by(Cohort, Country, Sex) |>
#     do(fitted_rates = apply_beard_model(.)) |>
#     ungroup()
# }

fit_all_gompertz <- function(data) {
  data |>
    group_by(Cohort, Country, Sex) |>
    do({
      model_result <- apply_gompertz_model(.)
      tibble(
        fitted_rates = list(model_result$fitted_mx),
        AIC = model_result$AIC,
        BIC = model_result$BIC
      )
    }) |>
    ungroup()
}



#### Makeham ####
makeham_get_init_params <- function(log_lik_fn, data, fn,
                                    a_values=seq(0.0001, 1, length.out = 20), 
                                    b_values=seq(0.0001, 1, length.out = 20),
                                    g_values=seq(-1, 1, length.out = 20)) {
  # Perform grid search
  grid_search_results <- expand.grid(a = a_values, b = b_values, g = g_values)
  
  # Evaluate log-likelihood for each combination of a and b
  grid_search_results$log_lik <- apply(grid_search_results, 1, function(row) {
    params <- c(row['a'], row['b'], row['g'])
    return(log_lik_fn(params, data, fn))
  })
  
  # Find the combination with the highest log-likelihood (least negative)
  best_params <- grid_search_results[which.max(grid_search_results$log_lik), ]
  
  # Extract the best initial values for a and b
  initial_params <- c(a = best_params$a, b = best_params$b, g = best_params$g)
  
  return(initial_params)
}

makeham_fn <- function(params, data) {
  a <- params[1]
  b <- params[2]
  g <- params[3]
  x <- data$Age
  
  ( g + (a*exp(b*( (x+0.5) - 80) )) ) 
}

makeham_get_grad <- function(params, data, fn, epsilon = 1e-8) {
  a <- params[1]
  b <- params[2]
  g <- params[3]
  
  x <- data$Age
  Dx <- data$Deaths
  Ex <- data$Exposure
  
  Q = b * (x + 0.5 - 80)
  P = exp(Q)
  N = a * P
  mux = g + N
  M = log(mux)
  G = Dx * M - Ex * mux
  L = sum(G)
  
  dLdG = rep(1, length(G))
  dGdmux = Dx / mux - Ex
  dmuxdN = 1 
  dmuxdg = 1
  dNda = P
  dNdP = a
  dPdQ = exp(Q)
  dQdb = x + 0.5 - 80
  
  grad_a = sum( (dLdG * dGdmux * dmuxdN * dNda) )
  grad_b = sum( (dLdG * dGdmux * dmuxdN * dNdP * dPdQ * dQdb) )
  grad_g = sum( (dLdG * dGdmux * dmuxdg) )
  
  grad <- c(grad_a, grad_b, grad_g) #* 10^-6 
  
  if (any(!is.finite(grad))) {
    print(params)
    print(grad)
    stop("Non-finite gradient value encountered")
  }
  
  return(grad)
}

makeham_get_mx <- function(a,b,g) {
  
  x <- 80:115
  Mx <- ( g + (a*exp(b*( (x+0.5) - 80) )) )
  Mx <- cbind(x, Mx)
  
  return(Mx)
}

fit_makeham <- function(init_param_fn, log_lik_fn, grad_fn, fn, data) {
  initial_params <- makeham_get_init_params(log_lik_fn = log_lik_fn, data = data, fn = fn)
  optim_result <- optim(par = initial_params, 
                        fn = function(params) log_lik_fn(params, data, fn),
                        gr = function(params) grad_fn(params, data, fn), 
                        method = "L-BFGS-B", 
                        lower = c(.Machine$double.eps, .Machine$double.eps), 
                        upper = c(5, 5),
                        control = list(fnscale=-1))
  # high_age <- case_when(
  #   data$Country == "fra" ~ fra_high_age,
  #   data$Country != "fra" ~ scand_high_age
  # )
  # high_age <- high_age[1]

  opt_params <- optim_result$par
  fitted_mx <- makeham_get_mx(optim_result$par[1], optim_result$par[2], optim_result$par[3])
  
  log_lik_opt <- log_lik_fn(opt_params, data, fn)
  num_params <- length(opt_params)
  
  n <- nrow(data)
  aic <- 2 * num_params - 2 * log_lik_opt
  bic <- log(n) * num_params - 2 * log_lik_opt
  # data <- data |>
  #   mutate(observed_rates = Deaths/Exposure)
  # observed <- data$observed_rates
  # fitted <- fitted_mx[, 2]
  # sse <- sum((observed - fitted)^2)
  
  results <- list(fitted_mx = fitted_mx, AIC = aic, BIC = bic)
  
  return(results)
}

apply_makeham_model <- function(data) {
  #fit_beard(
  result <- fit_makeham(
    init_param_fn = makeham_get_init_params, 
    log_lik_fn = log_lik_fn, 
    grad_fn = makeham_get_grad, 
    fn = makeham_fn, 
    data = data
  )
  
  return(result)
}

# fit_all_beard <- function(data) {
#   data |>
#     group_by(Cohort, Country, Sex) |>
#     do(fitted_rates = apply_beard_model(.)) |>
#     ungroup()
# }

fit_all_makeham <- function(data) {
  data |>
    group_by(Cohort, Country, Sex) |>
    do({
      model_result <- apply_makeham_model(.)
      tibble(
        fitted_rates = list(model_result$fitted_mx),
        AIC = model_result$AIC,
        BIC = model_result$BIC
      )
    }) |>
    ungroup()
}




#### Log-Quadratic ####
lq_get_init_params <- function(log_lik_fn, data, fn,
                                    a_values=seq(-10, 3, length.out = 20), 
                                    b_values=seq(-5, 10, length.out = 20),
                                    g_values=seq(-5, 1, length.out = 20)) {
  # Perform grid search
  grid_search_results <- expand.grid(a = a_values, b = b_values, g = g_values)
  
  # Evaluate log-likelihood for each combination of a and b
  grid_search_results$log_lik <- apply(grid_search_results, 1, function(row) {
    params <- c(row['a'], row['b'], row['g'])
    return(log_lik_fn(params, data, fn))
  })
  
  # Find the combination with the highest log-likelihood (least negative)
  best_params <- grid_search_results[which.max(grid_search_results$log_lik), ]
  
  # Extract the best initial values for a and b
  initial_params <- c(a = best_params$a, b = best_params$b, g = best_params$g)
  
  return(initial_params)
}

lq_fn <- function(params, data) {
  a <- params[1]
  b <- params[2]
  g <- params[3]
  x <- data$Age
  
  ( exp(a + b*( (x+0.5) - 80) + g*( (x+0.5) - 80)^2) ) 
}

lq_get_grad <- function(params, data, fn, epsilon = 1e-8) {
  a <- params[1]
  b <- params[2]
  g <- params[3]
  
  x <- data$Age
  Dx <- data$Deaths
  Ex <- data$Exposure
  
  Q = b * (x + 0.5 - 80)
  R = g * (x + 0.5 - 80)^2
  S = a + Q + R
  mux = exp(S)
  M = log(mux)
  G = Dx * M - Ex * mux
  L = sum(G)
  
  dLdG = rep(1, length(G))
  dGdmux = Dx / mux - Ex
  dmuxdS = exp(S) 
  dSda = 1
  dSdQ = 1
  dSdR = 1
  dRdg = (x + 0.5 - 80)^2
  dQdb = x + 0.5 - 80
  
  grad_a = sum( (dLdG * dGdmux * dmuxdS * dSda) )
  grad_b = sum( (dLdG * dGdmux * dmuxdS * dSdQ * dQdb) )
  grad_g = sum( (dLdG * dGdmux * dmuxdS * dSdR * dRdg) )
  
  grad <- c(grad_a, grad_b, grad_g) #* 10^-6 
  
  if (any(!is.finite(grad))) {
    print(params)
    print(grad)
    stop("Non-finite gradient value encountered")
  }
  
  return(grad)
}

lq_get_mx <- function(a,b,g) {
  
  x <- 80:115
  Mx <- ( exp(a + b*( (x+0.5) - 80) + g*( (x+0.5) - 80)^2) )
  Mx <- cbind(x, Mx)
  
  return(Mx)
}

fit_lq <- function(init_param_fn, log_lik_fn, grad_fn, fn, data) {
  initial_params <- init_param_fn(log_lik_fn = log_lik_fn, data = data, fn = fn)
  optim_result <- optim(par = initial_params, 
                        fn = function(params) log_lik_fn(params, data, fn),
                        gr = function(params) grad_fn(params, data, fn), 
                        method = "BFGS",
                        control = list(fnscale=-1))
  # high_age <- case_when(
  #   data$Country == "fra" ~ fra_high_age,
  #   data$Country != "fra" ~ scand_high_age
  # )
  # high_age <- high_age[1]

  opt_params <- optim_result$par
  fitted_mx <- lq_get_mx(optim_result$par[1], optim_result$par[2], optim_result$par[3])
  
  log_lik_opt <- log_lik_fn(opt_params, data, fn)
  num_params <- length(opt_params)
  
  n <- nrow(data)
  aic <- 2 * num_params - 2 * log_lik_opt
  bic <- log(n) * num_params - 2 * log_lik_opt
  # data <- data |>
  #   mutate(observed_rates = Deaths/Exposure)
  # observed <- data$observed_rates
  # fitted <- fitted_mx[, 2]
  # sse <- sum((observed - fitted)^2)
  
  results <- list(fitted_mx = fitted_mx, AIC = aic, BIC = bic)
  
  return(results)
}

apply_lq_model <- function(data) {
  result <- fit_lq(
    init_param_fn = lq_get_init_params, 
    log_lik_fn = log_lik_fn, 
    grad_fn = lq_get_grad, 
    fn = lq_fn, 
    data = data
  )
  
  return(result)
}

# fit_all_beard <- function(data) {
#   data |>
#     group_by(Cohort, Country, Sex) |>
#     do(fitted_rates = apply_beard_model(.)) |>
#     ungroup()
# }

fit_all_lq <- function(data) {
  data |>
    group_by(Cohort, Country, Sex) |>
    do({
      model_result <- apply_lq_model(.)
      tibble(
        fitted_rates = list(model_result$fitted_mx),
        AIC = model_result$AIC,
        BIC = model_result$BIC
      )
    }) |>
    ungroup()
}
