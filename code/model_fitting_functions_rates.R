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
  
  # if (!is.finite(log_lik)) {
  #   print(params)
  #   print(mux)
  #   stop("Non-finite log likelihood value encountered")
  # }
  
  return(-log_lik * 10^-6) # scale down according to HMD procedure 
}

##########################
# Model specific functions 
##########################

#### Kannisto ####
get_initial_params <- function(log_lik_fn, data, fn,
                               a_values=seq(0.01, .3, length.out = 10), 
                               b_values=seq(0.01, .3, length.out = 10)) {
  # Perform grid search
  grid_search_results <- expand.grid(a = a_values, b = b_values)
  
  # Evaluate log-likelihood for each combination of a and b
  grid_search_results$log_lik <- apply(grid_search_results, 1, function(row) {
    params <- c(row['a'], row['b'])
    return(log_lik_fn(params, data, fn))
  })
  
  # Find the combination with the highest log-likelihood (least negative)
  best_params <- grid_search_results[which.min(grid_search_results$log_lik), ]
  
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

# manual gradient function 
# kannisto_get_grad <- function(params, data, fn, epsilon = 1e-8) {
#   a <- params[1]
#   b <- params[2]
# 
#   x <- data$Age
#   Dx <- data$Deaths
#   Ex <- data$Exposure
#   g <- a * exp(b * (x + 0.5 - 80))
# 
#   grad_a <- sum((Dx * exp(b * (x + 0.5 - 80))) / (g * (1 + g)) -
#                   (Ex * exp(b * (x + 0.5 - 80))) / ((1 + g)^2))
#   grad_b <- sum((Dx * a * exp(b * (x + 0.5 - 80)) * (x + 0.5 - 80)) /
#                   (g * (1 + g)) -
#                   (Ex * a * exp(b * (x + 0.5 - 80)) * (x + 0.5 - 80)) / ((1 + g)^2))
# 
#   grad <- -c(grad_a, grad_b) * 10^-6
# 
#   if (any(!is.finite(grad))) {
#     print(params)
#     print(grad)
#     stop("Non-finite gradient value encountered")
#   }
# 
#   return(grad) # Return negative gradient
# }

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
  
  grad <- -c(grad_a, grad_b) * 10^-6 
  
  if (any(!is.finite(grad))) {
       print(params)
       print(grad)
       stop("Non-finite gradient value encountered")
     }

  return(grad)
}

kannisto_get_mx <- function(a,b) {
  
  x <- 80:110
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
                        upper = c(5, 5))
  fitted_mx <- kannisto_get_mx(optim_result$par[1], optim_result$par[2])
  
  return(fitted_mx)
}

apply_kannisto_model <- function(data) {
  fit_kannisto(
    init_param_fn = get_initial_params, 
    log_lik_fn = log_lik_fn, 
    grad_fn = kannisto_get_grad, 
    fn = kannisto_fn, 
    data = data
  )
}

fit_all_kannisto <- function(data) {
  data |>
    group_by(Year, country, Sex) |>
    do(fitted_rates = apply_kannisto_model(.)) |>
    ungroup()
}

#### Beard ####
beard_get_init_params <- function(log_lik_fn, data, fn,
                                  a_values=seq(0.001, 1, length.out = 20), 
                                  b_values=seq(0.001, 1, length.out = 20),
                                  d_values=seq(0.001, 1, length.out = 20)) {
  # Perform grid search
  grid_search_results <- expand.grid(a = a_values, b = b_values, d = d_values)
  
  # Evaluate log-likelihood for each combination of a and b
  grid_search_results$log_lik <- apply(grid_search_results, 1, function(row) {
    params <- c(row['a'], row['b'], row['d'])
    return(log_lik_fn(params, data, fn))
  })
  
  # Find the combination with the highest log-likelihood (least negative)
  best_params <- grid_search_results[which.min(grid_search_results$log_lik), ]
  
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
  
  grad <- -c(grad_a, grad_b, grad_d) * 10^-6 
  
  if (any(!is.finite(grad))) {
    print(params)
    print(grad)
    stop("Non-finite gradient value encountered")
  }
  
  return(grad)
}

beard_get_mx <- function(a,b,d) {
  
  x <- 80:110
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
                        upper = c(5, 5))
  fitted_mx <- beard_get_mx(optim_result$par[1], optim_result$par[2], optim_result$par[3])
  
  return(fitted_mx)
}

apply_beard_model <- function(data) {
  fit_beard(
    init_param_fn = beard_get_init_params, 
    log_lik_fn = log_lik_fn, 
    grad_fn = beard_get_grad, 
    fn = beard_fn, 
    data = data
  )
}

fit_all_beard <- function(data) {
  data |>
    group_by(Year, country, Sex) |>
    do(fitted_rates = apply_beard_model(.)) |>
    ungroup()
}






