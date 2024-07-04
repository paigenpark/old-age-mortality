##########################
# Generic functions 
##########################

## TODO
# This file is currently identical to the model_fitting_functions_rates.R file.
# Needs to be updated with logliklihood for counts, but I'm not sure how to set 
# that up yet (see pages 31 and 32 of Feehan 2018)

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

# Function to plot original and fitted rates
plot_original_vs_fitted <- function(original_df, smoothed_df, year, country, sex) {
  # Filter the original dataframe
  original_data <- original_df %>%
    filter(Year == year, country == country, Sex == sex)
  
  # Extract the fitted rates for the same group
  fitted_data <- smoothed_df %>%
    filter(Year == year, country == country, Sex == sex) %>%
    select(fitted_rates) %>%
    hoist(fitted_rates, c("x", "Mx"))
  
  fitted_data <- fitted_data %>%
    rename(Age = x, FittedRate = Mx)
  
  # Combine the original and fitted data
  combined_data <- original_data %>%
    left_join(fitted_data, by = c("Age"))
  
  # Plot using ggplot2
  ggplot(combined_data, aes(x = Age)) +
    geom_line(aes(y = Deaths / Exposure, color = "Original Rate")) +
    geom_line(aes(y = FittedRate, color = "Fitted Rate")) +
    labs(
      title = paste("Original vs Fitted Mortality Rates for", country, "in", year, "(", sex, ")"),
      x = "Age",
      y = "Mortality Rate"
    ) +
    scale_color_manual(name = "Rate Type", values = c("Original Rate" = "blue", "Fitted Rate" = "red")) +
    theme_minimal()
}

##########################
# Model specific functions 
##########################

#### Kannisto ####
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
  g <- a * exp(b * (x + 0.5 - 80))
  mux <- g / (1 + g)
  
  grad_a <- sum((Dx * exp(b * (x + 0.5 - 80))) / (g * (1 + g)) - 
                  (Ex * exp(b * (x + 0.5 - 80))) / ((1 + g)^2))
  grad_b <- sum((Dx * a * exp(b * (x + 0.5 - 80)) * (x + 0.5 - 80)) / 
                  (g * (1 + g)) - 
                  (Ex * a * exp(b * (x + 0.5 - 80)) * (x + 0.5 - 80)) / ((1 + g)^2))
  
  grad <- -c(grad_a, grad_b) * 10^-6
  
  if (any(!is.finite(grad))) {
    print(params)
    print(mux_plus_0.5)
    print(grad)
    stop("Non-finite gradient value encountered")
  }
  
  return(-c(grad_a, grad_b) * 10^-6) # Return negative gradient
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





