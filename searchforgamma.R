objective_function <- function(params) {
  alpha <- params[1]
  beta <- params[2]
  
  # Calculate mean and quantiles
  mean_calculated <- alpha / beta
  q25_calculated <- qgamma(0.25, shape = alpha, rate = beta)
  q75_calculated <- qgamma(0.75, shape = alpha, rate = beta)
  
  # Calculate differences from targets
  diff_mean <- (mean_calculated - 1.28)^2
  diff_q25 <- (q25_calculated - 1.19)^2
  diff_q75 <- (q75_calculated - 1.37)^2
  
  # Return sum of squared differences
  return(diff_mean + diff_q25 + diff_q75)
}

# Initial guesses for alpha and beta
initial_guess <- c(1, 0.5)

# Use optim to find the parameters that minimize the objective function
optim_result <- optim(initial_guess, objective_function, method = "L-BFGS-B", lower = c(0.001, 0.001))

# Print the optimized parameters
optim_result$par
