# Function to check stationarity for a single set of phi1, phi2
check_stationarity <- function(phi1, phi2) {
  # Construct the polynomial (replace if you get a different polynomial)
  polynomial <- c(-0.7 , -0.2, 1) 

  # Solve for the roots
  roots <- polyroot(polynomial)

  # Print the roots
  print(paste0("Roots: ", roots))

  # Check if roots are within the unit circle
  all_within_unit_circle <- all(abs(roots) > 1)

  # Print the result
  if (all_within_unit_circle) {
    print("Stationary")
  } else {
    print("Not Stationary")
  }
}




































phi1 <- -0.7
phi2 <- -0.2

# Characteristic equation (using polynomial coefficients)
coeff <- c(1, phi1, phi2) # Note: the signs of phi1 and phi2 are inverted in the characteristic equation

# Solve for the roots of the polynomial
roots <- polyroot(coeff)

# Check if roots are inside the unit circle (stationarity condition)
stationary <- all(abs(roots) > 1)

# Print the results
print(roots)
print(paste0("Stationary: ", stationary))

# Set seed for reproducibility
set.seed(123)  # Feel free to change the seed number

# Number of realizations
n_realizations <- 5

# Number of observations
n <- 200

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations) 
for (i in 1:n_realizations) {
  # Temporary storage for a single realization
  x <- numeric(n)

  # Generate the first two values (adjust if you have specific starting values)
  x[1] <- rnorm(1) 
  x[2] <- phi1 * x[1] + rnorm(1) 

  # Generate the rest of the series
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }

  simulations[, i] <- x
}

# Calculate variances of each realization
realization_variances <- apply(simulations, 2, var) 

# Print the variances
print(realization_variances)

# Plot the realizations
matplot(simulations, type = "l", col = 1:n_realizations, 
        main = "Realizations of AR(2) Process", xlab = "Time", ylab = "Value")


# Number of lags for ACF
nlag <- 30

# Function to calculate theoretical ACF (Yule-Walker implementation)
calculate_theoretical_acf <- function(phi1, phi2, nlag) {
  if (nlag < 0) {
    stop("nlag must be non-negative")
  }

  acf_values <- numeric(nlag + 1)

  acf_values[1] <- 1 
  if (nlag >= 1) {
    acf_values[2] <- phi1 / (1 - phi2) 
  }

  for (k in 3:(nlag + 1)) {
    acf_values[k] <- phi1 * acf_values[k - 1] + phi2 * acf_values[k - 2]
  }

  return(acf_values)
}

# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)


# Plot the theoretical ACF
plot(theoretical_acf, type = "h", lwd = 2, col = "blue",
     main = "Theoretical ACF of AR(2) Process", xlab = "Lag", ylab = "ACF")

# Add reference line at zero
abline(h = 0, col = "gray")


# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)

# Function to calculate empirical ACF
empirical_acf <- function(x, nlag) {
  acf(x, lag.max = nlag, type = "correlation", plot = FALSE)$acf 
}

# Calculate empirical ACFs for each realization
empirical_acfs <- apply(simulations, 2, empirical_acf, nlag = nlag)


# Calculate variances of the ACF at each lag
acf_variances <- apply(empirical_acfs, 1, var)  # Note: empirical_acfs might need to be transposed depending on its structure

# Print the variances of the ACF
print(acf_variances)

# Optionally, plot the variances of the ACF
plot(acf_variances, type = "o", pch = 20, col = "blue", xlab = "Lag", ylab = "Variance of ACF",
     main = "Variance of ACF at Each Lag")



# Determine layout for multiple plots
n_plots <- n_realizations  
plot_rows <- ceiling(sqrt(n_plots)) # Adjust if you want a different layout
plot_cols <- ceiling(n_plots / plot_rows) 

# Set up the plot with correct layout
par(mfrow = c(plot_rows, plot_cols)) 

# Confidence level for the ACF confidence intervals (95% confidence)
z <- 1.96
conf_int <- z / sqrt(n)

# Plot each ACF individually with confidence intervals
for (i in 1:n_realizations) {
  plot(empirical_acfs[, i], type = "h", lwd = 2, 
       main = paste0("Realization ", i), xlab = "Lag", ylab = "ACF")

  # Add theoretical ACF (with a different color)
  lines(theoretical_acf, col = "blue", lwd = 2)

  abline(h = 0, col = "gray") # Zero line for reference
  
  # Add confidence intervals
  abline(h = conf_int, col = "red", lty = 2) # Upper confidence interval
  abline(h = -conf_int, col = "red", lty = 2) # Lower confidence interval
}

# Reset to default single plotting area
par(mfrow = c(1, 1))


set.seed(123)

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations)
for (i in 1:n_realizations) {
  x <- numeric(n)
  x[1:2] <- rnorm(2)
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }
  simulations[, i] <- x
}

plot(NULL, xlim=c(1, nlag), ylim=c(-1, 1), xlab="Lag", ylab="ACF",
     main="Empirical ACF of AR(2) Realizations")
abline(h=0, col="gray")
colors <- rainbow(n_realizations)

for (i in 1:n_realizations) {
  # Ensure that you get 30 lags exactly
  acf_values <- acf(simulations[, i], lag.max=nlag, plot = FALSE)$acf[1:(nlag+1)]
  lines(1:nlag, acf_values[-1], type="o", col=colors[i], pch=20, lwd=2)
}

# Add confidence interval lines
conf_int <- qnorm((1 + 0.95) / 2) / sqrt(n)
abline(h = c(-conf_int, conf_int), col = "blue", lty = "dotted")

# Add a legend to the top right
legend("topright", legend=paste("Realization", 1:n_realizations), 
       col=colors, lty=1, lwd=2, cex=0.8)






































## 1.7 i

phi1 <- -0.2
phi2 <- 0.2

# Characteristic equation (using polynomial coefficients)
coeff <- c(1, -phi1, -phi2) # Note: the signs of phi1 and phi2 are inverted in the characteristic equation

# Solve for the roots of the polynomial
roots <- polyroot(coeff)

# Check if roots are inside the unit circle (stationarity condition)
stationary <- all(abs(roots) > 1)

# Print the results
print(roots)
print(paste0("Stationary: ", stationary)) 

# Set seed for reproducibility
set.seed(123)  # Feel free to change the seed number

# Number of realizations
n_realizations <- 5

# Number of observations
n <- 200

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations) 
for (i in 1:n_realizations) {
  # Temporary storage for a single realization
  x <- numeric(n)

  # Generate the first two values (adjust if you have specific starting values)
  x[1] <- rnorm(1) 
  x[2] <- phi1 * x[1] + rnorm(1) 

  # Generate the rest of the series
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }

  simulations[, i] <- x
}

# Calculate variances of each realization
realization_variances <- apply(simulations, 2, var) 

# Print the variances
print(realization_variances)

# Plot the realizations
matplot(simulations, type = "l", col = 1:n_realizations, 
        main = "Realizations of AR(2) Process", xlab = "Time", ylab = "Value")


# Number of lags for ACF
nlag <- 30

# Function to calculate theoretical ACF (Yule-Walker implementation)
calculate_theoretical_acf <- function(phi1, phi2, nlag) {
  if (nlag < 0) {
    stop("nlag must be non-negative")
  }

  acf_values <- numeric(nlag + 1)

  acf_values[1] <- 1 
  if (nlag >= 1) {
    acf_values[2] <- phi1 / (1 - phi2) 
  }

  for (k in 3:(nlag + 1)) {
    acf_values[k] <- phi1 * acf_values[k - 1] + phi2 * acf_values[k - 2]
  }

  return(acf_values)
}

# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)


# Plot the theoretical ACF
plot(theoretical_acf, type = "h", lwd = 2, col = "blue",
     main = "Theoretical ACF of AR(2) Process", xlab = "Lag", ylab = "ACF")

# Add reference line at zero
abline(h = 0, col = "gray")


# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)

# Function to calculate empirical ACF
empirical_acf <- function(x, nlag) {
  acf(x, lag.max = nlag, type = "correlation", plot = FALSE)$acf 
}

# Calculate empirical ACFs for each realization
empirical_acfs <- apply(simulations, 2, empirical_acf, nlag = nlag)


# Calculate variances of the ACF at each lag
acf_variances <- apply(empirical_acfs, 1, var)  # Note: empirical_acfs might need to be transposed depending on its structure

# Print the variances of the ACF
print(acf_variances)

# Optionally, plot the variances of the ACF
plot(acf_variances, type = "o", pch = 20, col = "blue", xlab = "Lag", ylab = "Variance of ACF",
     main = "Variance of ACF at Each Lag")



# Determine layout for multiple plots
n_plots <- n_realizations  
plot_rows <- ceiling(sqrt(n_plots)) # Adjust if you want a different layout
plot_cols <- ceiling(n_plots / plot_rows) 

# Set up the plot with correct layout
par(mfrow = c(plot_rows, plot_cols)) 

# Confidence level for the ACF confidence intervals (95% confidence)
z <- 1.96
conf_int <- z / sqrt(n)

# Plot each ACF individually with confidence intervals
for (i in 1:n_realizations) {
  plot(empirical_acfs[, i], type = "h", lwd = 2, 
       main = paste0("Realization ", i), xlab = "Lag", ylab = "ACF")

  # Add theoretical ACF (with a different color)
  lines(theoretical_acf, col = "blue", lwd = 2)

  abline(h = 0, col = "gray") # Zero line for reference
  
  # Add confidence intervals
  abline(h = conf_int, col = "red", lty = 2) # Upper confidence interval
  abline(h = -conf_int, col = "red", lty = 2) # Lower confidence interval
}

# Reset to default single plotting area
par(mfrow = c(1, 1))


set.seed(123)

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations)
for (i in 1:n_realizations) {
  x <- numeric(n)
  x[1:2] <- rnorm(2)
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }
  simulations[, i] <- x
}

plot(NULL, xlim=c(1, nlag), ylim=c(-1, 1), xlab="Lag", ylab="ACF",
     main="Empirical ACF of AR(2) Realizations")
abline(h=0, col="gray")
colors <- rainbow(n_realizations)

for (i in 1:n_realizations) {
  # Ensure that you get 30 lags exactly
  acf_values <- acf(simulations[, i], lag.max=nlag, plot = FALSE)$acf[1:(nlag+1)]
  lines(1:nlag, acf_values[-1], type="o", col=colors[i], pch=20, lwd=2)
}

# Add confidence interval lines
conf_int <- qnorm((1 + 0.95) / 2) / sqrt(n)
abline(h = c(-conf_int, conf_int), col = "blue", lty = "dotted")

# Add a legend to the top right
legend("topright", legend=paste("Realization", 1:n_realizations), 
       col=colors, lty=1, lwd=2, cex=0.8)










## 1.7 ii

phi1 <- 0.7
phi2 <- 0.2

# Characteristic equation (using polynomial coefficients)
coeff <- c(1, -phi1, -phi2) # Note: the signs of phi1 and phi2 are inverted in the characteristic equation

# Solve for the roots of the polynomial
roots <- polyroot(coeff)

# Check if roots are inside the unit circle (stationarity condition)
stationary <- all(abs(roots) < 1)

# Print the results
print(roots)
print(paste0("Stationary: ", stationary))

# Set seed for reproducibility
set.seed(123)  # Feel free to change the seed number

# Number of realizations
n_realizations <- 5

# Number of observations
n <- 200

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations) 
for (i in 1:n_realizations) {
  # Temporary storage for a single realization
  x <- numeric(n)

  # Generate the first two values (adjust if you have specific starting values)
  x[1] <- rnorm(1) 
  x[2] <- phi1 * x[1] + rnorm(1) 

  # Generate the rest of the series
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }

  simulations[, i] <- x
}

# Calculate variances of each realization
realization_variances <- apply(simulations, 2, var) 

# Print the variances
print(realization_variances)

# Plot the realizations
matplot(simulations, type = "l", col = 1:n_realizations, 
        main = "Realizations of AR(2) Process", xlab = "Time", ylab = "Value")


# Number of lags for ACF
nlag <- 30

# Function to calculate theoretical ACF (Yule-Walker implementation)
calculate_theoretical_acf <- function(phi1, phi2, nlag) {
  if (nlag < 0) {
    stop("nlag must be non-negative")
  }

  acf_values <- numeric(nlag + 1)

  acf_values[1] <- 1 
  if (nlag >= 1) {
    acf_values[2] <- phi1 / (1 - phi2) 
  }

  for (k in 3:(nlag + 1)) {
    acf_values[k] <- phi1 * acf_values[k - 1] + phi2 * acf_values[k - 2]
  }

  return(acf_values)
}

# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)


# Plot the theoretical ACF
plot(theoretical_acf, type = "h", lwd = 2, col = "blue",
     main = "Theoretical ACF of AR(2) Process", xlab = "Lag", ylab = "ACF")

# Add reference line at zero
abline(h = 0, col = "gray")


# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)

# Function to calculate empirical ACF
empirical_acf <- function(x, nlag) {
  acf(x, lag.max = nlag, type = "correlation", plot = FALSE)$acf 
}

# Calculate empirical ACFs for each realization
empirical_acfs <- apply(simulations, 2, empirical_acf, nlag = nlag)


# Calculate variances of the ACF at each lag
acf_variances <- apply(empirical_acfs, 1, var)  # Note: empirical_acfs might need to be transposed depending on its structure

# Print the variances of the ACF
print(acf_variances)

# Optionally, plot the variances of the ACF
plot(acf_variances, type = "o", pch = 20, col = "blue", xlab = "Lag", ylab = "Variance of ACF",
     main = "Variance of ACF at Each Lag")



# Determine layout for multiple plots
n_plots <- n_realizations  
plot_rows <- ceiling(sqrt(n_plots)) # Adjust if you want a different layout
plot_cols <- ceiling(n_plots / plot_rows) 

# Set up the plot with correct layout
par(mfrow = c(plot_rows, plot_cols)) 

# Confidence level for the ACF confidence intervals (95% confidence)
z <- 1.96
conf_int <- z / sqrt(n)

# Plot each ACF individually with confidence intervals
for (i in 1:n_realizations) {
  plot(empirical_acfs[, i], type = "h", lwd = 2, 
       main = paste0("Realization ", i), xlab = "Lag", ylab = "ACF")

  # Add theoretical ACF (with a different color)
  lines(theoretical_acf, col = "blue", lwd = 2)

  abline(h = 0, col = "gray") # Zero line for reference
  
  # Add confidence intervals
  abline(h = conf_int, col = "red", lty = 2) # Upper confidence interval
  abline(h = -conf_int, col = "red", lty = 2) # Lower confidence interval
}

# Reset to default single plotting area
par(mfrow = c(1, 1))


set.seed(123)

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations)
for (i in 1:n_realizations) {
  x <- numeric(n)
  x[1:2] <- rnorm(2)
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }
  simulations[, i] <- x
}

plot(NULL, xlim=c(1, nlag), ylim=c(-1, 1), xlab="Lag", ylab="ACF",
     main="Empirical ACF of AR(2) Realizations")
abline(h=0, col="gray")
colors <- rainbow(n_realizations)

for (i in 1:n_realizations) {
  # Ensure that you get 30 lags exactly
  acf_values <- acf(simulations[, i], lag.max=nlag, plot = FALSE)$acf[1:(nlag+1)]
  lines(1:nlag, acf_values[-1], type="o", col=colors[i], pch=20, lwd=2)
}

# Add confidence interval lines
conf_int <- qnorm((1 + 0.95) / 2) / sqrt(n)
abline(h = c(-conf_int, conf_int), col = "blue", lty = "dotted")

# Add a legend to the top right
legend("topright", legend=paste("Realization", 1:n_realizations), 
       col=colors, lty=1, lwd=2, cex=0.8)











## 1.7 iii

phi1 <- -0.8
phi2 <- 0.2

# Characteristic equation (using polynomial coefficients)
coeff <- c(1, -phi1, -phi2) # Note: the signs of phi1 and phi2 are inverted in the characteristic equation

# Solve for the roots of the polynomial
roots <- polyroot(coeff)

# Check if roots are inside the unit circle (stationarity condition)
stationary <- all(abs(roots) > 1)

# Print the results
print(roots)
print(paste0("Stationary: ", stationary))

# Set seed for reproducibility
set.seed(123)  # Feel free to change the seed number

# Number of realizations
n_realizations <- 5

# Number of observations
n <- 200

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations) 
for (i in 1:n_realizations) {
  # Temporary storage for a single realization
  x <- numeric(n)

  # Generate the first two values (adjust if you have specific starting values)
  x[1] <- rnorm(1) 
  x[2] <- phi1 * x[1] + rnorm(1) 

  # Generate the rest of the series
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }

  simulations[, i] <- x
}

# Calculate variances of each realization
realization_variances <- apply(simulations, 2, var) 

# Print the variances
print(realization_variances)

# Plot the realizations
matplot(simulations, type = "l", col = 1:n_realizations, 
        main = "Realizations of AR(2) Process", xlab = "Time", ylab = "Value")


# Number of lags for ACF
nlag <- 30

# Function to calculate theoretical ACF (Yule-Walker implementation)
calculate_theoretical_acf <- function(phi1, phi2, nlag) {
  if (nlag < 0) {
    stop("nlag must be non-negative")
  }

  acf_values <- numeric(nlag + 1)

  acf_values[1] <- 1 
  if (nlag >= 1) {
    acf_values[2] <- phi1 / (1 - phi2) 
  }

  for (k in 3:(nlag + 1)) {
    acf_values[k] <- phi1 * acf_values[k - 1] + phi2 * acf_values[k - 2]
  }

  return(acf_values)
}

# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)


# Plot the theoretical ACF
plot(theoretical_acf, type = "h", lwd = 2, col = "blue",
     main = "Theoretical ACF of AR(2) Process", xlab = "Lag", ylab = "ACF")

# Add reference line at zero
abline(h = 0, col = "gray")


# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)

# Function to calculate empirical ACF
empirical_acf <- function(x, nlag) {
  acf(x, lag.max = nlag, type = "correlation", plot = FALSE)$acf 
}

# Calculate empirical ACFs for each realization
empirical_acfs <- apply(simulations, 2, empirical_acf, nlag = nlag)


# Calculate variances of the ACF at each lag
acf_variances <- apply(empirical_acfs, 1, var)  # Note: empirical_acfs might need to be transposed depending on its structure

# Print the variances of the ACF
print(acf_variances)

# Optionally, plot the variances of the ACF
plot(acf_variances, type = "o", pch = 20, col = "blue", xlab = "Lag", ylab = "Variance of ACF",
     main = "Variance of ACF at Each Lag")



# Determine layout for multiple plots
n_plots <- n_realizations  
plot_rows <- ceiling(sqrt(n_plots)) # Adjust if you want a different layout
plot_cols <- ceiling(n_plots / plot_rows) 

# Set up the plot with correct layout
par(mfrow = c(plot_rows, plot_cols)) 

# Confidence level for the ACF confidence intervals (95% confidence)
z <- 1.96
conf_int <- z / sqrt(n)

# Plot each ACF individually with confidence intervals
for (i in 1:n_realizations) {
  plot(empirical_acfs[, i], type = "h", lwd = 2, 
       main = paste0("Realization ", i), xlab = "Lag", ylab = "ACF")

  # Add theoretical ACF (with a different color)
  lines(theoretical_acf, col = "blue", lwd = 2)

  abline(h = 0, col = "gray") # Zero line for reference
  
  # Add confidence intervals
  abline(h = conf_int, col = "red", lty = 2) # Upper confidence interval
  abline(h = -conf_int, col = "red", lty = 2) # Lower confidence interval
}

# Reset to default single plotting area
par(mfrow = c(1, 1))


set.seed(123)

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations)
for (i in 1:n_realizations) {
  x <- numeric(n)
  x[1:2] <- rnorm(2)
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }
  simulations[, i] <- x
}

plot(NULL, xlim=c(1, nlag), ylim=c(-1, 1), xlab="Lag", ylab="ACF",
     main="Empirical ACF of AR(2) Realizations")
abline(h=0, col="gray")
colors <- rainbow(n_realizations)

for (i in 1:n_realizations) {
  # Ensure that you get 30 lags exactly
  acf_values <- acf(simulations[, i], lag.max=nlag, plot = FALSE)$acf[1:(nlag+1)]
  lines(1:nlag, acf_values[-1], type="o", col=colors[i], pch=20, lwd=2)
}

# Add confidence interval lines
conf_int <- qnorm((1 + 0.95) / 2) / sqrt(n)
abline(h = c(-conf_int, conf_int), col = "blue", lty = "dotted")

# Add a legend to the top right
legend("topright", legend=paste("Realization", 1:n_realizations), 
       col=colors, lty=1, lwd=2, cex=0.8)





## 1.7 iv

phi1 <- 0.85
phi2 <- 0.2

# Characteristic equation (using polynomial coefficients)
coeff <- c(1, -phi1, -phi2) # Note: the signs of phi1 and phi2 are inverted in the characteristic equation

# Solve for the roots of the polynomial
roots <- polyroot(coeff)

# Check if roots are inside the unit circle (stationarity condition)
stationary <- all(abs(roots) > 1)

# Print the results
print(roots)
print(paste0("Stationary: ", stationary))

# Set seed for reproducibility
set.seed(123)  # Feel free to change the seed number

# Number of realizations
n_realizations <- 5

# Number of observations
n <- 200

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations) 
for (i in 1:n_realizations) {
  # Temporary storage for a single realization
  x <- numeric(n)

  # Generate the first two values (adjust if you have specific starting values)
  x[1] <- rnorm(1) 
  x[2] <- phi1 * x[1] + rnorm(1) 

  # Generate the rest of the series
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }

  simulations[, i] <- x
}

# Calculate variances of each realization
realization_variances <- apply(simulations, 2, var) 

# Print the variances
print(realization_variances)

# Plot the realizations
matplot(simulations, type = "l", col = 1:n_realizations, 
        main = "Realizations of AR(2) Process", xlab = "Time", ylab = "Value")


# Number of lags for ACF
nlag <- 30

# Function to calculate theoretical ACF (Yule-Walker implementation)
calculate_theoretical_acf <- function(phi1, phi2, nlag) {
  if (nlag < 0) {
    stop("nlag must be non-negative")
  }

  acf_values <- numeric(nlag + 1)

  acf_values[1] <- 1 
  if (nlag >= 1) {
    acf_values[2] <- phi1 / (1 - phi2) 
  }

  for (k in 3:(nlag + 1)) {
    acf_values[k] <- phi1 * acf_values[k - 1] + phi2 * acf_values[k - 2]
  }

  return(acf_values)
}

# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)


# Plot the theoretical ACF
plot(theoretical_acf, type = "h", lwd = 2, col = "blue",
     main = "Theoretical ACF of AR(2) Process", xlab = "Lag", ylab = "ACF")

# Add reference line at zero
abline(h = 0, col = "gray")


# Calculate theoretical ACF 
theoretical_acf <- calculate_theoretical_acf(phi1, phi2, nlag)

# Function to calculate empirical ACF
empirical_acf <- function(x, nlag) {
  acf(x, lag.max = nlag, type = "correlation", plot = FALSE)$acf 
}

# Calculate empirical ACFs for each realization
empirical_acfs <- apply(simulations, 2, empirical_acf, nlag = nlag)


# Calculate variances of the ACF at each lag
acf_variances <- apply(empirical_acfs, 1, var)  # Note: empirical_acfs might need to be transposed depending on its structure

# Print the variances of the ACF
print(acf_variances)

# Optionally, plot the variances of the ACF
plot(acf_variances, type = "o", pch = 20, col = "blue", xlab = "Lag", ylab = "Variance of ACF",
     main = "Variance of ACF at Each Lag")



# Determine layout for multiple plots
n_plots <- n_realizations  
plot_rows <- ceiling(sqrt(n_plots)) # Adjust if you want a different layout
plot_cols <- ceiling(n_plots / plot_rows) 

# Set up the plot with correct layout
par(mfrow = c(plot_rows, plot_cols)) 

# Confidence level for the ACF confidence intervals (95% confidence)
z <- 1.96
conf_int <- z / sqrt(n)

# Plot each ACF individually with confidence intervals
for (i in 1:n_realizations) {
  plot(empirical_acfs[, i], type = "h", lwd = 2, 
       main = paste0("Realization ", i), xlab = "Lag", ylab = "ACF")

  # Add theoretical ACF (with a different color)
  lines(theoretical_acf, col = "blue", lwd = 2)

  abline(h = 0, col = "gray") # Zero line for reference
  
  # Add confidence intervals
  abline(h = conf_int, col = "red", lty = 2) # Upper confidence interval
  abline(h = -conf_int, col = "red", lty = 2) # Lower confidence interval
}

# Reset to default single plotting area
par(mfrow = c(1, 1))


set.seed(123)

# Simulate the realizations
simulations <- matrix(0, nrow = n, ncol = n_realizations)
for (i in 1:n_realizations) {
  x <- numeric(n)
  x[1:2] <- rnorm(2)
  for (t in 3:n) {
    x[t] <- phi1 * x[t - 1] + phi2 * x[t - 2] + rnorm(1)
  }
  simulations[, i] <- x
}

plot(NULL, xlim=c(1, nlag), ylim=c(-1, 1), xlab="Lag", ylab="ACF",
     main="Empirical ACF of AR(2) Realizations")
abline(h=0, col="gray")
colors <- rainbow(n_realizations)

for (i in 1:n_realizations) {
  # Ensure that you get 30 lags exactly
  acf_values <- acf(simulations[, i], lag.max=nlag, plot = FALSE)$acf[1:(nlag+1)]
  lines(1:nlag, acf_values[-1], type="o", col=colors[i], pch=20, lwd=2)
}

# Add confidence interval lines
conf_int <- qnorm((1 + 0.95) / 2) / sqrt(n)
abline(h = c(-conf_int, conf_int), col = "blue", lty = "dotted")

# Add a legend to the top right
legend("topright", legend=paste("Realization", 1:n_realizations), 
       col=colors, lty=1, lwd=2, cex=0.8)
