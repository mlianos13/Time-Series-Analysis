library(ggplot2)
library(diagram)
library(emmeans)
library(dplyr)
library(MASS)
library(car)
library(multcomp)
library(readxl)
library(tidyverse)
library(fpp2)
library(tidyr)
library(numDeriv)
library(stats4)

# Load the data
data1 <- read.csv("data/rain1.csv")
data2 <- read.csv("data/rain2.csv")
data3 <- read.csv("data/rain3.csv")
data4 <- read.csv("data/rain4.csv")

# Adding a source identifier for each dataset
data1$Source <- 'Event 1'
data2$Source <- 'Event 2'
data3$Source <- 'Event 3'
data4$Source <- 'Event 4'

# Combine all data frames into one
combined_data <- rbind(data1, data2, data3, data4)

# Creating a new variable for the x-axis and scaling the rain input
combined_data$minutes <- as.numeric(as.character(combined_data$minutes))
combined_data$Scaled_u <- combined_data$u * 10  # Scale the rain input for better visibility




## 1.1 

# Plotting the data with facets
ggplot(combined_data, aes(x = minutes)) +
  geom_line(aes(y = Scaled_u, color = "Rain Input (u)"), size = 1) +
  geom_line(aes(y = y, color = "Basin Water Level (y)"), size = 1) +
  labs(title = "Rain Input and Basin Water Level for Four Rain Events",
       x = "Minutes",
       y = "Scaled Rain Input (1000 m3/min) and Water Level (100 m3)") +
  scale_color_manual(values = c("Rain Input (u)" = "#6a00ff", "Basin Water Level (y)" = "red")) +
  facet_wrap(~Source, scales = "free_x") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max(combined_data$Scaled_u, combined_data$y)), oob = scales::rescale_none)


## 2.1 - 2.2

# 3 Model parameters
a <- 0.045          # Transition rate
sigma1 <- 0.4995      # Standard deviation of system noise for each state
sigma2 <- 2.45        # Standard deviation of observation noise

# Transition matrix A
A <- matrix(c(1-a, 0,   0,   0,
               a, 1-a, 0,   0,
               0,   a, 1-a, 0,
               0,   0,   a, 0.98), byrow=TRUE, nrow=4)

# Input matrix B
B <- matrix(c(1, 0, 0, 0), ncol=1)

# Output matrix C
C <- matrix(c(0, 0, 0, 1), nrow=1)

# State-dependent system noise function
G <- function(X) {
  diag(sapply(X, function(x) sqrt(abs(x))))
}

ut <- data1$u

# Initial state vector
X <- c(0, 0, 0, 0)

# Number of time points
n <- nrow(data1)

# Vectors to save the output and states
y <- numeric(n)
states <- matrix(NA, nrow=n, ncol=4)

# Simulation loop
for (i in 1:n) {
  # System noise
  e1 <- mvrnorm(1, mu = rep(0, 4), Sigma = (sigma1^2) * diag(4))
  
  # State update
  X <- A %*% X + B * ut[i] + G(X) %*% e1
  
  # Save states
  states[i,] <- X
  
  # Observation with noise
  y[i] <- C %*% X + rnorm(1, mean = 0, sd = sigma2)
}

# plot of the results
plot(y, type = 'l', main = "Simulated Basin Water Level", xlab = "Time (minutes)", ylab = "Water Level")





## 2.4

ut <- data1$u  # Assuming 'u' column contains the input rain data
n <- nrow(data1)  # Number of time points

# Initial state vector
X0 <- c(0, 0, 0, 0)

# Define baseline parameters
a_base <- 0.045
sigma1_base <- 0.4995
sigma2_base <- 2.45

# Transition matrix A function
A_matrix <- function(a) {
  matrix(c(1-a, 0,   0,   0,
           a, 1-a, 0,   0,
           0,   a, 1-a, 0,
           0,   0,   a, 0.98), byrow=TRUE, nrow=4)
}

# Simulation function
run_simulation <- function(a, sigma1, sigma2, label, category) {
  A <- A_matrix(a)
  X <- X0
  y <- numeric(n)
  for (i in 1:n) {
    # System noise
    e1 <- mvrnorm(1, mu = rep(0, 4), Sigma = (sigma1^2) * diag(4))
    # State update
    X <- A %*% X + B * ut[i] + G(X) %*% e1
    # Observation with noise
    y[i] <- C %*% X + rnorm(1, mean = 0, sd = sigma2)
  }
  data.frame(Time = 1:n, WaterLevel = y, Label = label, Category = category)
}

# Parameters variations and multiple simulations
results <- data.frame()
for (a in c(0.01, 0.1)) {
  results <- rbind(results, run_simulation(a, sigma1_base, sigma2_base, paste("a =", a), "Transition Rate (a)"))
}
for (sigma1 in c(0.001, 1)) {
  results <- rbind(results, run_simulation(a_base, sigma1, sigma2_base, paste("sigma1 =", sigma1), "System Noise (sigma1)"))
}
for (sigma2 in c(0.1, 5)) {
  results <- rbind(results, run_simulation(a_base, sigma1_base, sigma2, paste("sigma2 =", sigma2), "Observation Noise (sigma2)"))
}

# Plotting the results with facets
ggplot(results, aes(x = Time, y = WaterLevel, color = Label)) +
  geom_line(size = 0.9) +  # Adjust the size as needed for line thickness
  facet_wrap(~Category, scales = "free_y", ncol = 1) +
  labs(title = "Impact of Varying Parameters on Basin Water Level Simulation",
       x = "Time (minutes)",
       y = "Simulated Water Level") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())



#--------------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------------------------------------------#


## 3.


# 3 Model parameters
data_list <- list(
    Event1 = read.csv("data/rain1.csv"),
    Event2 = read.csv("data/rain2.csv"),
    Event3 = read.csv("data/rain3.csv"),
    Event4 = read.csv("data/rain4.csv")
)
names(data_list) <- paste("Event", 1:4)

# Define the simulation function
run_simulation <- function(a, sigma1, sigma2, data) {
    n <- nrow(data)
    if (n == 0) {
        stop("Data is empty or not loaded correctly.")
    }
    X <- c(0, 0, 0, 0)  # Initial state vector
    y <- numeric(n)
    ut <- data$u
    A <- matrix(c(1-a, 0, 0, 0, a, 1-a, 0, 0, 0, a, 1-a, 0, 0, 0, a, 0.98), byrow=TRUE, nrow=4)
    B <- matrix(c(1, 0, 0, 0), ncol=1)
    C <- matrix(c(0, 0, 0, 1), nrow=1)
    
    for (i in 1:n) {
        e1 <- mvrnorm(1, mu = rep(0, 4), Sigma = (sigma1^2) * diag(4))
        X <- A %*% X + B * ut[i] + diag(sapply(X, function(x) sqrt(max(x, 0)))) %*% e1
        y[i] <- C %*% X + rnorm(1, mean = 0, sd = sigma2)
    }
    return(y)
}

# Set specific parameters
a = 0.0325
sigma1 = 0.001
sigma2 = 2.5

# Plotting
par(mfrow = c(2, 2))  # Set up the plotting area to have 2 rows and 2 columns
for (event_name in names(data_list)) {
    event_data <- data_list[[event_name]]
    simulated_y <- run_simulation(a, sigma1, sigma2, event_data)
    plot(simulated_y, type = 'l', main = paste(event_name, "\na =", a, ", sigma1 =", sigma1, ", sigma2 =", sigma2), xlab = "Time (minutes)", ylab = "Water Level")
}

# Reset the plotting parameters
par(mfrow = c(1, 1))


#--------------------------------------------------------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------------------------------------------------------#

## - 4 -###


##---------------------------------------------------------------------------------------##
##-----------------------------------data1-----------------------------------------------##
##---------------------------------------------------------------------------------------##





# Define parameters
n <- nrow(data1)  # Adjust this based on your actual data

# Model matrices initialization
A <- matrix(c(1-0.07, 0, 0, 0,
              0.07, 1-0.07, 0, 0,
              0, 0.07, 1-0.07, 0,
              0, 0, 0.07, 0.98), nrow=4, byrow=TRUE)
B <- matrix(c(1, 0, 0, 0), ncol=1)
C <- matrix(c(0, 0, 0, 1), nrow=1)
Sigma1 <- diag(rep(0.001, 4))  # Initial sigma1, can be modified per current state
Sigma2 <- 1.8^2
X0 <- matrix(c(0, 0, 0, 0), ncol=1)  # Ensure X0 is a column vector for matrix operations

# Data vectors
u <- data1$u  # Adjust 'data1' with your actual data frame name and column
y <- data1$y

# Simulation with corrected state-dependent noise
X <- X0
y_simulated <- numeric(n)
for (i in 2:n) {
    # State-dependent scaling matrix for the noise
    G <- diag(abs(X[,1]))  # Scaling factor based on the absolute values of the state
    noise <- mvrnorm(1, mu = rep(0, 4), Sigma = G * Sigma1)  # Generate scaled noise
    X <- A %*% X + B * u[i-1] + noise
    y_simulated[i] <- as.numeric(C %*% X + sqrt(Sigma2) * rnorm(1))
}

#  implementation considering the adjusted noise handling
negloglik <- function(prm){
    # Parameters dynamically updated within the function
    A[1, 1] <- A[2, 2] <- A[3, 3] <- 1 - prm["a"]
    A[2, 1] <- A[3, 2] <- A[4, 3] <- prm["a"]
    Sigma1_mod <- diag(rep(prm["sigma1"], 4))  # Adjusted noise level

    X <- X0
    SigmaX <- diag(rep(1000, 4))
    lik <- rep(NA, n)
    
    for (i in 1:n) {
        G <- diag(abs(X[,1]))  # State-dependent noise scaling
        Sigma1_scaled <- G * Sigma1_mod * G
        X_pred <- A %*% X + B * u[i]
        SigmaX_pred <- A %*% SigmaX %*% t(A) + Sigma1_scaled
        SigmaY <- C %*% SigmaX_pred %*% t(C) + prm["sigma2"]^2
        lik[i] <- dnorm(y[i], mean = as.numeric(C %*% X_pred), sd = sqrt(SigmaY))
        K <- SigmaX_pred %*% t(C) %*% solve(SigmaY)
        X <- X_pred + K %*% (y[i] - C %*% X_pred)
        SigmaX <- SigmaX_pred - K %*% C %*% SigmaX_pred
    }
    
    return(-sum(log(lik), na.rm = TRUE))
}


# Calculate likelihoods at optimized, lower, and upper bound parameters
optimized_params <- c(a = fit$par["a"], sigma1 = fit$par["sigma1"], sigma2 = fit$par["sigma2"])
lower_bound_params <- c(a = 0.01, sigma1 = 0.0001, sigma2 = 0.1)
upper_bound_params <- c(a = 0.1, sigma1 = 1, sigma2 = 5)

likelihood_optimized <- negloglik(optimized_params)
likelihood_lower_bound <- negloglik(lower_bound_params)
likelihood_upper_bound <- negloglik(upper_bound_params)

# Print the calculated likelihoods
print(paste("Likelihood at optimized parameters:", likelihood_optimized))
print(paste("Likelihood at lower bound parameters:", likelihood_lower_bound))
print(paste("Likelihood at upper bound parameters:", likelihood_upper_bound))


# Parameter optimization call
fit <- optim(c(a = 0.07, sigma1 = 0.001, sigma2 = 1.8), negloglik, method = "L-BFGS-B", lower = c(0.01, 0.0001, 0.1), upper = c(0.1, 1, 5))

##plot actual with simulated data
print("Optimized Parameters:")
print(fit$par)

# Update parameters based on optimization results
optimized_a = fit$par["a"]
optimized_sigma1 = fit$par["sigma1"]
optimized_sigma2 = fit$par["sigma2"]

# Update model matrices with optimized parameters
A[1, 1] <- A[2, 2] <- A[3, 3] <- 1 - optimized_a
A[2, 1] <- A[3, 2] <- A[4, 3] <- optimized_a
Sigma1_opt <- diag(rep(optimized_sigma1, 4))

# Rerun simulation with optimized parameters
X <- X0
y_optimized <- numeric(n)
for (i in 2:n) {
    G <- diag(abs(X[,1]))  # State-dependent scaling factor
    noise <- mvrnorm(1, mu = rep(0, 4), Sigma = G * Sigma1_opt)  # Generate scaled noise
    X <- A %*% X + B * u[i-1] + noise
    y_optimized[i] <- as.numeric(C %*% X + sqrt(optimized_sigma2) * rnorm(1))
}

# Plot the simulated water levels with optimized parameters
plot(y_optimized, type = "l", col = "red", xlab = "Time (minutes)", ylab = "Simulated Water Level (100 m3) with Optimized Parameters")
lines(y, col = "blue")  # Optionally plot the actual data for comparison
legend("topright", legend=c("Simulated", "Actual"), col=c("red", "blue"), lty=1, cex=0.8)



#Initial values: 'a' = 0.07, 'sigma1' = 0.001, 'sigma2' = 1.8


##---------------------------------------------------------------------------------------##
##-----------------------------------data2-----------------------------------------------##
##---------------------------------------------------------------------------------------##



# Define parameters
n <- nrow(data2)  # Adjust this based on your actual data

# Model matrices initialization
A <- matrix(c(1-0.1, 0, 0, 0,
              0.1, 1-0.1, 0, 0,
              0, 0.1, 1-0.1, 0,
              0, 0, 0.1, 0.98), nrow=4, byrow=TRUE)
B <- matrix(c(1, 0, 0, 0), ncol=1)
C <- matrix(c(0, 0, 0, 1), nrow=1)
Sigma1 <- diag(rep(0.001, 4))  # Initial sigma1, can be modified per current state
Sigma2 <- 1.8^2
X0 <- matrix(c(0, 0, 0, 0), ncol=1)  # Ensure X0 is a column vector for matrix operations

# Data vectors
u <- data2$u  # Adjust 'data1' with your actual data frame name and column
y <- data2$y

# Simulation with corrected state-dependent noise
X <- X0
y_simulated <- numeric(n)
for (i in 2:n) {
    # State-dependent scaling matrix for the noise
    G <- diag(abs(X[,1]))  # Scaling factor based on the absolute values of the state
    noise <- mvrnorm(1, mu = rep(0, 4), Sigma = G * Sigma1)  # Generate scaled noise
    X <- A %*% X + B * u[i-1] + noise
    y_simulated[i] <- as.numeric(C %*% X + sqrt(Sigma2) * rnorm(1))
}


# Kalman filter implementation considering the adjusted noise handling
negloglik <- function(prm){
    # Parameters dynamically updated within the function
    A[1, 1] <- A[2, 2] <- A[3, 3] <- 1 - prm["a"]
    A[2, 1] <- A[3, 2] <- A[4, 3] <- prm["a"]
    Sigma1_mod <- diag(rep(prm["sigma1"], 4))  # Adjusted noise level

    X <- X0
    SigmaX <- diag(rep(1000, 4))
    lik <- rep(NA, n)
    
    for (i in 1:n) {
        G <- diag(abs(X[,1]))  # Apply state-dependent noise scaling
        Sigma1_scaled <- G * Sigma1_mod * G  # Scale Sigma1 according to the state
        X_pred <- A %*% X + B * u[i]
        SigmaX_pred <- A %*% SigmaX %*% t(A) + Sigma1_scaled
        SigmaY <- C %*% SigmaX_pred %*% t(C) + prm["sigma2"]^2
        lik[i] <- dnorm(y[i], mean = as.numeric(C %*% X_pred), sd = sqrt(SigmaY))
        K <- SigmaX_pred %*% t(C) %*% solve(SigmaY)
        X <- X_pred + K %*% (y[i] - C %*% X_pred)
        SigmaX <- SigmaX_pred - K %*% C %*% SigmaX_pred
    }
    
    return(-sum(log(lik), na.rm = TRUE))
}

# Calculate likelihoods at optimized, lower, and upper bound parameters
optimized_params <- c(a = fit$par["a"], sigma1 = fit$par["sigma1"], sigma2 = fit$par["sigma2"])
lower_bound_params <- c(a = 0.01, sigma1 = 0.0001, sigma2 = 0.1)
upper_bound_params <- c(a = 0.1, sigma1 = 1, sigma2 = 5)

likelihood_optimized <- negloglik(optimized_params)
likelihood_lower_bound <- negloglik(lower_bound_params)
likelihood_upper_bound <- negloglik(upper_bound_params)

# Print the calculated likelihoods
print(paste("Likelihood at optimized parameters:", likelihood_optimized))
print(paste("Likelihood at lower bound parameters:", likelihood_lower_bound))
print(paste("Likelihood at upper bound parameters:", likelihood_upper_bound))


##plot actual with simulated data
# Parameter optimization call
fit <- optim(c(a = 0.07, sigma1 = 0.001, sigma2 = 1.8), negloglik, method = "L-BFGS-B", lower = c(0.01, 0.0001, 0.1), upper = c(0.1, 1, 5))
print("Optimized Parameters:")
print(fit$par)

# Update parameters based on optimization results
optimized_a = fit$par["a"]
optimized_sigma1 = fit$par["sigma1"]
optimized_sigma2 = fit$par["sigma2"]

# Update model matrices with optimized parameters
A[1, 1] <- A[2, 2] <- A[3, 3] <- 1 - optimized_a
A[2, 1] <- A[3, 2] <- A[4, 3] <- optimized_a
Sigma1_opt <- diag(rep(optimized_sigma1, 4))

# Rerun simulation with optimized parameters
X <- X0
y_optimized <- numeric(n)
for (i in 2:n) {
    G <- diag(abs(X[,1]))  # State-dependent scaling factor
    noise <- mvrnorm(1, mu = rep(0, 4), Sigma = G * Sigma1_opt)  # Generate scaled noise
    X <- A %*% X + B * u[i-1] + noise
    y_optimized[i] <- as.numeric(C %*% X + sqrt(optimized_sigma2) * rnorm(1))
}

# Plot the simulated water levels with optimized parameters
plot(y_optimized, type = "l", col = "red", xlab = "Time (minutes)", ylab = "Simulated Water Level (100 m3) with Optimized Parameters")
lines(y, col = "blue")  # Optionally plot the actual data for comparison
legend("topright", legend=c("Simulated", "Actual"), col=c("red", "blue"), lty=1, cex=0.8)


# 'a' = 0.1, 'sigma1' = 0.0001, 'sigma2' = 1.8

##---------------------------------------------------------------------------------------##
##-----------------------------------data3-----------------------------------------------##
##---------------------------------------------------------------------------------------##

# Define parameters
n <- nrow(data3)  # Adjust this based on your actual data

# Model matrices initialization
A <- matrix(c(1-0.1, 0, 0, 0,
              0.1, 1-0.1, 0, 0,
              0, 0.1, 1-0.1, 0,
              0, 0, 0.1, 0.98), nrow=4, byrow=TRUE)
B <- matrix(c(1, 0, 0, 0), ncol=1)
C <- matrix(c(0, 0, 0, 1), nrow=1)
Sigma1 <- diag(rep(0.001, 4))  # Initial sigma1, can be modified per current state
Sigma2 <- 2.8^2
X0 <- matrix(c(0, 0, 0, 0), ncol=1)  # Ensure X0 is a column vector for matrix operations

# Data vectors
u <- data3$u  # Adjust 'data1' with your actual data frame name and column
y <- data3$y

# Simulation with corrected state-dependent noise
X <- X0
y_simulated <- numeric(n)
for (i in 2:n) {
    # State-dependent scaling matrix for the noise
    G <- diag(abs(X[,1]))  # Scaling factor based on the absolute values of the state
    noise <- mvrnorm(1, mu = rep(0, 4), Sigma = G * Sigma1)  # Generate scaled noise
    X <- A %*% X + B * u[i-1] + noise
    y_simulated[i] <- as.numeric(C %*% X + sqrt(Sigma2) * rnorm(1))
}
# Kalman filter implementation considering the adjusted noise handling
negloglik <- function(prm){
    # Parameters dynamically updated within the function
    A[1, 1] <- A[2, 2] <- A[3, 3] <- 1 - prm["a"]
    A[2, 1] <- A[3, 2] <- A[4, 3] <- prm["a"]
    Sigma1_mod <- diag(rep(prm["sigma1"], 4))  # Adjusted noise level

    X <- X0
    SigmaX <- diag(rep(1000, 4))
    lik <- rep(NA, n)
    
    for (i in 1:n) {
        G <- diag(abs(X[,1]))  # Apply state-dependent noise scaling
        Sigma1_scaled <- G * Sigma1_mod * G  # Scale Sigma1 according to the state
        X_pred <- A %*% X + B * u[i]
        SigmaX_pred <- A %*% SigmaX %*% t(A) + Sigma1_scaled
        SigmaY <- C %*% SigmaX_pred %*% t(C) + prm["sigma2"]^2
        lik[i] <- dnorm(y[i], mean = as.numeric(C %*% X_pred), sd = sqrt(SigmaY))
        K <- SigmaX_pred %*% t(C) %*% solve(SigmaY)
        X <- X_pred + K %*% (y[i] - C %*% X_pred)
        SigmaX <- SigmaX_pred - K %*% C %*% SigmaX_pred
    }
    
    return(-sum(log(lik), na.rm = TRUE))
}

# Calculate likelihoods at optimized, lower, and upper bound parameters
optimized_params <- c(a = fit$par["a"], sigma1 = fit$par["sigma1"], sigma2 = fit$par["sigma2"])
lower_bound_params <- c(a = 0.01, sigma1 = 0.0001, sigma2 = 0.1)
upper_bound_params <- c(a = 0.1, sigma1 = 1, sigma2 = 5)

likelihood_optimized <- negloglik(optimized_params)
likelihood_lower_bound <- negloglik(lower_bound_params)
likelihood_upper_bound <- negloglik(upper_bound_params)

# Print the calculated likelihoods
print(paste("Likelihood at optimized parameters:", likelihood_optimized))
print(paste("Likelihood at lower bound parameters:", likelihood_lower_bound))
print(paste("Likelihood at upper bound parameters:", likelihood_upper_bound))

##plot actual with simulated data
# Parameter optimization call
fit <- optim(c(a = 0.1, sigma1 = 0.001, sigma2 = 2.8), negloglik, method = "L-BFGS-B", lower = c(0.01, 0.0001, 0.1), upper = c(0.1, 1, 5))
print("Optimized Parameters:")
print(fit$par)

# Update parameters based on optimization results
optimized_a = fit$par["a"]
optimized_sigma1 = fit$par["sigma1"]
optimized_sigma2 = fit$par["sigma2"]

# Update model matrices with optimized parameters
A[1, 1] <- A[2, 2] <- A[3, 3] <- 1 - optimized_a
A[2, 1] <- A[3, 2] <- A[4, 3] <- optimized_a
Sigma1_opt <- diag(rep(optimized_sigma1, 4))

# Rerun simulation with optimized parameters
X <- X0
y_optimized <- numeric(n)
for (i in 2:n) {
    G <- diag(abs(X[,1]))  # State-dependent scaling factor
    noise <- mvrnorm(1, mu = rep(0, 4), Sigma = G * Sigma1_opt)  # Generate scaled noise
    X <- A %*% X + B * u[i-1] + noise
    y_optimized[i] <- as.numeric(C %*% X + sqrt(optimized_sigma2) * rnorm(1))
}

# Plot the simulated water levels with optimized parameters
plot(y_optimized, type = "l", col = "red", xlab = "Time (minutes)", ylab = "Simulated Water Level (100 m3) with Optimized Parameters")
lines(y, col = "blue")  # Optionally plot the actual data for comparison
legend("topright", legend=c("Simulated", "Actual"), col=c("red", "blue"), lty=1, cex=0.8)


#Initial values: 'a' = 0.1, 'sigma1' = 0.001, 'sigma2' = 2.8

##---------------------------------------------------------------------------------------##
##-----------------------------------data4-----------------------------------------------##
##---------------------------------------------------------------------------------------##


# Define parameters
n <- nrow(data4)  # Adjust this based on your actual data

# Model matrices initialization
A <- matrix(c(1-0.0325, 0, 0, 0,
              0.0325, 1-0.0325, 0, 0,
              0, 0.0325, 1-0.0325, 0,
              0, 0, 0.0325, 0.98), nrow=4, byrow=TRUE)
B <- matrix(c(1, 0, 0, 0), ncol=1)
C <- matrix(c(0, 0, 0, 1), nrow=1)
Sigma1 <- diag(rep(0.5, 4))  # Initial sigma1, can be modified per current state
Sigma2 <- 2^2
X0 <- matrix(c(0, 0, 0, 0), ncol=1)  # Ensure X0 is a column vector for matrix operations

# Data vectors
u <- data4$u  # Adjust 'data1' with your actual data frame name and column
y <- data4$y

# Simulation with corrected state-dependent noise
X <- X0
y_simulated <- numeric(n)
for (i in 2:n) {
    # State-dependent scaling matrix for the noise
    G <- diag(abs(X[,1]))  # Scaling factor based on the absolute values of the state
    noise <- mvrnorm(1, mu = rep(0, 4), Sigma = G * Sigma1)  # Generate scaled noise
    X <- A %*% X + B * u[i-1] + noise
    y_simulated[i] <- as.numeric(C %*% X + sqrt(Sigma2) * rnorm(1))
}

# Kalman filter implementation considering the adjusted noise handling
negloglik <- function(prm){
    # Parameters dynamically updated within the function
    A[1, 1] <- A[2, 2] <- A[3, 3] <- 1 - prm["a"]
    A[2, 1] <- A[3, 2] <- A[4, 3] <- prm["a"]
    Sigma1_mod <- diag(rep(prm["sigma1"], 4))  # Adjusted noise level

    X <- X0
    SigmaX <- diag(rep(1000, 4))
    lik <- rep(NA, n)
    
    for (i in 1:n) {
        G <- diag(abs(X[,1]))  # Apply state-dependent noise scaling
        Sigma1_scaled <- G * Sigma1_mod * G  # Scale Sigma1 according to the state
        X_pred <- A %*% X + B * u[i]
        SigmaX_pred <- A %*% SigmaX %*% t(A) + Sigma1_scaled
        SigmaY <- C %*% SigmaX_pred %*% t(C) + prm["sigma2"]^2
        lik[i] <- dnorm(y[i], mean = as.numeric(C %*% X_pred), sd = sqrt(SigmaY))
        K <- SigmaX_pred %*% t(C) %*% solve(SigmaY)
        X <- X_pred + K %*% (y[i] - C %*% X_pred)
        SigmaX <- SigmaX_pred - K %*% C %*% SigmaX_pred
    }
    
    return(-sum(log(lik), na.rm = TRUE))
}

# Calculate likelihoods at optimized, lower, and upper bound parameters
optimized_params <- c(a = fit$par["a"], sigma1 = fit$par["sigma1"], sigma2 = fit$par["sigma2"])
lower_bound_params <- c(a = 0.01, sigma1 = 0.0001, sigma2 = 0.1)
upper_bound_params <- c(a = 0.1, sigma1 = 1, sigma2 = 5)

likelihood_optimized <- negloglik(optimized_params)
likelihood_lower_bound <- negloglik(lower_bound_params)
likelihood_upper_bound <- negloglik(upper_bound_params)

# Print the calculated likelihoods
print(paste("Likelihood at optimized parameters:", likelihood_optimized))
print(paste("Likelihood at lower bound parameters:", likelihood_lower_bound))
print(paste("Likelihood at upper bound parameters:", likelihood_upper_bound))

##plot actual with simulated data
# Parameter optimization call
fit <- optim(c(a = 0.0325, sigma1 = 0.5, sigma2 = 2), negloglik, method = "L-BFGS-B", lower = c(0.01, 0.0001, 0.1), upper = c(0.1, 1, 5))
print("Optimized Parameters:")
print(fit$par)

# Update parameters based on optimization results
optimized_a = fit$par["a"]
optimized_sigma1 = fit$par["sigma1"]
optimized_sigma2 = fit$par["sigma2"]

# Update model matrices with optimized parameters
A[1, 1] <- A[2, 2] <- A[3, 3] <- 1 - optimized_a
A[2, 1] <- A[3, 2] <- A[4, 3] <- optimized_a
Sigma1_opt <- diag(rep(optimized_sigma1, 4))

# Rerun simulation with optimized parameters
X <- X0
y_optimized <- numeric(n)
for (i in 2:n) {
    G <- diag(abs(X[,1]))  # State-dependent scaling factor
    noise <- mvrnorm(1, mu = rep(0, 4), Sigma = G * Sigma1_opt)  # Generate scaled noise
    X <- A %*% X + B * u[i-1] + noise
    y_optimized[i] <- as.numeric(C %*% X + sqrt(optimized_sigma2) * rnorm(1))
}

# Plot the simulated water levels with optimized parameters
plot(y_optimized, type = "l", col = "red", xlab = "Time (minutes)", ylab = "Simulated Water Level (100 m3) with Optimized Parameters")
lines(y, col = "blue")  # Optionally plot the actual data for comparison
legend("topright", legend=c("Simulated", "Actual"), col=c("red", "blue"), lty=1, cex=0.8)


#Initial values: 'a' = 0.0325, 'sigma1' = 0.5, 'sigma2' = 2

##---------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------##
##---------------------------------------------------------------------------------------##

