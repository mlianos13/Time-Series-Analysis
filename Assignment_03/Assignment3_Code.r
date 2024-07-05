library(ggplot2)
library(gridExtra)
library(dplyr)
library(TSA)
library(marima2)


## Source the files in the "functions" folder
files <- dir("functions",full.names=TRUE)
for(i in 1:length(files)) source(files[i])

# Load your data
X <- read.table("experiment1.csv", sep=",", header=TRUE)
X$t <- X$t - X$t[1]

##--------------------------------------------------------------------------------.1.--------------------------------------------##

## 1.1 

# Make plots
# Calculate maximum values for scaling
max_Ta <- max(X$Ta)
max_Tinner <- max(X$Tinner)
max_Touter <- max(X$Touter)
max_Pouter <- max(X$Pouter)

# Combined plot: Ta, Tinner, Touter
plot1 <- ggplot(X, aes(x = t)) +
    geom_line(aes(y = Ta, color = "Ta"), size = 0.8) +
    geom_line(aes(y = Tinner, color = "Tinner"), size = 0.8) +
    geom_line(aes(y = Touter, color = "Touter"), size = 0.8) +
    labs(x = "Time (s)", y = "Temperature (°C)") +
    theme_minimal() +
    scale_color_manual(values = c("Ta" = "black", "Tinner" = "blue", "Touter" = "#006600")) +
    ggtitle("Temperature Over Time")

# Combined plot: Pinner and Pouter with scaled Pouter
plot2 <- X %>%
    mutate(Pouter = ifelse(Pouter > 0, max_Pouter, 0)) %>% 
    ggplot(aes(x = t)) +
    geom_line(aes(y = Pinner, color = "Pinner"), size = 0.8) +
    geom_step(aes(y = Pouter, color = "Pouter"), size = 0.8) +
    labs(x = "Time (s)", y = "Temperature (°C)") +
    theme_minimal() +
    scale_color_manual(values = c("Pinner" = "red", "Pouter" = "#006600")) + 
    ggtitle("Pinner and Pouter")

# Tinner and Pinner with scaled Pinner
plot3 <- X %>%
    mutate(Pinner = ifelse(Pinner > 0, max_Tinner, 0)) %>% 
    ggplot(aes(x = t)) +
    geom_line(aes(y = Tinner, color = "Tinner"), size = 0.8) +
    geom_step(aes(y = Pinner, color = "Pinner"), size = 0.8) +
    labs(x = "Time (s)", y = "Temperature (°C)") +
    theme_minimal() +
    scale_color_manual(values = c("Tinner" = "black", "Pinner" = "red")) +
    ggtitle("Tinner and Pinner")

# Touter and Pouter with scaled Pouter
plot4 <- X %>%
    mutate(Pouter = ifelse(Pouter > 0, max_Touter, 0)) %>% 
    ggplot(aes(x = t)) +
    geom_line(aes(y = Touter, color = "Touter"), size = 0.8) +
    geom_step(aes(y = Pouter, color = "Pouter"), size = 0.8) +
    labs(x = "Time (s)", y = "Temperature (°C)") +
    theme_minimal() +
    scale_color_manual(values = c("Touter" = "black", "Pouter" = "#006600")) +
    ggtitle("Touter and Pouter")

# Arrange plots in a 4x1 grid
grid.arrange(plot1, plot2, plot3, plot4, nrow = 4, ncol = 1)

## 1.3

par(mfrow=c(3, 3))

# Calculate and plot CCFs for the specified variable pairs

# Row 1
ccf(X$Tinner, X$Ta, lag.max = 50, main="CCF(Tinner, Ta)")
ccf(X$Tinner, X$Touter, lag.max = 50, main="CCF(Tinner, Touter)")
ccf(X$Tinner, X$Pinner, lag.max = 50, main="CCF(Tinner, Pinner)")

# Row 2
ccf(X$Touter, X$Ta, lag.max = 50, main="CCF(Touter, Ta)")
ccf(X$Touter, X$Tinner, lag.max = 50, main="CCF(Touter, Tinner)")
ccf(X$Touter, X$Pinner, lag.max = 50, main="CCF(Touter, Pinner)")

# Row 3
ccf(X$Tinner, X$Pouter, lag.max = 50, main="CCF(Tinner, Pouter)")
ccf(X$Touter, X$Pouter, lag.max = 50, main="CCF(Touter, Pouter)")

# Leave the last plot empty
plot.new()

##-----------------------------------------------------------------------------------------------.2.----------------------------------##

## 2.1

# Use the lagdf to make lags
lagdf(X[ ,2:6], 1)

# Add the lags to X
maxlag <- 10
for(i in 1:maxlag){
  tmp <- lagdf(X[ ,2:6], i)
  names(tmp) <- paste0(names(tmp),".l",i)
  X <- cbind(X, tmp)
}

# ----------------------------------------------------------------
# Function for making a formula for lm
ARX("Tinner", c("Ta"), 1)
ARX("Tinner", c("Touter","Pouter", "Pinner"), 1:10)
ARX("Tinner", c("Touter", "Ta", "Pinner","Pouter"), 1:2)

# Fit and print
fit <- lm(ARX("Tinner", c("Ta"), 1), data=X)
summary(fit)


## 2.2

# Validation plot
validate(fit)

# CCF plot
ccfplot(fit, X)

AIC(fit)


## 2.3, 2.4

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:2), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:3), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:4), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:5), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:6), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:7), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)


fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:8), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:9), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)

fit <- lm(ARX("Tinner", c("Touter", "Pinner"), 1:10), data=X)
summary(fit)
# Validation plot
validate(fit)
# CCF plot
ccfplot(fit, X)
# AIC
AIC(fit)


##---------------------------------------------.3.-----------------------------##

## 3.1 

# ARMAX model 
# Note MARIMA don't need the lagged (it makes them internally)
X <- X[ ,1:6]
names(X)
fit <- marima("Tinner ~ AR(1) + Ta(1)", data=X)

summary(fit)

validate(fit)


## 3.2

# ARMAX order 1 with all the relevant inputs
fit <- marima("Tinner ~ AR(1) + Touter(1) + Pinner(1) + MA(1)", data=X)
summary(fit)
validate(fit)

# ARX order 1 with all the relevant inputs
fit <- marima("Tinner ~ AR(1) + Touter(1) + Pinner(1)", data=X)
summary(fit)
validate(fit)


# 3.3 ARMAX model selection

# 


fit <- marima("Tinner ~ AR(1:2) + Pinner(1:2) + MA(1:2)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:3) + Pinner(1:3) + MA(1:3)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:4) + Pinner(1:4) + MA(1:4)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:5) + Pinner(1:5) + MA(1:5)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:6) + Pinner(1:6) + MA(1:6)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:7) + Pinner(1:7) + MA(1:7)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:8) + Pinner(1:8) + MA(1:8)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:9) + Pinner(1:9) + MA(1:9)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:10) + Pinner(1:10) + MA(1:10)", data=X, penalty=2)
summary(fit)
validate(fit)




fit <- marima("Tinner ~ AR(1:2) + Pinner(1:2) + Ta(1:2) + MA(1:2)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:3) + Pinner(1:3) + Ta(1:3) + MA(1:3)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:4) + Pinner(1:4) + Ta(1:4) + MA(1:4)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:5) + Pinner(1:5) + Ta(1:5) + MA(1:5)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:6) + Pinner(1:6) + Ta(1:6) + MA(1:6)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:7) + Pinner(1:7) + Ta(1:7) + MA(1:7)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:8) + Pinner(1:8) + Ta(1:8) + MA(1:8)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:9) + Pinner(1:9) + Ta(1:9) + MA(1:9)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Tinner ~ AR(1:10) + Pinner(1:10) + Ta(1:10) + MA(1:10)", data=X, penalty=2)
summary(fit)
validate(fit)




## 3.5 

# Multi-step forecasts
val <- predict(fit, nstep=nrow(X)-1)
plot(X$Tinner)
lines(val$forecasts[1, ], type="l", col=2)


## 3.6
## armax-simulate-step-response
input <- "Pinner"
Xs <- X[ ,2:6]
#Xs[ ,2:6] <- 0
Xs$Pinner <- 0
Xs$Pouter <- 0
#
Xs$Tinner <- 20
Xs$Touter <- 20
Xs$Ta <- 20
Xs[-1:-10,input] <- Xs[-1:-10,input] + 100
#
val <- predict(fit, Xs, nstep=nrow(X)-1)
#plot(M$data[ ,input], type="l")
yhat <- val$forecasts[1, ]
se <- sqrt(val$pred.var[1,1, ])
plot(yhat, type="l")
lines(yhat[-1] - qnorm(0.975)*se, lty=2)
lines(yhat[-1] + qnorm(0.975)*se, lty=2)









input <- "Pinner"  # Ensure 'Pinner' is a valid input variable
Xs <- X[ ,2:6]     # Get appropriate columns
Xs[ ,2:6] <- 0     # Set other inputs to zero (or sensible values)

# Baseline temperature (modify as needed)
Xs$Tinner <- 20 
Xs$Touter <- 20
Xs$Ta <- 20

# Step Change
step_size <- 100
Xs[-1:-10, input] <- Xs[-1:-10, input] + step_size

# Simulation
val <- predict(fit, Xs, nstep=nrow(X)-1)  # Might need adjusting if X's length changed
yhat <- val$forecasts[1, ]
se <- sqrt(val$pred.var[1,1, ])

# Plot
plot(yhat, type="l", ylim = c(min(yhat) - 10, max(yhat) + 10)) # Adjust y-axis for visibility
lines(yhat[-1] - qnorm(0.975)*se, lty=2)
lines(yhat[-1] + qnorm(0.975)*se, lty=2)

# Extreme Step Change (Illustrative)
extreme_step_size <- 500  # Or a value that makes sense in your context
Xs[-1:-10, input] <- Xs[-1:-10, input] + extreme_step_size 

val_extreme <- predict(fit, Xs, nstep=nrow(X)-1)
yhat_extreme <- val_extreme$forecasts[1, ]

lines(yhat_extreme, type="l", col="red") # Plot the extreme response in a different color




##---------------------------------------------------------------------------------------------.4.----------------------------------##

## 4.1

fit <- marima("Touter ~ AR(1) + Tinner(1) + Pouter(1) + MA(1)", data=X, penalty=2)
summary(fit)
validate(fit)

fit <- marima("Touter ~ AR(1:2) + Tinner(1:2) + Pouter(1:2) + MA(1:2)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:3) + Tinner(1:3) + Pouter(1:3) + MA(1:3)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:4) + Tinner(1:4) + Pouter(1:4) + MA(1:4)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:5) + Tinner(1:5) + Pouter(1:5) + MA(1:5)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:6) + Tinner(1:6) + Pouter(1:6) + MA(1:6)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:7) + Tinner(1:7) + Pouter(1:7) + MA(1:7)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:8) + Tinner(1:8) + Pouter(1:8) + MA(1:8)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:9) + Tinner(1:9) + Pouter(1:9) + MA(1:9)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:10) + Tinner(1:10) + Pouter(1:10) + MA(1:10)", data=X, penalty=2)
summary(fit)
validate(fit)












fit <- marima("Touter ~ AR(1) + Pouter(1) + MA(1)", data=X, penalty=2)
summary(fit)
validate(fit)

fit <- marima("Touter ~ AR(1:2) + Pouter(1:2) + MA(1:2)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:3) + Pouter(1:3) + MA(1:3)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:4) + Pouter(1:4) + MA(1:4)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:5) + Pouter(1:5) + MA(1:5)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:6) + Pouter(1:6) + MA(1:6)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:7) + Pouter(1:7) + MA(1:7)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:8) + Pouter(1:8) + MA(1:8)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:9) + Pouter(1:9) + MA(1:9)", data=X, penalty=2)
summary(fit)
validate(fit)


fit <- marima("Touter ~ AR(1:10) + Pouter(1:10) + MA(1:10)", data=X, penalty=2)
summary(fit)
validate(fit)


### Exc.5 Multi-Output Prediction


##5.1)

# 2nd Iteration without Tinner in the 2nd output
fit <- marima("Tinner ~ AR(1:7) + Ta(1:7) + Pinner(1:7) + MA(1:7)",
                "Touter ~ AR(1:5) + Pouter(1:5) + Ta(1:5) + MA(1:5)", data=X,penalty=2)
summary(fit)
validate(fit)



##5.2)

# Decrease and incerase the order equally
fit <- marima("Tinner ~ AR(1:4) + Ta(1:4) + Pinner(1:4) + MA(1:4)",
                "Touter ~ AR(1:2) + Pouter(1:2) + Ta(1:2) + MA(1:2)", data=X,penalty=2)
summary(fit)
validate(fit)

fit <- marima("Tinner ~ AR(1:5) + Ta(1:5) + Pinner(1:5) + MA(1:5)",
                "Touter ~ AR(1:3) + Pouter(1:3) + Ta(1:3) + MA(1:3)", data=X,penalty=2)
summary(fit)
validate(fit)

fit <- marima("Tinner ~ AR(1:6) + Ta(1:6) + Pinner(1:6) + MA(1:6)",
                "Touter ~ AR(1:4) + Pouter(1:4) + Ta(1:4) + MA(1:4)", data=X,penalty=2)
summary(fit)
validate(fit)

fit <- marima("Tinner ~ AR(1:7) + Ta(1:7) + Pinner(1:7) + MA(1:7)",
                "Touter ~ AR(1:5) + Pouter(1:5) + Ta(1:5) + MA(1:5)", data=X,penalty=2)
summary(fit)
validate(fit)

fit <- marima("Tinner ~ AR(1:8) + Ta(1:8) + Pinner(1:8) + MA(1:8)",
                "Touter ~ AR(1:6) + Pouter(1:6) + Ta(1:6) + MA(1:6)", data=X,penalty=2)
summary(fit)
validate(fit)

fit <- marima("Tinner ~ AR(1:9) + Ta(1:9) + Pinner(1:9) + MA(1:9)",
                "Touter ~ AR(1:7) + Pouter(1:7) + Ta(1:7) + MA(1:7)", data=X,penalty=2)
summary(fit)
validate(fit)

fit <- marima("Tinner ~ AR(1:10) + Ta(1:10) + Pinner(1:10) + MA(1:10)",
                "Touter ~ AR(1:8) + Pouter(1:8) + Ta(1:8) + MA(1:8)", data=X,penalty=2)
summary(fit)
validate(fit)


## 5.3

fit <- marima("Tinner ~ AR(1:5) + Ta(1:5) + Pinner(1:5) + MA(1:5)",
                "Touter ~ AR(1:3) + Pouter(1:3) + Ta(1:3) + MA(1:3)", data=X,penalty=2)

# Multi-step forecasts

# ... (Assume your multi-output model is fitted in 'fit')

predictions <- predict(fit, nstep=nrow(X)-1) # Assuming a multi-output format
Tinner_pred <- predictions$forecasts[1,]
Touter_pred <- predictions$forecasts[2,]
X1 <- cbind(X, Tinner_pred, Touter_pred)  # Add predictions to X

# Assuming val$forecasts is a matrix (rows = timesteps, columns = outputs)
plot1 <- ggplot(X1, aes(x = t)) +
    geom_line(aes(y = Tinner, color = "Tinner"), size = 0.8) +
    geom_line(aes(y = Tinner_pred, color = "Tinner_pred"), size = 0.8) +
    labs(x = "Time (s)", y = "Temperature (°C)") +
    theme_minimal() +
    scale_color_manual(values = c("Tinner" = "blue", "Tinner_pred" = "red")) +
    ggtitle("Temperature Over Time")
plot(plot1)

plot2 <- ggplot(X1, aes(x = t)) +
    geom_line(aes(y = Touter, color = "Touter"), size = 0.8) +
    geom_line(aes(y = Touter_pred, color = "Touter_pred"), size = 0.8) +
    labs(x = "Time (s)", y = "Temperature (°C)") +
    theme_minimal() +
    scale_color_manual(values = c("Touter" = "#006600","Touter_pred" = "orange")) +
    ggtitle("Temperature Over Time")
plot(plot2)


## 5.4


## Pinner
input <- "Pinner"  # Ensure 'Pinner' is a valid input variable
Xs <- X[ ,2:6]     # Get appropriate columns
Xs[ ,2:6] <- 0     # Set other inputs to zero (or sensible values)

# Baseline temperature (modify as needed)
Xs$Tinner <- 20 
Xs$Touter <- 20
Xs$Ta <- 20

# Step Change
step_size <- 100
Xs[-1:-10, input] <- Xs[-1:-10, input] + step_size

# Simulation
predictions <- predict(fit, Xs, nstep=nrow(Xs)-1) # Assuming a multi-output format
Tinner_pred <- predictions$forecasts[1,]
Touter_pred <- predictions$forecasts[2,]

# Plot
df <- data.frame(time = 1:length(Tinner_pred), 
                 Ti = Tinner_pred, 
                 To = Touter_pred)

plot_pred <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = Ti, color = "Tinner"), size = 1) +
  geom_line(aes(y = To, color = "Touter"), size = 1) +
  labs(y = "Temperature", title = "Temperature Predictions") +
  theme_minimal()

plot(plot_pred)


# Pouter
input <- "Pouter"  # Ensure 'Pinner' is a valid input variable
Xs <- X[ ,2:6]     # Get appropriate columns
Xs[ ,2:6] <- 0     # Set other inputs to zero (or sensible values)

# Baseline temperature (modify as needed)
Xs$Tinner <- 20 
Xs$Touter <- 20
Xs$Ta <- 20

# Step Change
step_size <- 100
Xs[-1:-10, input] <- Xs[-1:-10, input] + step_size

# Simulation
predictions <- predict(fit, Xs, nstep=nrow(Xs)-1) # Assuming a multi-output format
Tinner_pred <- predictions$forecasts[1,]
Touter_pred <- predictions$forecasts[2,]

# Plot
df <- data.frame(time = 1:length(Tinner_pred), 
                 Ti = Tinner_pred, 
                 To = Touter_pred)

plot_pred <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = Ti, color = "Tinner"), size = 1) +
  geom_line(aes(y = To, color = "Touter"), size = 1) +
  labs(y = "Temperature", title = "Temperature Predictions") +
  theme_minimal()

plot(plot_pred)




## Pinner and Pouter
input1 <- "Pinner"  # Ensure 'Pinner' is a valid input variable
input <- "Pouter"  # Ensure 'Pinner' is a valid input variable
Xs <- X[ ,2:6]     # Get appropriate columns
Xs[ ,2:6] <- 0     # Set other inputs to zero (or sensible values)

# Baseline temperature (modify as needed)
Xs$Tinner <- 20 
Xs$Touter <- 20
Xs$Ta <- 20

# Step Change
step_size <- 100
step_size1 <- 100
Xs[-1:-10, input] <- Xs[-1:-10, input] + step_size
Xs[-1:-10, input1] <- Xs[-1:-10, input1] + step_size1

# Simulation
predictions <- predict(fit, Xs, nstep=nrow(Xs)-1) # Assuming a multi-output format
Tinner_pred <- predictions$forecasts[1,]
Touter_pred <- predictions$forecasts[2,]

# Plot
df <- data.frame(time = 1:length(Tinner_pred), 
                 Ti = Tinner_pred, 
                 To = Touter_pred)

plot_pred <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = Ti, color = "Tinner"), size = 1) +
  geom_line(aes(y = To, color = "Touter"), size = 1) +
  labs(y = "Temperature", title = "Temperature Predictions") +
  theme_minimal()

plot(plot_pred)