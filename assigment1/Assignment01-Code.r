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

data <- read_excel("assigment1/DST_BIL54_train.xlsx")
colnames(data)[1] <- "Category"

dt1 <- filter(data, Category == "Drivmidler i alt")
dt1 <- dt1[, -1]
# Create the time variable
months <- seq(as.Date("2018-01-01"), as.Date("2022-11-01"), by="month")
# months <- seq(as.Date("2018-01-01"), as.Date("2022-11-01"), by="month") creates a sequence of dates from January 2018 to November 2022, one for each month.

x <- 2018 + (as.numeric(format(months, "%m")) - 1) / 12 

dt <- pivot_longer(dt1, 
                          cols = everything(),
                          names_to = "Time", 
                          values_to = "Value")

dt$Time <- as.Date(months)

# Now, let's plot
ggplot(dt, aes(x = Time, y = Value)) + 
  geom_line() + 
  labs(title = "Training Data over Time", x = "Time", y = "Value") +
  theme_minimal()



##--------------------------------------------------------------------------------------------------------------------------##

## 2. OLS Model


nparams <- 2 # number of parameters
n <- length(dt$Time) # number of observations
X <- cbind(1, dt$Time) # X is the "design matrix"
print(X)
y <- cbind(dt$Value) # y is vector with observations:
print(y)

## 2.1 Estimate the parameters β0 and β1 using the training data (OLS model). Describe how you estimate the parameters.

# Estimate the parameters β0 and β1 using the training data (OLS model). Describe how you estimate the parameters.
# The OLS estimate of the parameters is given by βˆ = (X'X)^(-1)X'y
# where X is the "design matrix" and y is the vector with observations.

# to estimate parameters we solve the "normal equations":
(OLS <- solve(t(X)%*%X)%*%t(X)%*%y)
# these are the parameter estimates!
theta_0 <- OLS[1]
theta_1 <- OLS[2]


## 2.2 Present the values of the parameter estimates βˆ0 and βˆ1 and their estimated standard errors ˆσβˆ0 and ˆσβˆ1

# now we compute y_hat values for all observations:
yhat_ols <- X%*%OLS
print(yhat_ols)

# plot:
ggplot(dt, aes(x=Time, y=Value)) +
  geom_point() + 
  geom_line(aes(y=yhat_ols), col="red", size=.5) 

# we will now calculate the standard errors on the parameters beta_0 and beta_1:

# first compute residuals:
e_ols <- y - yhat_ols
print(e_ols)

# calculate sum of squared residuals:
RSS_ols <- t(e_ols)%*%e_ols
print(RSS_ols)

# calculate sigma^2:
sigma2_ols <- as.numeric(RSS_ols/(n - nparams)) # = "Residual variance"
sqrt(sigma2_ols) # = "Residual standard error"

# calculate variance-covariance matrix of _parameters_:
V_ols <- sigma2_ols * solve(t(X) %*% X)
print(V_ols)

# the variances of the parameters are the values in the diagonal:
diag(V_ols)

# and the standard errors are given by:
(sqrt(diag(V_ols))) 

se_theta_0 <- (sqrt(diag(V_ols)))[1]
se_theta_1 <- (sqrt(diag(V_ols)))[2]

# now we have both point estimates and standard errors:
# intercept:
theta_0
se_theta_0
# slope:
theta_1
se_theta_1

## 2.3 Make a forecast for the next 12 months - i.e., compute predicted values with corresponding prediction intervals. Present these values in a table.

# now we use the model for predictions on future timepoints
# we use the timepoints from the testdata:
Xtest <- cbind(1, dt$Time)
print(Xtest)

# compute predictions 
y_pred <- Xtest%*%OLS
print(y_pred)

# compute prediction variance-covariance matrix:
Vmatrix_pred <- sigma2_ols*(1+(Xtest%*%solve(t(X)%*%X))%*%t(Xtest))

# the variances of individual predictions are in the diagonal of the matrix above

# compute "prediction intervals"
y_pred_lwr <- y_pred - 1.96*sqrt(diag(Vmatrix_pred))
y_pred_upr <- y_pred + 1.96*sqrt(diag(Vmatrix_pred))

## 2.4 Plot the fitted model together with the training data and the forecasted values (also plot the prediction intervals of the forecasted values).

# plot forecast:
ggplot(dt, aes(x=Time, y=Value)) +
  geom_point() + 
  geom_line(aes(y=yhat_ols), col="red", size=.5) +
  geom_point(data=dt, aes(x=Time,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=dt, aes(x=Time,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.2, fill="red")

# plot WITH test data:
ggplot(dt, aes(x=Time, y=Value)) +
  geom_point() + 
  geom_line(aes(y=yhat_ols), col="red", size=.5) +
  geom_point(data=dt, aes(x=Time,y=y_pred), col="red", size=.5) +
  geom_ribbon(data=dt, aes(x=Time,ymin=y_pred_lwr, ymax=y_pred_upr), inherit.aes=FALSE, alpha=0.2, fill="red") +
  geom_point(data=dt, aes(x=Time,y=Value), col="blue", size=.5)



# Calculate residuals
e_ols <- y - yhat_ols

# Plot residuals against fitted values
ggplot() +
  geom_point(data = data.frame(Fitted = yhat_ols[,1], Residuals = e_ols[,1]), aes(x = Fitted, y = Residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted Values", x = "Fitted Values", y = "Residuals")

# Plot residuals against Time (or any other relevant predictor variable)
ggplot() +
  geom_point(data = data.frame(Time = dt$Time, Residuals = e_ols[,1]), aes(x = Time, y = Residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Time", x = "Time", y = "Residuals")

# Q-Q plot of residuals to assess normality
qqnorm(e_ols[,1])
qqline(e_ols[,1])

# Shapiro-Wilk test for normality of residuals
shapiro.test(e_ols[,1])
