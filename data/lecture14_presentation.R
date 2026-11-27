# Lecture 14 presentation script: runs the glmnet code from slides 13-14 of the 
# lecture slides posted for Wednesday, November 18.

# Run with the working directory set to whichever folder contains 
# "fredmd_transformed.csv". 

# Needs glmnet: install.packages("glmnet")

# ---- Slide 13: forming the data ----
library(tidyverse)

# Target: next month's IP growth. Each row is dated by the predictor month.
fm <- read_csv("fredmd_transformed.csv") |> mutate(IP_next = lead(INDPRO))

X <- fm |> select(-date, -IP_next) |> as.matrix()   # predictors: all series this month
y <- fm$IP_next                                     # target: IP growth next month

# Training sample (targets 1960-1979) and validation sample (targets 1980-1989)
train   <- fm$date < as.Date("1979-12-01")
val     <- fm$date >= as.Date("1979-12-01") & fm$date < as.Date("1989-12-01")
X_train <- X[train, ]
y_train <- y[train]
X_val   <- X[val, ]
y_val   <- y[val]

# ---- Slide 14: estimating and forecasting ----
library(glmnet)

# Estimate the lasso for a grid of lambda values (alpha = 0 would be ridge)
grid <- 10^seq(2, -3, length.out = 60)
fit <- glmnet(X_train, y_train, alpha = 1, lambda = grid)

# Forecast the validation sample with one value of lambda (s = 0.1)
fcst <- predict(fit, X_val, s = 0.1)
sqrt(mean((y_val - fcst)^2))   # validation RMSE

# The parameters at that lambda (the lasso sets many to exactly zero)
coef(fit, s = 0.1)
