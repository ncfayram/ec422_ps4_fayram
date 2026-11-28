library(tidyverse)
library(fpp3)
library(BVAR)

#load data
liquor_raw <- read_csv(
  "data/MRTSSM4453USN.csv",
  na = ".",
  col_types = cols(
    observation_date = col_date(format = "%Y-%m-%d"),
    MRTSSM4453USN = col_double()
  )
)

#build tibble
liquor <- liquor_raw |>
  transmute(
    Month = yearmonth(observation_date),
    Sales = MRTSSM4453USN
  ) |>
  as_tsibble() |>
  filter(Month <= yearmonth("2019 Dec"))

has_gaps(liquor)


#autoplot
autoplot(liquor)

gg_season(liquor, y = Sales) +
  geom_line() +
  labs(
    title = "Liquor Sales by Month",
    x = "Month",
    y = "Sales in Millions of USD"
  ) +
  theme_minimal()
## The seasonal plot shows an overall upward trend over the coarse of the year. It seems to be a consistent pattern with just the level growing each year. There is an initial crest in liquor sales around the summer months of July and August, and another spike at the end of each year in December.

liquor_fit <- liquor |>
  model(arima = ARIMA(log(Sales)),
        ets   = ETS(log(Sales)))

liquor_fit |> select(arima) |> report()
liquor_fit |> select(ets) |> report()
#The chosen ARIMA is a (1,0,2)(2,1,1). This means that the first lag has significant autocorrelation, and the q value of 2 means shocks effect the series for 2 periods after. The series includes 2 seasonal lags of the differenced series, there is 1 seasonal difference applied, and the model includes 1 seasonal lag of white noise shock.
#The chosen ETS is an (A,A,A), meaning that the error, trend, and seasonal components are all additive and don't interact.

liquor_fit |> select(arima) |> gg_tsresiduals()
#The residuals look like they have pretty substantial autocorrelation but they do appear to average to zero.

liquor_fit |>
  select(arima) |>
  augment() |>
  features(.innov, ljung_box, lag = 24, dof = 6)

#Horse Race
liquor_cv <- liquor |>
  stretch_tsibble(.init = 120, .step = 1) |>
  filter(.id < max(.id))

liquor_cv_fit <- liquor_cv |>
  model(snaive = SNAIVE(Sales),
        ets    = ETS(log(Sales)),
        arima  = ARIMA(log(Sales)))

saveRDS(liquor_cv_fit, "results/liquor_cv_fit.rds")
#There will be 215 forecast origins
liquor_cv_fc <- liquor_cv_fit |>
forecast(h = 1)

liquor_cv_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))
#The first forecast Occurs in Jan 2002 and the final forecast occurs in November 2019.
#The ARIMA model forecasts better. It has a lower measure of both RMSE and MAE than the ETS and benchmark models.

#Diebold-Mariano Test: ARIMA vs. ETS
liquor_errs <- liquor_cv_fc |>
  as_tibble() |>
  left_join(as_tibble(liquor) |> select(Month, actual = Sales),
            by = "Month") |>
  mutate(e = actual - .mean) |>
  select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e)

forecast::dm.test(liquor_errs$arima, liquor_errs$ets, h = 1, power = 2)
#DM stat = -2.9062, p-value = 0.004042

#Diebold-Mariano Test: ARIMA vs. SNAIVE
liquor_errs2 <- liquor_cv_fc |>
  as_tibble() |>
  left_join(as_tibble(liquor) |> select(Month, actual = Sales),
            by = "Month") |>
  mutate(e = actual - .mean) |>
  select(Month, .model, e) |>
  pivot_wider(names_from = .model, values_from = e)

forecast::dm.test(liquor_errs2$arima, liquor_errs2$snaive, h = 1, power = 2)
#DM stat = -10.755, p-value < 2.2e-16
#These DM tests show that the ARIMA model performs significantly better than the ETS model, confirming what we saw in the horse race. ARIMA also fits better than the SNAIVE model by a dramatic amount. Knowing this, I would use the ARIMA model on this dataset.

#Forecast Combination
liquor_combo_fc <- liquor_cv_fit |>
  mutate(combo = (arima + ets) / 2) |>
  forecast(h = 1)

liquor_combo_fc |>
  accuracy(liquor, measures = list(RMSE = RMSE, MAE = MAE))


cor(liquor_errs$arima, liquor_errs$ets)
#The errors are highly correlated, meaning combining them won't necessarily provide more information to get a better fit.

#Part 4
qd_raw <- read_csv("data/fredqd_2026-09.csv", show_col_types = FALSE) |>
  slice(-(1:2))

growth <- function(x) 100 * difference(log(x))

qd <- qd_raw |>
  mutate(date = as.Date(sasdate, "%m/%d/%Y")) |>
  transmute(date,
            gdp    = growth(GDPC1),     # real GDP
            cons   = growth(PCECC96),   # real consumption
            inv    = growth(GPDIC1),    # real investment
            emp    = growth(PAYEMS),    # payroll employment
            infl   = growth(GDPCTPI),   # GDP price index
            unrate = difference(UNRATE),     # unemployment rate
            ffr    = difference(FEDFUNDS),   # federal funds rate
            gs10   = difference(GS10)) |>    # 10-year Treasury yield
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

y <- qd |> select(-date) |> as.matrix()

nrow(y)
head(y)

#TS-CV
N <- nrow(y)
origins <- 120:(N - 1)
actual <- y[121:N, ]
#there will be 119 estimates for each variable
#first forecast in Dec 1988, last forecast Sept 2019

#OLS Vars
var_ols_fc <- function(z, p) {
  e <- embed(z, p + 1)            # columns: y_t, then its lags 1 to p
  M <- ncol(z)                    # number of variables
  B <- coef(lm(e[, 1:M] ~ e[, -(1:M)]))
  x <- c(1, t(z[nrow(z):(nrow(z) - p + 1), ]))   # constant, latest p lags
  x %*% B
}

fc_var1 <- NULL
fc_var4 <- NULL
for (T in origins) {
  fc_var1 <- rbind(fc_var1, var_ols_fc(y[1:T, ], 1))
  fc_var4 <- rbind(fc_var4, var_ols_fc(y[1:T, ], 4))
}

rmse_by_var <- function(f) sqrt(colMeans((actual - f)^2))

rmse_by_var(fc_var1)
rmse_by_var(fc_var4)
#The VAR(1) forecasts better than the VAR(4) for every variable except inflation. Because VAR(1) is a special case of VAR(4), it is able to fit the specific parts of the model without having to include the noisier parts that would also be included in VAR(4).

#BVARs
bvar_cv <- function(lambda) {
  pr <- bv_priors(hyper = "lambda",
                  mn = bv_mn(lambda = bv_lambda(mode = lambda, sd = 0.4,
                                                min = 0.999 * lambda,
                                                max = 1.001 * lambda),
                             b = 0))
  fc <- NULL
  for (T in origins) {
    fit_T <- bvar(y[1:T, ], lags = 4, priors = pr,
                  n_draw = 2000, n_burn = 1000, verbose = FALSE)
    x <- c(1, t(y[T:(T - 3), ]))         # constant and 4 lags of each variable
    fc <- rbind(fc, x %*% coef(fit_T))   # forecast of y_{T+1}
  }
  fc
}

set.seed(422)
fc_bvar_03 <- bvar_cv(0.3)
saveRDS(fc_bvar_03, "results/fc_bvar_03.rds")

set.seed(422)
fc_bvar_005 <- bvar_cv(0.05)

var_results <- tibble(
  Variable = colnames(y),
  `VAR(1)` = rmse_by_var(fc_var1),
  `VAR(4)` = rmse_by_var(fc_var4),
  `BVAR(4), lambda = 0.3` = rmse_by_var(fc_bvar_03),
  `BVAR(4), lambda = 0.05` = rmse_by_var(fc_bvar_005)
)

saveRDS(var_results, "results/var_results.rds")

var_results |> knitr::kable(digits = 3)
#VAR(4) performed better for 7/8 variables than the BVAR(4), with lambda = 0.3. The difference is likely due to inaccuracies in the prior men which biased the estimates incorrectly.
#It performed better than VAR(1) for 3/8 variables.
# The RMSE grows when we move lambda from 0.3 to 0.05. This tells us that heavy shrinkage past a certain point is detrimental to our estimates.
# The OLS is still a fair predictor of the data generating process, because this comparison tells us that the prior does not accurately bias the coefficient estimates.

#DM-Test GDP
e_var4 <- actual[, "gdp"] - fc_var4[, 1]
e_bvar <- actual[, "gdp"] - fc_bvar_03[, 1]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)
#The p-value is 0.03 meaning we can reject the null hypothesis that VAR(4) and BVAR(4), lambda = 0.3 are equally accurate predictors of gdp.
# This confirms what I saw in part 4.3b that VAR(4) is a more accurate predictor than BVAR(4).
e_var4 <- actual[, "inv"] - fc_var4[, 3]
e_bvar <- actual[, "inv"] - fc_bvar_03[, 3]

forecast::dm.test(e_var4, e_bvar, h = 1, power = 2)
#The p-value is 0.1065 meaning that we fail to reject the null hypothesis that VAR(4) and BVAR(4), lambda = 0.3 are equally accurate predictors of inv.

md_raw <- read_csv("data/fredmd_2026-09.csv", show_col_types = FALSE)
tcode <- unlist(md_raw[1, -1])
md <- md_raw[-1, ] |> mutate(date = as.Date(sasdate, "%m/%d/%Y"))

transform_series <- function(x, code) {
  d <- function(v) c(NA, diff(v))
  switch(as.character(code),
         "1" = x,                      # level
         "2" = d(x),                   # first difference
         "3" = d(d(x)),                # second difference
         "4" = log(x),                 # log
         "5" = 100 * d(log(x)),        # growth rate
         "6" = 100 * d(d(log(x))),     # change in the growth rate
         "7" = d(x / dplyr::lag(x) - 1))
}

fm <- map2_dfc(md[names(tcode)], tcode, transform_series) |>
  mutate(date = md$date, .before = 1) |>
  filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))

fm <- fm[, colSums(is.na(fm)) == 0]

dim(fm)
# INDPRO and PAYEMS are transformed as growth rates, CPIAUCSL is transformed to changes in growth rates, and FEDFUNDS is first differenced. The difference in transformations for INDPRO and FEDFUNDS is likely due to the fact that the federal funds rate is already measured as a percentage whereas industy production is measured as a quantity.

#Principal components
Z <- fm |> select(-date) |> as.matrix()

pc <- prcomp(scale(Z))

summary(pc)$importance[, 1:5]
#The first 3 components make up 29.45% of the variation in the series.

F1 <- pc$x[, 1] * sign(cor(pc$x[, 1], fm$INDPRO))

tibble(date = fm$date, F1 = F1) |>
  ggplot(aes(date, F1)) +
  geom_line() +
  labs(title = "First Factor from FRED-MD", x = NULL, y = NULL) +
  theme_minimal()

sort(abs(cor(pc$x[, 1], Z))[1, ], decreasing = TRUE)[1:8]
sort(abs(cor(pc$x[, 2], Z))[1, ], decreasing = TRUE)[1:8]
#The first factor seems to capture something closely related to production.
#The second series seems to be related to consumption expenditures.

#Horse Race
ip <- fm$INDPRO
n <- nrow(Z)
init <- which(fm$date == as.Date("1989-12-01"))

e_ar <- NULL
e_fa <- NULL
for (T in init:(n - 1)) {
  F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:3]         # factors using data through T
  d <- data.frame(y_next = ip[2:T], y = ip[1:(T - 1)], F_T[1:(T - 1), ])
  new <- data.frame(y = ip[T], t(F_T[T, ]))
  e_ar <- c(e_ar, ip[T + 1] - predict(lm(y_next ~ y, d), new))
  e_fa <- c(e_fa, ip[T + 1] - predict(lm(y_next ~ ., d), new))
}

c(RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)))
#That would be cheating because it uses information that wouldn't have been available at the time we are forecasting from.

direct_race <- function(series, h = 12, r = 3) {
  lv <- md |>
    mutate(log_x = log(.data[[series]]),
           g = 100 * (log_x - lag(log_x))) |>
    filter(date >= as.Date("1960-03-01"), date <= as.Date("2019-12-01"))
  yh <- 100 * (lead(lv$log_x, h) - lv$log_x)      # growth over the next h months
  g  <- lv$g                                      # growth over the latest month
  
  e_ar <- NULL
  e_fa <- NULL
  for (T in init:(n - h)) {
    F_T <- prcomp(scale(Z[1:T, ]))$x[, 1:r, drop = FALSE]
    i <- 1:(T - h)
    d <- data.frame(yh = yh[i], g = g[i], F_T[i, , drop = FALSE])
    new <- data.frame(g = g[T], F_T[T, , drop = FALSE])
    e_ar <- c(e_ar, yh[T] - predict(lm(yh ~ g, d), new))
    e_fa <- c(e_fa, yh[T] - predict(lm(yh ~ ., d), new))
  }
  tibble(series, forecasts = length(e_ar),
         RMSE_AR = sqrt(mean(e_ar^2)), RMSE_FA = sqrt(mean(e_fa^2)),
         MAE_AR = mean(abs(e_ar)), MAE_FA = mean(abs(e_fa)),
         DM_pvalue = forecast::dm.test(e_ar, e_fa, h = h, power = 2,
                                       varestimator = "bartlett")$p.value)
}

direct_race("CPIAUCSL")
#The factors do improve the models. The RMSE is decreased by 10.43%.The difference is statistically significant given that the DM p-value is below 0.05.