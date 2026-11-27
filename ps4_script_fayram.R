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