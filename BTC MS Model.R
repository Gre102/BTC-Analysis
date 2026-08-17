# Libraries ----

library(quantmod)
library(dplyr)
library(ggplot2)
library(moments)
library(zoo)
library(car)
library(FinTS)
library(tseries)
library(strucchange)
library(rugarch)
library(MSwM)

# BTC Data ----
## 1. Data Pull ----

btc_raw <- getSymbols(
  "BTC-USD",
  src = "yahoo",
  from = as.Date("2016-01-04"),
  to = as.Date("2026-07-31"),
  auto.assign = FALSE
)

# Convert xts object to tibble
raw_combined <- tibble::tibble(
  Date = as.Date(zoo::index(btc_raw)),
  BTC = as.numeric(Cl(btc_raw))
) |>
  arrange(Date)


## 2. Log Returns ----

full_log <- raw_combined |>
  mutate(
    log_BTC = log(BTC),
    log_return = log(BTC / lag(BTC))
  )


## 3. Data Split ----

break_date <- as.Date("2024-01-11")

# Full pre/post samples

full_pre <- full_log |>
  filter(Date < break_date)

full_post <- full_log |>
  filter(Date >= break_date)


# Balanced pre/post samples
# Use same number of observations on each side of ETF date

n_balanced <- min(
  sum(!is.na(full_pre$log_return)),
  sum(!is.na(full_post$log_return))
)


balanced_pre <- full_pre |>
  filter(!is.na(log_return)) |>
  slice_tail(n = n_balanced)

balanced_post <- full_post |>
  filter(!is.na(log_return)) |>
  slice_head(n = n_balanced)


balanced_full <- bind_rows(
  balanced_pre,
  balanced_post
) |>
  arrange(Date)

# Price Plot ----

ggplot(raw_combined, aes(x = Date, y = BTC)) +
  geom_line() +
  labs(
    title = "Bitcoin Price",
    x = "Date",
    y = "BTC/USD"
  ) +
  theme_minimal()


# Log Return Plot ----

ggplot(full_log, aes(x = Date, y = log_return)) +
  geom_line() +
  labs(
    title = "Bitcoin Daily Log Returns",
    x = "Date",
    y = "Log Return"
  ) +
  theme_minimal()


# Descriptive Statistics ----

summary_data <- bind_rows(
  full_log |>
    select(Date, log_return) |>
    mutate(Period = "Full"),
  
  full_pre |>
    select(Date, log_return) |>
    mutate(Period = "Full_PreETF"),
  
  full_post |>
    select(Date, log_return) |>
    mutate(Period = "Full_PostETF"),
  
  balanced_pre |>
    select(Date, log_return) |>
    mutate(Period = "Bal_PreETF"),
  
  balanced_post |>
    select(Date, log_return) |>
    mutate(Period = "Bal_PostETF")
)


summary_table <- summary_data |>
  group_by(Period) |>
  summarise(
    N = sum(!is.na(log_return)),
    Mean = mean(log_return, na.rm = TRUE),
    Median = median(log_return, na.rm = TRUE),
    SD = sd(log_return, na.rm = TRUE),
    Min = min(log_return, na.rm = TRUE),
    Q05 = quantile(log_return, 0.05, na.rm = TRUE),
    Q25 = quantile(log_return, 0.25, na.rm = TRUE),
    Q75 = quantile(log_return, 0.75, na.rm = TRUE),
    Q95 = quantile(log_return, 0.95, na.rm = TRUE),
    Max = max(log_return, na.rm = TRUE),
    Skewness = skewness(log_return, na.rm = TRUE),
    Kurtosis = kurtosis(log_return, na.rm = TRUE)
  )

summary_table

# Graphs ----
## 1.Histogram ----

ggplot(summary_data, aes(x = log_return)) +
  geom_histogram(bins = 50) +
  facet_wrap(~ Period, scales = "free_y") +
  labs(
    title = "Distribution of Bitcoin Log Returns",
    x = "Daily Log Return",
    y = "Frequency"
  ) +
  theme_minimal()


## 2.Density ----

ggplot(summary_data, aes(x = log_return)) +
  geom_density() +
  facet_wrap(~ Period) +
  labs(
    title = "Density of Bitcoin Log Returns",
    x = "Daily Log Return",
    y = "Density"
  ) +
  theme_minimal()


## 3.Pre/Post Histogram ----

summary_data |>
  filter(Period %in% c("Bal_PreETF", "Bal_PostETF")) |>
  ggplot(aes(
    x = log_return,
    fill = Period
  )) +
  geom_histogram(
    bins = 50,
    alpha = 0.5,
    position = "identity"
  ) +
  labs(
    title = "Bitcoin Return Distribution: Pre- vs Post-ETF",
    x = "Daily Log Return",
    y = "Frequency",
    fill = "Period"
  ) +
  theme_minimal()

## 4.QQ Plot ----

summary_data |>
  filter(Period %in% c("Bal_PreETF", "Bal_PostETF")) |>
  ggplot(aes(sample = log_return)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Period) +
  labs(
    title = "QQ Plots of Bitcoin Returns",
    x = "Theoretical Quantiles",
    y = "Sample Quantiles"
  ) +
  theme_minimal()

## 5.Annualized 30d Realized Vol ----


full_log <- full_log |>
  mutate(
    vol_30d = rollapply(
      log_return,
      width = 30,
      FUN = sd,
      fill = NA,
      align = "right",
      na.rm = TRUE
    ),
    vol_30d_annualized = vol_30d * sqrt(365)
  )

ggplot(full_log, aes(x = Date, y = vol_30d_annualized)) +
  geom_line() +
  geom_vline(
    xintercept = as.Date("2024-01-11"),
    linetype = "dashed"
  ) +
  labs(
    title = "Bitcoin 30-Day Rolling Volatility",
    x = "Date",
    y = "Annualized Volatility"
  ) +
  theme_minimal()

# Distribution Tests ----
## 1. Brown–Forsythe test ----
# check whether different groups have equal variances. Uses median instead of mean
# to avoid noise from outliers 
# H0 is that bothvariances are equal on both periods. 

bf_data <- summary_data |>
  filter(Period %in% c("Bal_PreETF", "Bal_PostETF")) |>
  filter(!is.na(log_return))

leveneTest(
  log_return ~ Period,
  data = bf_data,
  center = median
)


bf_425 <- bind_rows(
  garch_pre |> mutate(Period = "Pre_ETF"),
  garch_post |> mutate(Period = "Post_ETF")
)

car::leveneTest(
  log_return ~ Period,
  data = bf_425,
  center = median
)

## 2. Welch test  ----
# Checking mean equality 

welch_test <- t.test(
  log_return ~ Period,
  data = bf_data,
  var.equal = FALSE
)

welch_test


## 3. KS test  ----

pre_returns <- bf_data |>
  filter(Period == "Bal_PreETF") |>
  pull(log_return)

post_returns <- bf_data |>
  filter(Period == "Bal_PostETF") |>
  pull(log_return)

ks_test <- ks.test(
  pre_returns,
  post_returns
)

ks_test


# Diagnostic tests ----
## 1. ADF Test ----


#Test on log 
adf.test(
  full_log$log_BTC,
  alternative = "stationary"
)

#Test on log return
adf.test(
  na.omit(full_log$log_return),
  alternative = "stationary"
)

## 2. ACF /PACF Test ----
acf(
  na.omit(full_log$log_return),
  lag.max = 30,
  main = "ACF of BTC Daily Log Returns"
)

## 3. ACF Squared Returns ----

acf(
  na.omit(full_log$log_return^2),
  lag.max = 30,
  main = "ACF of Squared BTC Log Returns"
)

## 4. Ljung Box Test  ----

#returns 

Box.test(
  na.omit(full_log$log_return),
  lag = 20,
  type = "Ljung-Box"
)

#Squared returns 

Box.test(
  na.omit(full_log$log_return)^2,
  lag = 20,
  type = "Ljung-Box"
)

## 5. Arch LM test  ----

for (L in c(5, 10, 12, 20)) {
  
  print(
    FinTS::ArchTest(
      na.omit(full_log$log_return),
      lags = L,
      demean = TRUE
    )
  )
}


# Structural break analysis ----

## 1. Bai Perron (unknown)----
# Algorithm searches through possible combinations of dates and finds the 
# segmentation that best fits the data, essentially minimizing the 
# total sum of squared residuals across the regimes.

bp_data <- full_log |>
  filter(!is.na(log_return)) |>
  arrange(Date) |>
  mutate(return_sq = log_return^2)


#Bai Perron on mean 
bp_mean <- breakpoints(
  log_return ~ 1,
  data = bp_data,
  h = 0.15
)


#Bai Perron on Squared Returns 

bp_vol <- breakpoints(
  return_sq ~ 1,
  data = bp_data,
  h = 0.1
)


summary(bp_mean)
plot(bp_mean)

summary(bp_vol)
plot(bp_vol)


bp_data$Date[breakpoints(bp_mean)$breakpoints]
bp_data$Date[breakpoints(bp_vol)$breakpoints]


#RUnning on amended timeframe to avoid major FTX impact in 2022

bp_data_post2022 <- full_log |>
  filter(!is.na(log_return)) |>
  arrange(Date) |>
  mutate(return_sq = log_return^2) |> 
  filter(Date >= as.Date("2022-11-12"))


bp_vol_postftx <- breakpoints(
  return_sq ~ 1,
  data = bp_data_post2022,  
  h = 0.10
)

summary(bp_vol_postftx)
bp_data_post2022$Date[breakpoints(bp_vol_postftx)$breakpoints]


## 2. Chow Test (known date)----

etf_date <- as.Date("2024-01-11")

# Chow test on  return
chow_mean <- sctest(
  log_return ~ 1,
  type = "Chow",
  point = etf_obs,
  data = bp_data_post2022
)
print(chow_mean)

# Chow test on squared returns 
chow_var <- sctest(
  return_sq ~ 1,
  type = "Chow",
  point = etf_obs,
  data = bp_data_post2022
)
print(chow_var)

## 3. CUMSUM on variance----
# CUSUM on variance

cusum_var <- efp(return_sq ~ 1, 
                 type = "OLS-CUSUM", 
                 data = bp_data_post2022)
plot(cusum_var, main = "CUSUM: Squared Returns (post-FTX)")


# GARCH ----
## 1. GARCH(1,1)  ----

garch_data_post2022 <- full_log |>
  filter(
    Date > as.Date("2022-11-11"),
    !is.na(log_return)
  ) |>
  arrange(Date) |> 
  mutate(
    ETF_dummy = if_else(
      Date >= as.Date("2024-01-11"),
      1,
      0
    )
  )

garch_spec <- ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

garch_fit <- ugarchfit(
  spec = garch_spec,
  data = garch_data_post2022$log_return,
  solver = "hybrid"
)

garch_fit #Summary of fit 
coef(garch_fit)
persistence(garch_fit)
infocriteria(garch_fit)

std_resid <- as.numeric(
  residuals(garch_fit, standardize = TRUE)
)
Box.test(
  std_resid,
  lag = 20,
  type = "Ljung-Box"
)
Box.test(
  std_resid^2,
  lag = 20,
  type = "Ljung-Box"
)
FinTS::ArchTest(
  std_resid,
  lags = 12,
  demean = TRUE
)


### 1.1 ETF Dummy ----

etf_xreg <- matrix(
  garch_data_post2022$ETF_dummy,
  ncol = 1
)

garch_etf_spec <- ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1),
    external.regressors = etf_xreg
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

garch_etf_fit <- ugarchfit(
  spec = garch_etf_spec,
  data = garch_data_post2022$log_return,
  solver = "hybrid"
)

garch_etf_fit

### 1.2 Pre/Post Garch ----


garch_pre <- garch_data_post2022 |>
  filter(ETF_dummy == 0)

n_regime <- nrow(garch_pre)

garch_post <- garch_data_post2022 |>
  filter(ETF_dummy == 1) |>
  slice_head(n = n_regime)


garch_pre_fit <- ugarchfit(
  spec = garch_spec,
  data = garch_pre$log_return,
  solver = "hybrid"
)

garch_post_fit <- ugarchfit(
  spec = garch_spec,
  data = garch_post$log_return,
  solver = "hybrid"
)
garch_pre_fit
garch_post_fit

## 2. GJR-GARCH(1,1)  ----

gjr_spec <- ugarchspec(
  variance.model = list(
    model = "gjrGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

gjr_fit <- ugarchfit(
  spec = gjr_spec,
  data = garch_data_post2022$log_return,
  solver = "hybrid"
)

gjr_fit
infocriteria(garch_fit)
infocriteria(gjr_fit)

### 1.1 ETF Dummy ----
gjr_etf_spec <- ugarchspec(
  variance.model = list(
    model = "gjrGARCH",
    garchOrder = c(1, 1),
    external.regressors = etf_xreg
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

gjr_etf_fit <- ugarchfit(
  spec = gjr_etf_spec,
  data = garch_data_post2022$log_return,
  solver = "hybrid"
)

gjr_etf_fit



### 1.2 Pre/Post Garch ----
gjr_spec <- ugarchspec(
  variance.model = list(
    model = "gjrGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

# Pre-ETF
gjr_pre_fit <- ugarchfit(
  spec = gjr_spec,
  data = garch_pre$log_return,
  solver = "hybrid"
)

# Post-ETF
gjr_post_fit <- ugarchfit(
  spec = gjr_spec,
  data = garch_post$log_return,
  solver = "hybrid"
)

gjr_pre_fit
gjr_post_fit


### 1.3 LR Test ----

#I made pre and post have the same number of rows so I am reuniting for the LR test
garch_balanced <- bind_rows(
  garch_pre,
  garch_post
) |>
  arrange(Date)

nrow(garch_balanced)
nrow(garch_pre) + nrow(garch_post)

gjr_pooled_bal_fit <- ugarchfit(
  spec = gjr_spec,
  data = garch_balanced$log_return,
  solver = "hybrid"
)
ll_pooled <- likelihood(gjr_pooled_bal_fit)

ll_separate <-
  likelihood(gjr_pre_fit) +
  likelihood(gjr_post_fit)

LR <- 2 * (ll_separate - ll_pooled)

LR

p_LR <- pchisq(
  LR,
  df = 6,
  lower.tail = FALSE
)

p_LR


# Markov-switching mean/variance model ----
ms_base <- lm(
  log_return ~ 1,
  data = garch_data_post2022
)

set.seed(123)

ms_fit <- msmFit(
  ms_base,
  k = 2,
  p = 0,
  sw = c(TRUE, TRUE),
  control = list(
    parallel = FALSE
  )
)

summary(ms_fit)

high_state <- 2
smooth_prob <- ms_fit@Fit@smoProb[-1, , drop = FALSE]

ms_data <- garch_data_post2022 |>
  mutate(
    Prob_HighVol = smooth_prob[, 2]
  )


ms_summary <- ms_data |>
  mutate(
    Period = if_else(
      Date < as.Date("2024-01-11"),
      "Pre_ETF",
      "Post_ETF"
    )
  ) |>
  group_by(Period) |>
  summarise(
    N = n(),
    Avg_HighVol_Prob = mean(Prob_HighVol),
    HighVol_Share = mean(Prob_HighVol > 0.5)
  )

ms_summary

#balanced dates 

ms_balanced <- ms_data |>
  filter(Date %in% c(garch_pre$Date, garch_post$Date)) |>
  mutate(
    Period = if_else(
      Date < as.Date("2024-01-11"),
      "Pre_ETF",
      "Post_ETF"
    )
  )

ms_balanced_summary <- ms_balanced |>
  group_by(Period) |>
  summarise(
    N = n(),
    Avg_HighVol_Prob = mean(Prob_HighVol),
    HighVol_Share = mean(Prob_HighVol > 0.5)
  )

ms_balanced_summary

