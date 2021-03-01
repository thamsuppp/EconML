library(tidyverse)
library(xts)
library(tseries)
library(vars)
library(dynlm)
library(forecast)
library(TTR)

library(tidyquant)

install.packages('tidyquant')

setwd('C:/Users/isaac/Documents/Isaac/2020-21 Penn Junior/Winter Break 2020/Economic Research PGC/VAR Coding')

# Read Data from CSV
df <- read.csv('fred_vars.csv')

# Take % change of GDP and CPI
df <- df %>% mutate(logGDP = log(GDPC1),
                    logCPI = log(CPIAUCSL)) 
df$GDP_Growth <- diff.xts(df$logGDP, lag = 1)
df$Inflation <- diff.xts(df$logCPI, lag = 1)

# Stationarity testing
for (i in 2:ncol(df)) {
  print(colnames(df)[i])
  print(adf.test(df[2:nrow(df),i]))
}
# Result: GDP_Growth and Unemployment stationary, GS1 and Inflation slightly non-stationary

# Remove all the other variables we don't need
# ** Remove last two quarters because of COVID
df <- df[2:(nrow(df)-2),c(1,8,9,4,3)]
df_copy <- data.frame(df)
row.names(df) <- NULL

# Convert to time series
time_series <- ts(df[,2:ncol(df)], start = c(1952, 7), frequency = 4)
plot(time_series, main = 'Quarterly US Macro Variables')

time_series

# Initial VAR -------------------------------------------------------------

# Getting the VAR matrix
varmat <- as.matrix(df[, 2:ncol(df)])

# Model selection of VAR lag
var_selection <- VARselect(varmat, lag.max = 12, type = 'const')
var_selection$selection

# Fitting the VAR
var <- VAR(varmat, p = 2, type = 'const', season = NULL, exog = NULL)
summary(var)

var$datamat

var_preds <- predict(var, varmat)
var_residuals <- residuals(var)
cor(var_residuals)


# In-sample VAR via OLS ---------------------------------------------------

### Create the lags explicitly and implementing VAR via OLS (for Cross-validation)
df2 <- df %>% mutate(GDP_Growth_lag1 = lag(GDP_Growth),
              GDP_Growth_lag2 = lag(GDP_Growth, 2),
              Inflation_lag1 = lag(Inflation),
              Inflation_lag2 = lag(Inflation, 2),
              UNRATE_lag1 = lag(UNRATE),
              UNRATE_lag2 = lag(UNRATE, 2),
              GS1_lag1 = lag(GS1),
              GS1_lag2 = lag(GS1, 2)) %>% .[3:nrow(.), ]
row.names(df2) <- NULL

# # Testing to see if OLS is equivalent to VAR function
# preds <- data.frame(matrix(ncol = 4, nrow = nrow(df2)))
# names(preds) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')
# all_residuals <- data.frame(matrix(ncol = 4, nrow = nrow(df2)))
# names(all_residuals) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')
# 
regressors <- c('GDP_Growth_lag1', 'Inflation_lag1','UNRATE_lag1', 'GS1_lag1',
                 'GDP_Growth_lag2',  'Inflation_lag2', 'UNRATE_lag2', 'GS1_lag2')
 regressors <- paste(regressors, collapse = ' + ')
# 
# gdp_reg <- lm(paste('GDP_Growth ~', regressors), df2)
# inflation_reg <- lm(paste('Inflation ~', regressors), df2)
# unemp_reg <- lm(paste('UNRATE ~', regressors), df2)
# gs1_reg <- lm(paste('GS1 ~', regressors), df2)
# 
# preds[, 1] <- predict(gdp_reg, df2[, c(6,8,10,12,7,9,11,13)])
# preds[, 2] <- predict(inflation_reg, df2[, c(6,8,10,12,7,9,11,13)])
# preds[, 3] <- predict(unemp_reg, df2[, c(6,8,10,12,7,9,11,13)])
# preds[, 4] <- predict(gs1_reg, df2[, c(6,8,10,12,7,9,11,13)])
# 
# all_residuals <- varmat[3:nrow(varmat), ] - preds
# cor(all_residuals)
# 
# # (Equivalent method of doing the same thing) Testing to see if both residuals are identical
# all_residuals[, 1] <- residuals(gdp_reg)
# all_residuals[, 2] <- residuals(inflation_reg)
# all_residuals[, 3] <- residuals(unemp_reg)
# all_residuals[, 4] <- residuals(gs1_reg)
# cor(all_residuals)


# Leave-one-out VAR -------------------------------------------------------

loo_preds <- data.frame(matrix(ncol = 4, nrow = nrow(df2)))
names(loo_preds) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')

for (i in 1:nrow(df2)) {
  
  # Leave one out of the matrix and estimate the VAR regression
  df_loo <- df2[-i, 2:ncol(df2)]
  left_out <- df2[i, 2:ncol(df2)]
  
  gdp_reg <- lm(paste('GDP_Growth ~', regressors), df_loo)
  inflation_reg <- lm(paste('Inflation ~', regressors), df_loo)
  unemp_reg <- lm(paste('UNRATE ~', regressors), df_loo)
  gs1_reg <- lm(paste('GS1 ~', regressors), df_loo)
  
  # Get out of sample predictions for the left out observation

  loo_preds[i, 1] <- predict(gdp_reg, left_out[5:ncol(left_out)])
  loo_preds[i, 2] <- predict(inflation_reg, left_out[5:ncol(left_out)])
  loo_preds[i, 3] <- predict(unemp_reg, left_out[5:ncol(left_out)])
  loo_preds[i, 4] <- predict(gs1_reg, left_out[5:ncol(left_out)])
  
}

# Get the LOO Residuals
loo_residuals <- varmat[3:nrow(varmat), ]  - loo_preds
cor(loo_residuals)


# K-fold Cross Validation -------------------------------------------------

# Hyperparameter of number of folds
k = round(nrow(df2) / 8) # to get quarter blocks

kfold_preds <- data.frame(matrix(ncol = 4, nrow = nrow(df2)))
names(kfold_preds) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')

# Get the fold cuts 
fold_cuts <- seq(from = 1, to = nrow(df2), length.out = k+1) %>% sapply(function(x) round(x, 0))

for (fold in 1:k) {
  # Exclude the observations from the fold_cuts[fold] -th index to the fold_cuts[fold+1] -th index
  df_fold <- df2[-(fold_cuts[fold]:fold_cuts[fold+1]), 2:ncol(df2)]
  left_out <- df2[(fold_cuts[fold]:fold_cuts[fold+1]), 2:ncol(df2)]
  
  gdp_reg <- lm(paste('GDP_Growth ~', regressors), df_fold)
  inflation_reg <- lm(paste('Inflation ~', regressors), df_fold)
  unemp_reg <- lm(paste('UNRATE ~', regressors), df_fold)
  gs1_reg <- lm(paste('GS1 ~', regressors), df_fold)
  
  kfold_preds[fold_cuts[fold]:fold_cuts[fold+1], 1] <- predict(gdp_reg, left_out[,5:ncol(left_out)])
  kfold_preds[fold_cuts[fold]:fold_cuts[fold+1], 2] <- predict(inflation_reg, left_out[,5:ncol(left_out)])
  kfold_preds[fold_cuts[fold]:fold_cuts[fold+1], 3] <- predict(unemp_reg, left_out[,5:ncol(left_out)])
  kfold_preds[fold_cuts[fold]:fold_cuts[fold+1], 4] <- predict(gs1_reg, left_out[,5:ncol(left_out)])
  
}

kfold_residuals <- varmat[3:nrow(varmat), ] - kfold_preds 
cor(kfold_residuals)

# Expanding Window Cross Validation  --------------------------------------

# Hyperparameter of Number of periods (i.e. predict nfor next n_periods)
n_periods = 8
min_periods = 50
xp_preds <- data.frame(matrix(ncol = 4, nrow = nrow(df2)))
names(xp_preds) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')

# Get the window intervals
#** START FROM AT LEAST 50 observations
window_cuts <- seq(from = min_periods, to = nrow(df2), by = n_periods)
if(window_cuts[length(window_cuts)] != nrow(df2)) {
  window_cuts <- c(window_cuts, nrow(df2))
}

for (i in 1:(length(window_cuts) - 1)) {
  
  # Expanding window to train VAR model on
  df_window <- df2[1:window_cuts[i], 2:ncol(df2)]
  # next n_periods periods to make predictions for
  df_predict <- df2[(window_cuts[i]+1):(window_cuts[i+1]), 2:ncol(df2)]
  
  gdp_reg <- lm(paste('GDP_Growth ~', regressors), df_window)
  inflation_reg <- lm(paste('Inflation ~', regressors), df_window)
  unemp_reg <- lm(paste('UNRATE ~', regressors), df_window)
  gs1_reg <- lm(paste('GS1 ~', regressors), df_window)
  
  xp_preds[(window_cuts[i]+1):(window_cuts[i+1]), 1] <- predict(gdp_reg, df_predict[,5:ncol(df_predict)])
  xp_preds[(window_cuts[i]+1):(window_cuts[i+1]), 2] <- predict(inflation_reg, df_predict[,5:ncol(df_predict)])
  xp_preds[(window_cuts[i]+1):(window_cuts[i+1]), 3] <- predict(unemp_reg, df_predict[,5:ncol(df_predict)])
  xp_preds[(window_cuts[i]+1):(window_cuts[i+1]), 4] <- predict(gs1_reg, df_predict[,5:ncol(df_predict)])
  print(i)
}

xp_residuals <- varmat[3:nrow(varmat), ] - xp_preds 
# Restrict xp_residuals to after first n periods
xp_residuals <- xp_residuals[(min_periods+1):nrow(xp_residuals), ]
cor(xp_residuals)


# Comparing Residuals -----------------------------------------------------

# Resiudals are: var_residuals, loo_residuals, kfold_resiudals, xp_residuals (size is 9 smaller than others)


serial_test <-serial.test(var)
plot(serial_test)

# Plotting the Squared Residuals

library(ggplot2)

par(mfrow=c(2,2))

ggplot(aes(x = 1:nrow(var_residuals), y = GDP_Growth), data = data.frame(var_residuals %>% apply(1:2, function(x) abs(x)))) +
  geom_line() +
  geom_smooth() +
  xlab('Index')

ggplot(aes(x = 1:nrow(var_residuals), y = Inflation), data = data.frame(var_residuals %>% apply(1:2, function(x) abs(x)))) +
  geom_line() +
  geom_smooth() +
  xlab('Index')

ggplot(aes(x = 1:nrow(var_residuals), y = Inflation), data = data.frame(var_residuals %>% apply(1:2, function(x) abs(x)))) +
  geom_line() +
  geom_smooth() +
  xlab('Index')

all_residuals$index = rownames(all_residuals)

residuals_long <- all_residuals %>% pivot_longer(-index, names_to = 'variable',)

ggplot() +
  geom_line(data = data.frame(var_residuals) %>% mutate(index = as.numeric(rownames(.))), aes(x = index, y = Inflation)) +
  geom_line(data = data.frame(loo_residuals) %>% mutate(index = as.numeric(rownames(.))), aes(x = index, y = Inflation)) +
  geom_line(data = data.frame(var_residuals) %>% mutate(index = as.numeric(rownames(.))), aes(x = index, y = Inflation)) +
  geom_line(data = data.frame(var_residuals) %>% mutate(index = as.numeric(rownames(.))), aes(x = index, y = Inflation)) +
  



png('Plots/In-sample GDP.png', width = 900, height = 600)
ggtsdisplay(var_residuals[,1], plot.type = 'histogram',
            main = 'In-sample GDP')
dev.off()

png('Plots/In-sample Inflation.png', width = 900, height = 600)
ggtsdisplay(var_residuals[,2], plot.type = 'histogram',
            main = 'In-sample Inflation')
dev.off()

png('Plots/In-sample Unemployment.png', width = 900, height = 600)
ggtsdisplay(var_residuals[,3], plot.type = 'histogram',
            main = 'In-sample Unemployment')
dev.off()

png('Plots/In-sample IR.png', width = 900, height = 600)
ggtsdisplay(var_residuals[,4], plot.type = 'histogram',
            main = 'In-sample IR')
dev.off()

png('Plots/Leave-one-out GDP.png', width = 900, height = 600)
ggtsdisplay(loo_residuals[,1], plot.type = 'histogram',
              main = 'Leave-one-out GDP')
dev.off()

png('Plots/Leave-one-out Inflation.png', width = 900, height = 600)
ggtsdisplay(loo_residuals[,2], plot.type = 'histogram',
            main = 'Leave-one-out Inflation')
dev.off()

png('Plots/Leave-one-out Unemployment.png', width = 900, height = 600)
ggtsdisplay(loo_residuals[,3], plot.type = 'histogram',
            main = 'Leave-one-out Unemployment')
dev.off()

png('Plots/Leave-one-out IR.png', width = 900, height = 600)
ggtsdisplay(loo_residuals[,4], plot.type = 'histogram',
            main = 'Leave-one-out IR')
dev.off()

png('Plots/k-fold CV GDP.png', width = 900, height = 600)
ggtsdisplay(kfold_residuals[,1], plot.type = 'histogram',
            main = 'k-fold CV GDP')
dev.off()

png('Plots/k-fold CV Inflation.png', width = 900, height = 600)
ggtsdisplay(kfold_residuals[,2], plot.type = 'histogram',
            main = 'k-fold CV Inflation')
dev.off()

png('Plots/k-fold CV Unemployment.png', width = 900, height = 600)
ggtsdisplay(kfold_residuals[,3], plot.type = 'histogram',
            main = 'k-fold CV Unemployment')
dev.off()

png('Plots/k-fold CV IR.png', width = 900, height = 600)
ggtsdisplay(kfold_residuals[,4], plot.type = 'histogram',
            main = 'k-fold CV IR')
dev.off()

png('Plots/Expanding Window GDP.png', width = 900, height = 600)
ggtsdisplay(xp_residuals[,1], plot.type = 'histogram',
            main = 'Expanding Window GDP')
dev.off()

png('Plots/Expanding Window Inflation.png', width = 900, height = 600)
ggtsdisplay(xp_residuals[,2], plot.type = 'histogram',
            main = 'Expanding Window Inflation')
dev.off()

png('Plots/Expanding Window Unemployment.png', width = 900, height = 600)
ggtsdisplay(xp_residuals[,3], plot.type = 'histogram',
            main = 'Expanding Window Unemployment')
dev.off()

png('Plots/Expanding Window IR.png', width = 900, height = 600)
ggtsdisplay(xp_residuals[,4], plot.type = 'histogram',
            main = 'Expanding Window IR')
dev.off()



