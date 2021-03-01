library(tidyverse)
library(xts)
library(tseries)
library(vars)
library(dynlm)
library(forecast)
library(TTR)
library(lindia)

wd = 'C:/Users/isaac/Documents/Isaac/2020-21 Penn Junior/Winter Break 2020/Economic Research PGC/VAR Coding/'
setwd(wd)
plot_path = paste(wd, 'Plots/OOS Investigation/', sep = '')

# PARAMETERS
var_p = 4
# For block cross validation
block_size = 8
# For rolling window
window_size = 100
n_periods = 8


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

year_index <- seq(from = 1953.5, to = 2020, by = 0.25)

# Initial VAR -------------------------------------------------------------

# Getting the data to input to VAR
input_mat <- as.matrix(df[, 2:ncol(df)])

# Model selection of VAR lag
var_selection <- VARselect(input_mat, lag.max = 12, type = 'const')
var_selection$selection

# Fitting the VAR
var <- VAR(input_mat, p = var_p, type = 'const', season = NULL, exog = NULL)
summary(var)

var_preds <- predict(var, input_mat)
var_residuals <- data.frame(residuals(var))
cor(var_residuals)

# Get the processed VAR matrix - to do OLS on
var_mat <- var$datamat

var_mat_year_index <- year_index[4:length(year_index)]

# Dynamically create the regressors for var_mat
regressors <- colnames(var_mat) %>% .[5:(length(.)-1)] %>% paste(., collapse = ' + ')

# In-sample VAR via OLS ---------------------------------------------------

# Testing to see if OLS is equivalent to VAR function
preds <- data.frame(matrix(ncol = 4, nrow = nrow(var_mat)))
names(preds) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')

ols_residuals <- data.frame(matrix(ncol = 4, nrow = nrow(var_mat)))
names(ols_residuals) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')

gdp_reg <- lm(paste('GDP_Growth ~', regressors), var_mat)
inflation_reg <- lm(paste('Inflation ~', regressors), var_mat)
unemp_reg <- lm(paste('UNRATE ~', regressors), var_mat)
gs1_reg <- lm(paste('GS1 ~', regressors), var_mat)

preds[, 1] <- predict(gdp_reg, var_mat[,5:ncol(var_mat)])
preds[, 2] <- predict(inflation_reg, var_mat[,5:ncol(var_mat)])
preds[, 3] <- predict(unemp_reg, var_mat[,5:ncol(var_mat)])
preds[, 4] <- predict(gs1_reg, var_mat[,5:ncol(var_mat)])

ols_residuals <-var_mat[,1:4] - preds
cor(ols_residuals)

# Leave-one-out VAR -------------------------------------------------------
# Wrapper function for LOO VAR - (generalized for lags, not generalized for variables)
leave_one_out_fn <- function(var_mat) {
  
  loo_preds <- data.frame(matrix(ncol = 4, nrow = nrow(var_mat)))
  names(loo_preds) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')
  
  for (i in 1:nrow(var_mat)) {
    
    # Leave one out of the matrix and estimate the VAR regression
    df_loo <- var_mat[-i, ]
    left_out <- var_mat[i, ]
    
    # Run the regression (exclude the other contemporaneous variables)
    gdp_reg <- lm(paste('GDP_Growth ~', regressors), df_loo)
    inflation_reg  <- lm(paste('Inflation ~', regressors), data = df_loo)
    unemp_reg  <- lm(paste('UNRATE ~', regressors), data = df_loo)
    gs1_reg  <- lm(paste('GS1 ~', regressors), data = df_loo)
    
    # Get out of sample predictions for the left out observation
    
    loo_preds[i, 1] <- predict(gdp_reg, left_out[5:ncol(left_out)])
    loo_preds[i, 2] <- predict(inflation_reg, left_out[5:ncol(left_out)])
    loo_preds[i, 3] <- predict(unemp_reg, left_out[5:ncol(left_out)])
    loo_preds[i, 4] <- predict(gs1_reg, left_out[5:ncol(left_out)])
  }
  
  loo_residuals <- var_mat[,1:4]  - loo_preds
  return(loo_residuals)
}

# K-fold Cross Validation -------------------------------------------------
kfold_fn <- function(var_mat, block_size = 8) {
  
  kfold_preds <- data.frame(matrix(ncol = 4, nrow = nrow(var_mat)))
  names(kfold_preds) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')
  # Get the fold cuts (some blocks may be block_size+1 due to rounding)
  fold_cuts <- seq(from = 1, to = nrow(var_mat), length.out = ceiling(nrow(var_mat) / block_size)) %>% sapply(function(x) round(x, 0))
  
  for (fold in 1:(length(fold_cuts)-1)) {
    # Exclude the observations from the fold_cuts[fold] -th index to the fold_cuts[fold+1] -th index
    df_fold <- var_mat[-(fold_cuts[fold]:fold_cuts[fold+1]), ]
    left_out <- var_mat[(fold_cuts[fold]:fold_cuts[fold+1]), ]
    
    gdp_reg <- lm(paste('GDP_Growth ~', regressors), df_fold)
    inflation_reg <- lm(paste('Inflation ~', regressors), df_fold)
    unemp_reg <- lm(paste('UNRATE ~', regressors), df_fold)
    gs1_reg <- lm(paste('GS1 ~', regressors), df_fold)
    
    kfold_preds[fold_cuts[fold]:fold_cuts[fold+1], 1] <- predict(gdp_reg, left_out[,5:ncol(left_out)])
    kfold_preds[fold_cuts[fold]:fold_cuts[fold+1], 2] <- predict(inflation_reg, left_out[,5:ncol(left_out)])
    kfold_preds[fold_cuts[fold]:fold_cuts[fold+1], 3] <- predict(unemp_reg, left_out[,5:ncol(left_out)])
    kfold_preds[fold_cuts[fold]:fold_cuts[fold+1], 4] <- predict(gs1_reg, left_out[,5:ncol(left_out)])
    
  }
  
  kfold_residuals <- var_mat[,1:4] - kfold_preds 
  return(kfold_residuals)
}

# Rolling Window Cross Validation -----------------------------------------
rolling_window_fn <- function(Var_mat, window_size = 100, n_periods = 8) {
  
  roll_preds <- data.frame(matrix(ncol = 4, nrow = nrow(var_mat)))
  names(roll_preds) <- c('GDP_Growth', 'Inflation', 'UNRATE', 'GS1')
  
  window_cuts <- seq(from = window_size, to = nrow(var_mat), by = n_periods)
  if(window_cuts[length(window_cuts)] != nrow(var_mat)) {
    window_cuts <- c(window_cuts, nrow(var_mat))
  }
  
  for (i in 1:(length(window_cuts) - 1)) {
    
    window_start_id <- window_cuts[i]-window_size+1
    
    # Rolling window to train VAR model on
    df_window <- var_mat[window_start_id:window_cuts[i], ]
    # next _ periods to make predictions for
    df_predict <- var_mat[(window_cuts[i]+1):(window_cuts[i+1]), ]
    
    gdp_reg <- lm(paste('GDP_Growth ~', regressors), df_window)
    inflation_reg <- lm(paste('Inflation ~', regressors), df_window)
    unemp_reg <- lm(paste('UNRATE ~', regressors), df_window)
    gs1_reg <- lm(paste('GS1 ~', regressors), df_window)
    
    roll_preds[(window_cuts[i]+1):(window_cuts[i+1]), 1] <- predict(gdp_reg, df_predict[,5:ncol(df_predict)])
    roll_preds[(window_cuts[i]+1):(window_cuts[i+1]), 2] <- predict(inflation_reg, df_predict[,5:ncol(df_predict)])
    roll_preds[(window_cuts[i]+1):(window_cuts[i+1]), 3] <- predict(unemp_reg, df_predict[,5:ncol(df_predict)])
    roll_preds[(window_cuts[i]+1):(window_cuts[i+1]), 4] <- predict(gs1_reg, df_predict[,5:ncol(df_predict)])
  
  }
  
  # Note: There are no predictions for the first 100 due to the rolling_window
  # Index of roll_residuals ar 101:263
  roll_residuals <- var_mat[(window_size+1):nrow(var_mat), 1:4] - roll_preds[(window_size+1):nrow(roll_preds), ]
  return(roll_residuals)
}
  

# Processing Residuals -----------------------------------------------------
# Residuals are: var_residuals (=ols_residuals), loo_residuals, kfold_residuals, roll_residuals (size is 100 smaller than others)


# EXECUTION CODE
# loo_residuals <- leave_one_out_fn(var_mat)
# kfold_residuals <- kfold_fn(var_mat, block_size = 8)
# roll_residuals <- rolling_window_fn(var_mat, window_size = 100, n_periods = 8)

# Combine residuals together into a list
# residuals_list = list(var_residuals, loo_residuals, kfold_residuals, roll_residuals)
# names(residuals_list) = c('All', 'Leave One Out', 'Block', 'Rolling Window')

residuals_list = list(kfold_fn(var_mat, block_size = 4),
                      kfold_fn(var_mat, block_size = 8),
                      kfold_fn(var_mat, block_size = 12),
                      kfold_fn(var_mat, block_size = 16),
                      kfold_fn(var_mat, block_size = 20))
names(residuals_list) <- c(4, 8, 12, 16, 20)

# Combine residuals into one df, adding the type and index columns
all_residuals = data.frame()
for (exp_type in names(residuals_list)) {
  all_residuals = all_residuals %>% rbind(residuals_list[[exp_type]] %>% mutate(type = exp_type, index = as.numeric(rownames(.))))
}

residuals_long <- all_residuals %>% pivot_longer(-c('type', 'index'), names_to = 'variable', values_to = 'residual') %>% data.frame()
# Get the Year column from Index
residuals_long$year <- sapply(residuals_long$index, function(x) var_mat_year_index[x])

# Set the factor levels
residuals_long$type <- factor(residuals_long$type, levels = c(4,8,12,16,20))

# Plots of Residuals over Time --------------------------------------------
# Plots: Residuals, Abs Residuals, Abs Residuals with SMA, Sq Residuals, Sq Residuals with SMA

# Comparison plot of the residuals over time
png(paste(plot_path, 'Residuals.png',sep = ''), width = 1800, height = 1200, res = 200)
ggplot(data = residuals_long, aes(x = year, y = residual, color = type)) +
  geom_line() +
  ggtitle('OLS Residuals for Different OOS Experiments') +
  facet_grid(rows = vars(variable), scales = 'free') +
  scale_x_continuous(breaks = seq(1955, 2020, by = 5))
dev.off()

# Abs of residuals over time
png(paste(plot_path, 'Abs Residuals.png',sep = ''), width = 1800, height = 1200, res = 200)
ggplot(data = residuals_long, aes(x = year, y = abs(residual), color = type)) +
  geom_line() +
  ggtitle('OLS Abs Residuals for Different OOS Experiments') +
  facet_grid(rows = vars(variable), scales = 'free') +
  scale_x_continuous(breaks = seq(1955, 2020, by = 5))
dev.off()

# SMA of square of residuals over time 
residuals_long <- residuals_long %>% group_by(type, variable) %>% mutate(sma = SMA(abs(residual), n = 20)) %>% data.frame()
png(paste(plot_path, 'Abs Residuals SMA.png',sep = ''), width = 1800, height = 1200, res = 200)
ggplot(data = residuals_long, aes(x = year, y = sma, color = type)) +
  geom_line() +
  ggtitle('OLS Square Residuals for Different OOS Experiments') +
  facet_grid(rows = vars(variable), scales = 'free') +
  scale_x_continuous(breaks = seq(1955, 2020, by = 5))
dev.off()

# Square of residuals over time - Variance
png(paste(plot_path, 'Squared Residuals.png',sep = ''), width = 1800, height = 1200, res = 200)
ggplot(data = residuals_long, aes(x = year, y = sapply(residual, function(x) x^2), color = type)) +
  geom_line() +
  ggtitle('OLS Squared Residuals for Different OOS Experiments') +
  facet_grid(rows = vars(variable), scales = 'free') +
  scale_x_continuous(breaks = seq(1955, 2020, by = 5))
dev.off()

residuals_long <- residuals_long %>% group_by(type, variable) %>% mutate(sma = SMA(residual ^2, n = 20)) %>% data.frame()
png(paste(plot_path, 'Squared Residuals SMA.png',sep = ''), width = 1800, height = 1200, res = 200)
ggplot(data = residuals_long, aes(x = year, y = sma, color = type)) +
  geom_line() +
  ggtitle('OLS Square Residuals for Different OOS Experiments') +
  facet_grid(rows = vars(variable), scales = 'free') +
  scale_x_continuous(breaks = seq(1955, 2020, by = 5))
dev.off()
  


# Residual Histograms -----------------------------------------------------

# Function to Plot the residuals' histogram
residual_histogram <- function(g, var_name, exp_name) {
  h = hist(g, freq = FALSE, breaks = 10, xlab = NULL, main = paste(var_name, ' - ', exp_name))
  lines(density(g), col = 'blue', lwd = 2)
  x = seq(min(g), max(g), length = 100)
  y = dnorm(x, mean = mean(g), sd = sd(g))
  lines(x, y, col = 'black', lwd = 2)
}

par(mfrow = c(2,2), mai = c(0.4, 0.7, 0.5, 0.3))

# Draw the residual histogram plot for all the experiments
for (i in 1:ncol(residuals_list[[exp_type]])) {
  for (exp_type in names(residuals_list)) {
    residual_histogram(residuals_list[[exp_type]][,i], colnames(residuals_list[[exp_type]])[i], exp_type)
  }
}

# Rolling Cross-Correlations ----------------------------------------------
# (Not generalized for variable)

# Rolling window to calculate cross-correlations with
corr_window = 20

# Function to calculate rollign correlations
# Input : residuals_long, Outpu: corr_df_long
rolling_corr_fn <- function(residuals_long, corr_window = 40) {

  # Rolling correlation matrix
  corr_df = data.frame(matrix(ncol = 8, nrow = 0))
  colnames(corr_df) = c('type', 'index', 'GDP-Inf', 'GDP-Unemp', 'GDP-GS1', 'Inf-Unemp', 'Inf-GS1', 'Unemp-GS1')
  
  for (i in corr_window:max(residuals_long$index)) {
    for (exp_type in unique(residuals_long$type)) {
      window_df = residuals_long %>% filter(index <= i & index > (i - corr_window) & type == exp_type)
      
      if (nrow(window_df) == corr_window * 4) {
        # Pivot wider
        window_df_wide = window_df[,c('index', 'variable', 'residual')] %>% pivot_wider(names_from = variable,
                                                                            values_from = residual) %>% data.frame()
        # Get correlation matrix - assign it to the corr_arr
        corr_mat = cor(window_df_wide[, 2:ncol(window_df_wide)])
        
        out = data.frame(exp_type, i, corr_mat[1,2], corr_mat[1,3], corr_mat[1,4], corr_mat[2,3], corr_mat[2,4], corr_mat[3,4])
        colnames(out) = colnames(corr_df)
        corr_df = rbind(corr_df, out)
      
      }
    }
  }
  # Plot GDP-inflation (1,2), GDP-UNRATE (1,3), GDP-GS1 (1,4), Inflation-UNRATE (2,3), Inflation-GS1 (2,4), UNRATE-GS1 (3,4)
  corr_df_long = corr_df %>% pivot_longer(-c('type', 'index'), names_to = 'pair', values_to = 'correlation') %>% drop_na()
  corr_df_long$year <- sapply(corr_df_long$index, function(x) var_mat_year_index[x])
  
  return(corr_df_long)

}

corr_df_long <- rolling_corr_fn(residuals_long, corr_window = 40)

png(paste(plot_path, 'Rolling CrossCorrelation varying Blocksize.png', sep = ''), width = 1800, height = 1200, res = 200)
ggplot(corr_df_long, aes(x = year, y = correlation, col = type)) +
  geom_line() +
  facet_wrap(~ pair, ncol = 2) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  scale_x_continuous(breaks = seq(1955, 2020, by = 5)) +
  ggtitle('Rolling Cross-Correlations Varying Block Size')
dev.off()



