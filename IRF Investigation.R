wd <- 'C:/Users/isaac/Documents/Isaac/2020-21 Penn Junior/Winter Break 2020/Economic Research PGC/EconML/'
setwd(wd)

library(tidyverse)
library(vars)
library(lpirfs)

# Getting Residuals and Variables -----------------------------------------

data_path = paste(wd, 'ML Results/', sep = '')
experiment_name = 'oob'

# Load the different variables' residuals into different list objects
# Changed from err to err_all
err_list <- list()
rsq_out_list <- list()
load(paste(data_path, 'macrotoy_gdph11mar_oob.Rdata', sep = ''))
err_list$gdp <- err_all
rsq_out_list$gdp <- Rsq_out
load(paste(data_path, 'macrotoy_unemph11mar_oob.Rdata', sep = ''))
err_list$unemp <- err_all
rsq_out_list$unemp <- Rsq_out
load(paste(data_path, 'macrotoy_infh11mar_oob.Rdata', sep = ''))
err_list$inf <- err_all
rsq_out_list$inf <- Rsq_out

# Update this when new algorithms come
algorithm_names <- c('ARRF', 'FAARRF', 'AR', 'FAAR', 'lasso', 'ridge', 'single tree', 'random forest', 'RF tuned', 'GBM tuned')

# Get the residuals
residuals <- data.frame()

# Get all the residuals for different variables into one dataframe
for (variable in names(err_list)) {
  error = err_list[[variable]]
  error = error[2:(length(algorithm_names) + 1), ] %>% t() %>% data.frame()
  colnames(error) <- algorithm_names
  error$index <- rownames(error)
  error$variable <- variable
  residuals <- rbind(residuals, error)
}

# Take negative of residuals to get the true residual (stop gap before I change the actual residual code!)
residuals[, 1:10] = -residuals[, 1:10]

### Get the actual out-of-sample y values

data_files = c('macrotoy_gdph1', 'macrotoy_unemph1', 'macrotoy_infh1')
oos_y = data.frame(index = 1:ncol(err))
all_y = data.frame(index = 1:ncol(err_all))

for (data_file in data_files) {
  load(paste(wd, data_file, '.Rdata', sep=''))
  N <- dim(newtrain)[1]
  Var<-dim(newtrain)[2]
  all_data =as.data.frame(newtrain)
  rm(newtrain)
  
  # Doing a 70-30 train-test split
  train_index <- 1:round(0.7*nrow(all_data)) # sample(1:nrow(newtrain), N*0.7)
  train_data <- all_data[train_index,]
  test_data <- all_data[-train_index,]
  form <- formula(y ~ .)
  Xtest <- model.matrix(form, data=test_data)[,-1]
  Ytest <- test_data$y
  oos_y[data_file] = Ytest
  all_y[data_file] = all_data$y
  
}
colnames(oos_y) = c('index', 'gdp', 'unemp', 'inf')


# Local Projections -------------------------------------------------------

# What we need: oos_y (from above) - or all_y (if using OOB data)

# Legnth of irf_df: 1080 = 3 shock_vars x 10 algos x 12 horizons x 3 estimates
irf_df = data.frame()
var_names = c('gdp', 'unemp', 'inf')

# Option  to include all the variables as controls for the local projection or not
include_all_controls = T

if (include_all_controls == T) {

  # Version where you include all variables as controls (in endog_data)
  for (shock_var in var_names) {
    residuals_var = residuals %>% filter(variable == shock_var)

    # Loop through the residuals for all the algos
    for (algo in colnames(residuals_var)[1:10]) {
      
      shock = data.frame(residuals_var[algo])
      endog_data = all_y[, 2:4] # USED TO BE oos_y
      
      # Estimate the local projection linear model
      results_lin_iv = lp_lin_iv(endog_data = endog_data,
                                 lags_endog_lin = 2, ### HYPERPARAM TO CHANGE
                                 shock = shock,
                                 trend = 0,
                                 confint = 1.96,
                                 hor = 12) ### HYPERPARAM TO CHANGE
      
      mean = results_lin_iv$irf_lin_mean %>% t() %>% data.frame() %>% mutate(estimate = 'mean', index = as.numeric(rownames(.)))
      lcl = results_lin_iv$irf_lin_low %>% t() %>% data.frame() %>% mutate(estimate = 'lower', index = as.numeric(rownames(.)))
      ucl = results_lin_iv$irf_lin_up %>% t() %>% data.frame() %>% mutate(estimate = 'upper', index = as.numeric(rownames(.)))
      irf_var_algo_df = rbind(mean, lcl, ucl)
      colnames(irf_var_algo_df) = c(var_names, 'estimate', 'index')
      irf_var_algo_df$algo = algo
      irf_var_algo_df$shock_var = shock_var
      irf_df = rbind(irf_df, irf_var_algo_df)
    }
  }
  
  # Length: 3240 = 3 shock_vars x 10 algos x 12 horizons x 3 target variables x 3 estimates
  irf_df_long = irf_df %>% pivot_longer(-c('index', 'algo', 'shock_var', 'estimate'),
                                        names_to = 'target_var',
                                        values_to = 'irf') %>%
    pivot_wider(names_from = 'estimate',
                values_from = 'irf') %>% data.frame()
  
} else {
  # Version where you only include the shock and target variable
  
  for (shock_var in var_names) {
    residuals_var = residuals %>% filter(variable == shock_var)

    # Loop through the residuals for all the algos
    for (algo in colnames(residuals_var)[1:10]) {
      for (target_var in var_names) {
        
        shock = data.frame(residuals_var[algo])
        target = data.frame(all_y[, target_var]) # USED TO BE oos_y
        
        # Estimate the local projection linear model
        results_lin_iv = lp_lin_iv(endog_data = target,
                                   lags_endog_lin = 2, ### HYPERPARAM TO CHANGE
                                   shock = shock,
                                   trend = 0,
                                   confint = 1.96,
                                   hor = 12) ### HYPERPARAM TO CHANGE
        
        results_lin_iv = list(results_lin_iv$irf_lin_low,
                              results_lin_iv$irf_lin_mean,
                              results_lin_iv$irf_lin_up)
        irf_out_names = c('lower', 'mean', 'upper')
        names(results_lin_iv) = irf_out_names
        
        results_lin_iv  = results_lin_iv %>% lapply(function(x) {x %>% t() %>% data.frame()}) %>%
          map2(irf_out_names, function(x, y) {
            x$index = as.numeric(rownames(x));
            colnames(x) = c('irf', 'index'); x$estimate = y; return(x)}) %>%
          bind_rows()
        results_lin_iv$target_var = target_var
        results_lin_iv$algo = algo
        results_lin_iv$shock_var = shock_var
        
        irf_df = rbind(irf_df, results_lin_iv)
      }
    }
  }
  
  # Length: 3240 = 3 shock_vars x 10 algos x 12 horizons x 3 target variables x 3 estimates
  irf_df_long = irf_df %>%
    pivot_wider(names_from = 'estimate',
                values_from = 'irf') %>% data.frame()
}


# Plotting ----------------------------------------------------------------

# Plots with the confidence intervals

# Individual Plots
# for (shock_variable in var_names) {
# 
#   g = irf_df_long %>% filter(shock_var == shock_variable,
#                          algo %in% c('ARRF', 'FAARRF', 'AR', 'FAAR', 'RF tuned', 'GBM tuned'))  %>%
#     ggplot(aes(x = index, y = mean, color = algo)) +
#       geom_ribbon(aes(ymin = lower, ymax = upper, fill = algo), linetype = 'dotted', alpha = 0.1) +
#       geom_line(size = 1) +
#       facet_wrap (~target_var, scales = 'free', dir = 'v') +
#       scale_color_manual(values = c('red', 'blue', 'green3', 'yellow3', 'orange', 'slategrey', 'saddlebrown')) +
#       scale_fill_manual(values = c('red', 'blue', 'green3', 'yellow3', 'orange', 'slategrey', 'saddlebrown')) +
#       ggtitle(paste('IRF for shock to', shock_variable))
# 
#   print(g)
# 
# }

irf_df_long$shock_var = paste(irf_df_long$shock_var, 'shock')
irf_df_long$target_var = paste(irf_df_long$target_var, 'target')

png(paste('Plots/ML/IRF_', experiment_name, '_withCI.png', sep = ''), width = 1800, height = 1500, res = 200)
irf_df_long %>% filter(algo %in% c('ARRF', 'FAARRF', 'AR', 'FAAR', 'RF tuned', 'GBM tuned'))  %>%
  ggplot(aes(x = index, y = mean, color = algo)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = algo), linetype = 'dotted', alpha = 0.1) +
  geom_line() +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = 'dashed') +
  facet_wrap (~shock_var + target_var, scales = 'free') +
  scale_color_manual(values = c('red', 'blue', 'green3', 'yellow3', 'orange2', 'slategrey', 'saddlebrown')) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) +
  ggtitle('Impulse Response Functions (Local Projection)')
dev.off()


# Plots without the confidence intervals

# Individual Plots
# 
# for (shock_variable in var_names) {
#   
#   g = irf_df_long %>% filter(shock_var == shock_variable, 
#                              algo %in% c('ARRF', 'FAARRF', 'AR', 'FAAR', 'RF tuned', 'GBM tuned'))  %>%
#     ggplot(aes(x = index, y = mean, color = algo)) +
#     geom_line(size = 1) +
#     facet_wrap (~target_var, scales = 'free', dir = 'v') +
#     scale_color_manual(values = c('red', 'blue', 'green3', 'yellow3', 'orange', 'slategrey', 'saddlebrown')) +
#     ggtitle(paste('IRF for shock to', shock_variable))
#   
#   print(g)
#   
# }

png(paste('Plots/ML/IRF_', experiment_name, '_noCI.png', sep = ''), width = 1800, height = 1500, res = 200)
irf_df_long %>% filter(algo %in% c('ARRF', 'FAARRF', 'AR', 'FAAR', 'RF tuned', 'GBM tuned'))  %>%
  ggplot(aes(x = index, y = mean, color = algo)) +
  geom_line() +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = 'dashed') +
  facet_wrap (~shock_var + target_var, scales = 'free') +
  scale_color_manual(values = c('red', 'blue', 'green3', 'yellow3', 'orange2', 'slategrey', 'saddlebrown')) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) +
  ggtitle('Impulse Response Functions (Local Projection)')
dev.off()



# Generating IRFs from VAR ------------------------------------------------


# Length of irf_df: 1080 = 3 shock_vars x 10 algos x 12 horizons x 3 estimates
var_irf_df = data.frame()
var_names = c('gdp', 'unemp', 'inf')

for (shock_var in var_names) {
  
  residuals_var = residuals %>% filter(variable == shock_var)

  # Loop through the residuals for all the algos
  for (algo in c('ARRF', 'FAARRF', 'AR', 'FAAR', 'RF tuned', 'GBM tuned')) {

    shock = data.frame(residuals_var[algo])
    endog_data = all_y[, 2:4] # USED TO BE oos_y
    
    var_data = cbind(shock, endog_data)
    colnames(var_data) = c('shock', var_names)
    
    # Estimate VAR model
    var_model = VAR(var_data, p = 2, type = 'const')
    #plot(var_model)
    
    var_irf_algo_df = data.frame()
    
    irf_out = irf(var_model, impulse = 'shock', 
                  response = var_names, 
                  ortho = T,
                  n.ahead = 12, runs = 100)
    irf_out_names = c('lower', 'mean', 'upper')
    irf_out = list(irf_out$Lower$shock, irf_out$irf$shock, irf_out$Upper$shock)
    names(irf_out) = irf_out_names
    
    var_irf_algo_df = irf_out %>% lapply(data.frame) %>% 
      # Remove first result (since it's 0)
      lapply(function(x) {x = x[2:nrow(x), ]; rownames(x) = NULL; return(x)}) %>%
      # Mutate another column storing n period ahead
      lapply(function(x) {x$index = as.numeric(rownames(x)); return(x)}) %>%
      # Mutate another column storing name of the estimate
      map2(irf_out_names, function(x, y) {x$estimate = y; return(x)}) %>%
      # Combine all the dataframes together into one dataframe
      bind_rows()
    
    var_irf_algo_df$algo = algo
    var_irf_algo_df$shock_var = shock_var
    var_irf_df = rbind(var_irf_df, var_irf_algo_df)
  }
}

# var_irf_df and irf_df are directly comparable

# Length: 3240 = 3 shock_vars x 10 algos x 12 horizons x 3 target variables x 3 estimates
var_irf_df_long = var_irf_df %>% pivot_longer(-c('index', 'algo', 'shock_var', 'estimate'),
                                      names_to = 'target_var',
                                      values_to = 'irf') %>%
  pivot_wider(names_from = 'estimate',
              values_from = 'irf') %>% data.frame()


var_irf_df_long$shock_var = paste(var_irf_df_long$shock_var, 'shock')
var_irf_df_long$target_var = paste(var_irf_df_long$target_var, 'target')

png(paste('Plots/ML/IRF_VAR_', experiment_name, '_withCI.png', sep = ''), width = 1800, height = 1500, res = 200)
var_irf_df_long %>% filter(algo %in% c('ARRF', 'FAARRF', 'AR', 'FAAR', 'RF tuned', 'GBM tuned'))  %>%
  ggplot(aes(x = index, y = mean, color = algo)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = algo), linetype = 'dotted', alpha = 0.1) +
  geom_line() +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = 'dashed') +
  facet_wrap (~shock_var + target_var, scales = 'free') +
  scale_color_manual(values = c('red', 'blue', 'green3', 'yellow3', 'orange2', 'slategrey', 'saddlebrown')) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) +
  ggtitle('Impulse Response Functions (VAR)')
dev.off()


png(paste('Plots/ML/IRF_VAR_', experiment_name, '_noCI.png', sep = ''), width = 1800, height = 1500, res = 200)
var_irf_df_long %>% filter(algo %in% c('ARRF', 'FAARRF', 'AR', 'FAAR', 'RF tuned', 'GBM tuned'))  %>%
  ggplot(aes(x = index, y = mean, color = algo)) +
  geom_line() +
  geom_hline(aes(yintercept = 0), colour = "black", linetype = 'dashed') +
  facet_wrap (~shock_var + target_var, scales = 'free') +
  scale_color_manual(values = c('red', 'blue', 'green3', 'yellow3', 'orange2', 'slategrey', 'saddlebrown')) +
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12)) +
  ggtitle('Impulse Response Functions (VAR)')
dev.off()




