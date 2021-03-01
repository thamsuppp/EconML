library(tidyverse)
library(forecast)
library(corrplot)

wd <- 'C:/Users/isaac/Documents/Isaac/2020-21 Penn Junior/Winter Break 2020/Economic Research PGC/VAR Coding/'
setwd(wd)

data_path = paste(wd, 'ML Results/', sep = '')

# Load the different variables' residuals into different list objects
err_list <- list()
rsq_out_list <- list()

load(paste(data_path, 'macrotoy_gdph122feb_nowcast.Rdata', sep = ''))
err_list$gdp <- err
rsq_out_list$gdp <- Rsq_out
load(paste(data_path, 'macrotoy_unemph122feb_nowcast.Rdata', sep = ''))
err_list$unemp <- err
rsq_out_list$unemp <- Rsq_out
load(paste(data_path, 'macrotoy_infh122feb_nowcast.Rdata', sep = ''))
err_list$inf <- err
rsq_out_list$inf <- Rsq_out

err_list2 <- list()
rsq_out_list2 <- list()

load(paste(data_path, 'macrotoy_gdph1model_15feb.Rdata', sep = ''))
err_list2$gdp <- err
rsq_out_list2$gdp <- Rsq_out
load(paste(data_path, 'macrotoy_unemph1model_15feb.Rdata', sep = ''))
err_list2$unemp <- err
rsq_out_list2$unemp <- Rsq_out
load(paste(data_path, 'macrotoy_infh1model_15feb.Rdata', sep = ''))
err_list2$inf <- err
rsq_out_list2$inf <- Rsq_out




# Update this when new algorithms come
algorithm_names <- c('ARRF', 'FAARRF', 'AR', 'FAAR', 'lasso', 'ridge', 'single tree', 'random forest', 'RF tuned', 'GBM tuned')

# Getting the predictions for each algorithm
predictions <- data.frame()

for (variable in names(err_list)) {
  
  err = err_list[[variable]][1:(length(algorithm_names) + 1), ]
  # Add the y value to the errors to get the prediction
  preds = err %>% apply(2, function(x) x[1] + x) %>% .[2:(length(algorithm_names) + 1), ] %>% t() %>% data.frame()
  
  # Join the actual values to the preds dataframe
  preds = cbind(err[1, ], preds)
  colnames(preds) <- c('Actual', algorithm_names)
  preds$index <- rownames(preds)
  preds$variable <- variable
  predictions <- rbind(predictions, preds)
}

predictions_wide = predictions %>% pivot_longer(-c(index, variable),
                                    values_to = 'pred',
                                    names_to = 'algorithm') %>%
  data.frame()

predictions_wide$algorithm <- factor(predictions_wide$algorithm)
predictions_wide$variable <- factor(predictions_wide$variable)
predictions_wide$index <- as.numeric(predictions_wide$index)

# Plot the out-of-sample predictions
png('Plots/ML/Predictions Nowcast 22feb.png', width = 1800, height = 1200, res = 200)
ggplot(data = predictions_wide, aes(x = index, y = pred, color = algorithm)) +
  geom_line() +
  ggtitle('OOS Predictions for different algorithms') +
  facet_grid(rows = vars(variable), scales = 'free') +
  scale_color_manual(values = c('black', 'darkred', 'blue', 'green3', 'yellow3', 'orange2', 'slategrey', 'saddlebrown', 'mediumorchid3', 'lightcoral', 'turquoise'))
dev.off()


# Investigating the out of sample residuals of the ML algorithms
# err, Rsq_out is the input to this script

residuals <- data.frame()

# Get all the residuals for different variables into one dataframe
for (variable in names(err_list)) {
  err = err_list[[variable]]
  err = err[2:(length(algorithm_names) + 1), ] %>% t() %>% data.frame()
  colnames(err) <- algorithm_names
  err$index <- rownames(err)
  err$variable <- variable
  residuals <- rbind(residuals, err)
}

residuals_wide <- residuals %>% pivot_longer(-c(index, variable),
                                             values_to = 'residual',
                                             names_to = 'algorithm') %>%
  data.frame()

residuals_wide$algorithm <- factor(residuals_wide$algorithm)
residuals_wide$variable <- factor(residuals_wide$variable)
residuals_wide$index <- as.numeric(residuals_wide$index)


# Comparing the residuals over time 
png('Plots/ML/Residuals Nowcast 22feb.png', width = 1800, height = 1200, res = 200)
ggplot(data = residuals_wide, aes(x = index, y = residual, color = algorithm)) +
  geom_line() +
  ggtitle('OOS Residuals for different algorithms') +
  facet_grid(rows = vars(variable), scales = 'free') +
  scale_color_manual(values = c('darkred', 'blue', 'green3', 'yellow3', 'orange2', 'slategrey', 'saddlebrown', 'mediumorchid3', 'lightcoral', 'turquoise'))
dev.off()

# Looking at the Rsquared
rsq_out_df = data.frame(gdph1)



# Getting the cross correlation
for (var in names(err_list)) {
  residuals_var = residuals %>% filter(variable == var) %>% .[,1:10]
  print(cor(residuals_var))
  par(mfrow = c(1,1))
  corrplot(cor(residuals_var), method = 'number', tl.col = 'black', main = paste('Correlation of Residuals Across Algos - ', var))


  # Persistence of residuals (ACF)
  # par(mfrow = c(2,3))
  # 
  # # Get all except ridge and single tree
  # for (col in c(1,2,3,6,7,8)) {
  #   acf(residuals_var[, col], main = paste(colnames(residuals)[col], var))
  # }
}

# Getting the OOS MSE 
mse_list <- rep(0, length(algorithm_names))
names(mse_list) = algorithm_names
mae_list <- rep(0, length(algorithm_names))
names(mae_list) = algorithm_names

for (algo in algorithm_names) {
  mse_list[algo] = sapply( (predictions[algo] - predictions['Actual']) ^ 2, mean)
  mae_list[algo] = sapply( (predictions[algo] - predictions['Actual']), function(x) mean(abs(x)))

}





