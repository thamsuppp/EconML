wd <- 'C:/Users/isaac/Documents/Isaac/2020-21 Penn Junior/Winter Break 2020/Economic Research PGC/VAR Coding/'
setwd(wd)



preprocessing_fn <- function(data_file) {
  
  load(paste(wd, data_file, '.Rdata', sep=''))
  
  ## split data
  N <- dim(newtrain)[1]
  Var<-dim(newtrain)[2]
  
  newtrain=as.data.frame(newtrain)
  
  # Doing a 70-30 train-test split
  train_index <- 1:round(0.7*nrow(newtrain)) # sample(1:nrow(newtrain), N*0.7)
  train_data <- newtrain[train_index,]
  test_data <- newtrain[-train_index,]
  
  form <- formula(y ~ .)
  # Create the design matrix, removing the first observation
  X <- model.matrix(form, data=train_data)[,-1]
  Y <- train_data$y
  
  Xtest <- model.matrix(form, data=test_data)[,-1]
  Ytest <- test_data$y
  
}




data_files <- c('macrotoy_gdph1', 'macrotoy_unemph1', 'macrotoy_infh1',
                'macrotoy_gdph2', 'macrotoy_unemph2', 'macrotoy_infh2')

for (data_file in data_files) {

  preprocessing_fn(data_file)
  run_models_fn(data_file)
  print(paste('Done with ', data_file))
  
}

install.packages('MacroRF')
