

wd <- 'C:/Users/isaac/Documents/Isaac/2020-21 Penn Junior/Winter Break 2020/Economic Research PGC/VAR Coding/'
setwd(wd)

## Function Libraries
library(glmnet)
library(rpart)
library(ranger)
library(mboost)
library(gbm)
library(EZtune)
library(ISLR)
library(caret)
library(dplyr)
library(tgp)
library(knockoff)
library(earth)
library(tuneRanger)
library(party)
library(readxl)
library(pracma)
library(tidyverse)
library(mlbench)
library(AppliedPredictiveModeling)
library(MacroRF)


# R2 Helper Functions -----------------------------------------------------


# Function to calculate R2
r.squared = function(pred, actual) {
  rss <- sum((actual-pred)^2)
  tss <- sum((actual-mean(actual))^2)
  rsq=1-rss/tss
  #if(rsq< -0.5){rsq=-0.5}
  return(rsq)
}

r.squared.fair = function(pred, actual,y.is) {
  rss <- sum((actual-pred)^2)
  tss <- sum((actual-mean(y.is))^2)
  rsq=1-rss/tss
  #if(rsq< -0.5){rsq=-0.5}
  return(rsq)
}

r.squared.fair.2 = function(pred, actual,y.is,na.rm=TRUE) {
  if(sum(is.na(pred))>0){print(paste('There were',sum(is.na(pred)),'NAs',sep=''))}
  rss <- mean((actual-pred)^2,na.rm=na.rm)
  tss <- mean((actual-mean(y.is))^2,na.rm=na.rm)
  rsq=1-rss/tss
  #if(rsq< -0.5){rsq=-0.5}
  return(rsq)
}


# Function to Run Models --------------------------------------------------

data_file = 'macrotoy_gdph1'
nowcast = F

# Run Basic machinery.R to get X, Y, Xtest, Ytest

run_models_fn <- function(data_file, nowcast = F, save_name = 'model') {


  ###########  Preprocessing #####################
  
  # Preprocessing
  load(paste(wd, data_file, '.Rdata', sep=''))
  
  all_data = as.data.frame(newtrain)
  rm(newtrain)
  # If nowcasting instead of forecasting, then replace y with L_0y and delete L_0y
  if (nowcast == T) {
    all_data$y = all_data$L_0y
    all_data = all_data[,!(names(all_data) %in% c('L_0y'))]
    
  }
  
  ## split data
  N <- dim(all_data)[1]
  Var<-dim(all_data)[2]
  
  # Doing a 70-30 train-test split
  train_index <- 1:round(0.7*nrow(all_data)) # sample(1:nrow(all_data), N*0.7)
  train_data <- all_data[train_index,]
  test_data <- all_data[-train_index,]
  
  form <- formula(y ~ .)
  # Create the design matrix, removing the first observation
  Xtrain <- model.matrix(form, data=train_data)[,-1]
  Ytrain <- train_data$y
  
  Xtest <- model.matrix(form, data=test_data)[,-1]
  Ytest <- test_data$y
  
  dataset = substr(data_file, start = 10, stop = nchar(data_file))

  set.seed(2020)
  p <- data.frame()
  Rsq_in <- data.frame()
  Rsq_out <- data.frame()
  Roos_list=list()
  Rin_list=list()
  
  dsds = strtoi(Sys.getenv("V")) #GLOBAL FROM THE HPCC DECIDING WHICH DATASET FOR WHICH CORE
  if(is.na(dsds)){dsds=2}
  DS=dsds
  #for(DS in dsds){ #21:22,19:20,12:15,1:4,16:17 unique
  

  
  ##----------------------------------------------------------------------------------------------------------------------##
  # Models
  ##----------------------------------------------------------------------------------------------------------------------##
  
  # Instantiate all errors (including in-sample)
  all_err = array()
  
  # Instantiate error matrix
  err=array(NA,dim=c(65,length(Ytest)))
  m=1
  #err[m,1,1:length(Y)]=Y
  
  err[m,1:length(Ytest)]=Ytest
  
  print(dsds)
  print(dim(train_data))
  
  
  ########### MACRO RANDOM FOREST #####################
  
  datamat = data.matrix(all_data)
  
  # Include both the training and test data - hence use all_data
  # ARRF: y.pos = 1, x.pos = 2:7, s.pos = 50:600
  arrf = MRF(datamat, y.pos = 1, x.pos = 2:7, S.pos = 2:ncol(datamat), 
             oos.pos = (round(0.7*nrow(datamat))+1):nrow(datamat) , B = 50)
  
  reg.out.error = arrf$pred - Ytest
  Rsq_out[dataset, "ARRF"] <- r.squared.fair(arrf$pred, Ytest, Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=reg.out.error
  
  print('Done with ARRF')
  print(Sys.time())
  
  # FAARRF: y.pos = 1, x.pos = 2:7, 15:18 (F1/F2 for 2 lags), s.pos = 50:600
  
  if (nowcast == T) {
    faarrf = MRF(datamat, y.pos = 1, x.pos = c(2:7,14:17), S.pos = 2:ncol(datamat), 
                 oos.pos = (round(0.7*nrow(datamat))+1):nrow(datamat) , B = 50)
  } else {
    faarrf = MRF(datamat, y.pos = 1, x.pos = c(2:7,15:18), S.pos = 2:ncol(datamat), 
                 oos.pos = (round(0.7*nrow(datamat))+1):nrow(datamat) , B = 50)
  }
  
  reg.out.error = faarrf$pred - Ytest
  Rsq_out[dataset, "FAARRF"] <- r.squared.fair(faarrf$pred, Ytest, Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=reg.out.error
  
  # mrf_preds <- data.frame(arrf$pred, faarrf$pred, Ytest)
  # mrf_preds_long = mrf_preds %>% mutate(index = as.numeric(rownames(.))) %>% 
  #   pivot_longer(-index, names_to = 'model', values_to = 'pred') %>% data.frame()
  # 
  # ggplot(mrf_preds_long, aes(x = index, y = pred, color = model)) +
  #   geom_line(size = 1) +
  #   ggtitle('MRF OOS Predictions vs Actual')
  
  print('Done with FAARRF')
  print(Sys.time())
  
  ########### MACRO BLOCK #####################
  
  # AR model on the past 6 lags
  if (nowcast == T) {
    reg <- lm(y~ L_1y+L_2y +L_3y + L_4y+L_5y , data=train_data)
  } else {
    reg <- lm(y~L_0y + L_1y+L_2y +L_3y + L_4y+L_5y , data=train_data)
  }
  
  reg.out.error <- predict(reg, newdata=test_data)-Ytest
  
  Rsq_in[dataset, "AR"] <- r.squared.fair(predict.lm(reg, train_data), Ytrain, Ytrain)
  Rsq_out[dataset, "AR"] <- r.squared.fair(predict.lm(reg, test_data), Ytest, Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=reg.out.error
    
  # FAAR model - adding the Macro Factors
  # reg <- lm(y~L_0y + L_1y+L_2y +L_3y + L_4y+L_5y+
  #             L0_F_UK1+L0_F_UK2+L1_F_UK1+L1_F_UK2, data=train_data)
  
  if (nowcast == T) {
    reg <- lm(y ~ L_1y+L_2y +L_3y + L_4y+L_5y +
                L0_F1 + L0_F2 + L1_F1 + L1_F2, data = train_data)
  } else {
    reg <- lm(y ~ L_0y + L_1y+L_2y +L_3y + L_4y+L_5y +
                L0_F1 + L0_F2 + L1_F1 + L1_F2, data = train_data)
  }
  
  reg.out.error <- predict(reg, newdata=test_data)-Ytest
  Rsq_in[dataset, "FA-AR"] <- r.squared.fair(predict.lm(reg, train_data), Ytrain, Ytrain)
  Rsq_out[dataset, "FA-AR"] <- r.squared.fair(predict.lm(reg, test_data), Ytest, Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=reg.out.error
  
  print('Done with Macro Block')
  
  ########### BENCHMARKS BLOCK ####################
  
  
  # Lasso Regression
  # (Perform k-fold cross-validation to get lambda hyp for lasso, type.measure refers to loss to use for cross-validation)
  lasso.lambda <- cv.glmnet(Xtrain, Ytrain, alpha=1, nfolds=10, type.measure = "deviance")$lambda.1se
  # Run lasso
  lasso <- glmnet(Xtrain, Ytrain, lambda=lasso.lambda)
  # Get oos errors
  lasso.out.error <- predict(lasso, Xtest)-Ytest
  
  Rsq_in[dataset, "lasso"] <- r.squared.fair(predict(lasso, Xtrain), Ytrain, Ytrain)
  Rsq_out[dataset, "lasso"] <- r.squared.fair(predict(lasso, Xtest), Ytest,Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=lasso.out.error
  
  print('Done with Lasso')
  
  # Ridge regression
  ridge.lambda <- cv.glmnet(Xtrain, Ytrain, alpha=0, nfolds=10, type.measure = "deviance")$lambda.1se
  ridge <- glmnet(Xtrain, Ytrain, lambda=ridge.lambda)
  ridge.out.error <- predict(ridge, Xtest)-Ytest
  Rsq_in[dataset, "ridge"] <- r.squared.fair(predict(ridge, Xtrain), Ytrain, Ytrain)
  Rsq_out[dataset, "ridge"] <- r.squared.fair(predict(ridge, Xtest), Ytest,Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=ridge.out.error
  
  print('Done with Ridge')
  
  
  # single tree
  tuned_mars <- caret::train(
    x = Xtrain,
    y = as.vector(Ytrain),
    method = "rpart",
    metric = "RMSE", #interaction.depth = 5, #shrinkage=0.001,
    trControl = trainControl(method = "cv", number = 10) #,verbose=FALSE #,
    #tuneGrid = hyper_grid
  )
  # Tuned RPART tree model
  single=tuned_mars$finalModel
  
  #single <-rpart(form,data=train_data, minsplit=10, cp=0.01)
  single.out.error <- predict(single, as.data.frame(Xtest))-Ytest
  Rsq_in[dataset, "single tree"] <- r.squared.fair(predict(single, as.data.frame(Xtrain)), Ytrain, Ytrain)
  Rsq_out[dataset, "single tree"] <- r.squared.fair(predict(single, as.data.frame(Xtest)), Ytest,Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=single.out.error
  
  # Random Forest - untuned
  rf <- ranger(form, data=train_data, mtry=ncol(train_data)/3, importance='permutation',num.trees =500, min.node.size=3,num.threads = 1)
  rf.out.error <- predict(rf, test_data)$prediction-Ytest
  Rsq_in[dataset, "random forest"] <- r.squared.fair(predict(rf, train_data)$prediction, Ytrain, Ytrain)
  Rsq_out[dataset, "random forest"] <- r.squared.fair(predict(rf, test_data)$prediction, Ytest, Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=rf.out.error
  
  print('Done with RF Untuned')
  
  # Random Forest - after tuning
  hyper_grid <-  expand.grid(mtry = round(ncol(Xtrain)*c(0.2,0.33,0.5,1)),
                             min.node.size=c(1,3,5,10,20,40,60,100),
                             splitrule='variance') #,
  
  tuned_mars <- caret::train(
    x = Xtrain,
    y = as.vector(Ytrain),
    method = "ranger",
    metric = "RMSE", #interaction.depth = 5, 
    #shrinkage=0.001,
    trControl = trainControl(method = "cv", number = 10),verbose=FALSE, #,
    tuneGrid = hyper_grid
  )
  
  suggested.degree=as.numeric(tuned_mars$bestTune[2])
  mod = tuned_mars$finalModel
  out.error <- predict(mod, as.data.frame(Xtest))$predictions -Ytest
  Rsq_in[dataset, 'RF tuned'] <- r.squared.fair(predict(mod, as.data.frame(Xtrain))$predictions, Ytrain, Ytrain)
  Rsq_out[dataset, 'RF tuned'] <- r.squared.fair(predict(mod, as.data.frame(Xtest))$predictions, Ytest, Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=out.error
  
  print('Done with RF Tuned')
  
  # GBM Tuned
  hyper_grid <-  expand.grid(interaction.depth = c(5),
                             n.trees = c(50,100,200,300,500,750,100,1250,1500,2000,2500), 
                             shrinkage = seq(.001, .01,.1),
                             n.minobsinnode = 10)
  
  tuned_mars <- caret::train(
    x = Xtrain,
    y = as.vector(Ytrain),
    method = "gbm",
    metric = "RMSE", #interaction.depth = 5, 
    #shrinkage=0.001,
    trControl = trainControl(method = "cv", number = 10),verbose=FALSE, #,
    tuneGrid = hyper_grid
  )
  
  #?
  suggested.degree=as.numeric(tuned_mars$bestTune[2])
  mod = tuned_mars$finalModel
  out.error <- predict(mod, as.data.frame(Xtest),n.trees = mod$n.trees)-Ytest
  Rsq_out[dataset, 'GBM tuned'] <- r.squared.fair(predict(mod, as.data.frame(Xtrain),n.trees = mod$n.trees), Ytrain, Ytrain)
  Rsq_out[dataset, 'GBM tuned'] <- r.squared.fair(predict(mod, as.data.frame(Xtest),n.trees = mod$n.trees), Ytest, Ytrain)
  m=m+1
  err[m,1:length(Ytest)]=out.error
  
  print('Done with GBM Tuned')
  
  
  # Settings for the plots
  par(las=2) # make label text perpendicular to axis
  par(mar=c(5,8,4,2)) # increase y-axis margin.
  par(mfrow=c(1,1))
  
  # Save Path
  path = paste(wd, 'ML Results/')
  
  # Make a matrix of the OOS Rsquared 
  png(file=paste(path,data_file,save_name,'.png',sep=''))
  barplot(as.matrix(Rsq_out)[dataset,],horiz=TRUE,col=1:ncol(Rsq_in),las=2)
  title(paste(dataset,'R^2 out',sep=' - '))
  dev.off()
  
  # Save the data for multiple datasets together
  filename=paste(path,data_file,save_name,'.RData',sep='') #,yymmdd ,'_v',1,'
  save(Rsq_in,Rsq_out,err,dataset,Ytest,file=filename)

}



#,'macrotoy_gdph2', 'macrotoy_unemph2', 'macrotoy_infh2'

data_files <- c('macrotoy_unemph1', 'macrotoy_infh1')
save_name = '22feb_nowcast'

for (data_file in data_files) {
  
  run_models_fn(data_file, save_name = save_name, nowcast = T)
  print(paste('Done with ', data_file))
  print(Sys.time())
  
}


