#Library
library(dr)
library(glmnet)
library(SISIR)
library(plsRglm)
library(AICcmodavg)
library(ggplot2)
library(doParallel)
library(boot)
#library(glmnet)
library(DMwR)
library(caret)
###Set seed------------------------------------
set.seed(1)
setwd("C:/Users/nealf/OneDrive/My R Work and Data/TCGNewdata")

###Searching Program---------------------------------------
####AICc
fastsearch.aicc <- function(StartPoint){
  
  #Model setting
  StartModel <- StartPoint$StartModelIndex
  Multiplier <- StartPoint$Multiplier
  Seq_Length <- StartPoint$Seq
  Cut_Length <- StartPoint$Cut
  threshold <- StartPoint$threshold
  distribution <- StartPoint$distribution
  z <- StartPoint$database
  
  #sigmoid function with multiplied by n
  Sigmoid <- function(x,n){ return(1/(1+exp(-n*x)))}
  
  #conditional probablity based on AICc
  ConditionalProbabilityAICc <- function(glm1, glm2, n){ return(Sigmoid(AICc(glm2)-AICc(glm1),n))}
  
  #Model selection 
  ModelChoice1 <- function(ModelSeletionMatrix, threshold, Cut_Length){
    CutModelMatrix <- ModelSeletionMatrix[(nrow(ModelSeletionMatrix) - Cut_Length + 1):nrow(ModelSeletionMatrix),]
    return(colSums(CutModelMatrix)/Cut_Length > threshold)
  }
  
  #Gibbs sampler
  GibbsSampler.fs <- function(StartModel, Seq_Length, Multiplier, ConditionalProbability, data, distribution){
    
    #Start
    database <- data
    n <- Multiplier
    StartIndex <- StartModel
    StartModel <- c(1, StartIndex)
    ModelSeletionMatrix <- StartModel
    
    #Loop
    j <- 0
    while(j < Seq_Length){
      for(i in 1:ncol(database) - 1){
        TempIndex <- StartIndex
        TempDual <- TempIndex
        TempDual[i] <- 1 - TempIndex[i]
        
        TempCheck <- c(1,TempIndex) == 1
        TempCheckDual <- c(1, TempDual) == 1
        
        if(sum(!(TempIndex == 0)) == 0){
          y <- database[,TempCheck]
          Model_One <- glm(y ~ 1,family = distribution)
        }else{
          Model_One <- glm(y ~ .,family = distribution, data = database[,TempCheck])
        }
        
        if(sum(!(TempDual == 0)) == 0){
          y <- database[,TempCheckDual]
          Model_Two <- glm(y ~ 1,family = distribution)
        }else{
          Model_Two <- glm(y ~ .,family = distribution, data = database[,TempCheckDual])
        }
        
        #Random jump
        BinomProb <- ConditionalProbability(Model_One, Model_Two, n)
        
        TossUp <- rbinom(1,1,BinomProb)
        if(TossUp == 1){
          TempIndex[i] <- TempIndex[i]
        }
        else{
          TempIndex[i] <- 1 - TempIndex[i]
        }
        StartIndex <- TempIndex
      }
      ModelSeletionMatrix <- rbind(ModelSeletionMatrix, c(1,StartIndex))
      j = j + 1
    }
    return(ModelSeletionMatrix)
  }
  
  #Model Selection
  ModelSeletionMatrix <- GibbsSampler.fs(StartModel, Seq_Length, Multiplier, ConditionalProbabilityAICc, z, distribution)
  SelectedModel <- ModelChoice1(ModelSeletionMatrix, threshold, Cut_Length) 
  return(SelectedModel)
  
}
AICcselect <- function(SelectedModel){ return(AICc(glm(y~., data = z[,SelectedModel], family = distribution)))}

####Occurance--------------------------
occurance <- 0.2
###Parallel
ncores <- detectCores() - 1

###Random start model -----------------------------------------------------
Multiplier <- 5
Seq_length <- 300
Cut_length <- 100
Loop <- 200
threshold_num <- 0.8
distribution <- binomial

###12 hours ahead---------------------------------
TC_data_12h_AR <- read.csv("TC_data_12h_AR_Matched.csv")
TC_data_12h_AR_NonID <- TC_data_12h_AR[,- (1:7)]

###First step LASSO search------------------------------
##Set up parameter and store results
i <- 1
VarChose_AR_12h <- vector()
Result_AR_12h <- vector()
Result_AR_12h_validation <- vector()
Metric_AR_12h <- vector()
Metric_AR_12h_validation <- vector()
##Cross-validation set up
folds <- createFolds(TC_data_12h_AR_NonID$y, k = 10, list = TRUE, returnTrain = FALSE)

while(i <= 10)
{
  ###make a partition
  index <- folds[[i]]
  TC_test <- TC_data_12h_AR_NonID[index,]
  TC_train <- TC_data_12h_AR_NonID[-index,]
  ##SMOTE method to overcome the inbalanced data
  TC_train$y <- as.factor(TC_train$y)
  TC_AR_12h_New <- SMOTE(y~., TC_train)
  ##Apply LASSO to do feature selection
  X_AR_12h <- as.matrix(TC_AR_12h_New[,2:150])
  Y_AR_12h <- as.matrix(TC_AR_12h_New[,1])
  ##LASSO
  LASSO.cv_AR_12h <- cv.glmnet(X_AR_12h, Y_AR_12h, family = c("binomial"))
  ##Results of LASSO
  coefs <- drop(predict(LASSO.cv_AR_12h, type = "coef"))
  ##Results of variable selection
  VarChose_AR_12h <- rbind(VarChose_AR_12h, as.numeric(coefs[-1]!=0))
  ##Prediction
  X_AR_12h <- as.matrix(TC_test[,2:150])
  Y_AR_12h <- as.matrix(TC_test[,1])
  Y_AR_12h_Pred <- predict(LASSO.cv_AR_12h, newx = X_AR_12h, type = c("class") )
  Y_AR_12h_Pred <- as.numeric(Y_AR_12h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_12h == 1 & Y_AR_12h_Pred == 1)
  FN <- sum(Y_AR_12h == 1 & Y_AR_12h_Pred == 0)
  FP <- sum(Y_AR_12h == 0 & Y_AR_12h_Pred == 1)
  TN <- sum(Y_AR_12h == 0 & Y_AR_12h_Pred == 0)
  Result_AR_12h_test <- cbind(TP,FN,FP,TN)
  Result_AR_12h <- rbind(Result_AR_12h, Result_AR_12h_test)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_12h_test <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_12h <- rbind(Metric_AR_12h, Metric_AR_12h_test)
  ##Validation on the train set
  X_AR_12h <- as.matrix(TC_train[,2:150])
  Y_AR_12h <- as.matrix(TC_train[,1])
  Y_AR_12h_Pred <- predict(LASSO.cv_AR_12h, newx = X_AR_12h, type = c("class") )
  Y_AR_12h_Pred <- as.numeric(Y_AR_12h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_12h == 1 & Y_AR_12h_Pred == 1)
  FN <- sum(Y_AR_12h == 1 & Y_AR_12h_Pred == 0)
  FP <- sum(Y_AR_12h == 0 & Y_AR_12h_Pred == 1)
  TN <- sum(Y_AR_12h == 0 & Y_AR_12h_Pred == 0)
  Result_AR_12h_train <- cbind(TP,FN,FP,TN)
  Result_AR_12h_validation <- rbind(Result_AR_12h_validation, Result_AR_12h_train)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_12h_train <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_12h_validation <- rbind(Metric_AR_12h_validation, Metric_AR_12h_train)
  ##Looper
  #print(i)
  i <- i+1
}
Result_AR_12h
Result_AR_12h_validation
Metric_AR_12h
Metric_AR_12h_validation

###Selected data-----------------------------------
s <- colSums(VarChose_AR_12h)/nrow(VarChose_AR_12h)
#barplot(s)
VarSel_AR_12h <- c(1,s) > occurance
TC_data_12h_AR_Select <- TC_data_12h_AR_NonID[VarSel_AR_12h]

###SMOTE sample
TC_data_12h_AR_Select$y <- as.factor(TC_data_12h_AR_Select$y)
TC_AR_12h_SMOTE <- SMOTE(y~., TC_data_12h_AR_Select)
#cor(TC_data_12h_SPO_Select)

###Fast Gibbs search------------------------------------------------------------------
z <- TC_AR_12h_SMOTE
p <- dim(z)[2]-1
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel_12h_AR <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel_12h_AR <- SelectedLModel_12h_AR[,colSums(SelectedLModel_12h_AR)!=1]

LModel_12h_AR_table <- apply(SelectedLModel_12h_AR,2,AICcselect)
table(LModel_12h_AR_table)

z_12h_AR <- z[,SelectedLModel_12h_AR[,which.min(LModel_12h_AR_table)]]
LModel_12h_AR <- glm(y~., data = z_12h_AR, family = distribution)
summary(LModel_12h_AR)
stopCluster(cl)


##Prediction
z1_12h_AR <- TC_data_12h_AR_Select[,SelectedLModel_12h_AR[,which.min(LModel_12h_AR_table)]]
Y_AR_12h_Pred_Prob <- predict(LModel_12h_AR, newdata = z1_12h_AR[,-1], type = "response")
Y_AR_12h_Pred <- as.numeric(Y_AR_12h_Pred_Prob > 0.5)
Y_AR_12h <- z1_12h_AR[,1]
##confusion matrix for binary classification
TP <- sum(Y_AR_12h == 1 & Y_AR_12h_Pred == 1)
FN <- sum(Y_AR_12h == 1 & Y_AR_12h_Pred == 0)
FP <- sum(Y_AR_12h == 0 & Y_AR_12h_Pred == 1)
TN <- sum(Y_AR_12h == 0 & Y_AR_12h_Pred == 0)
Result_AR_12h_Gibbs <- cbind(TP,FN,FP,TN)
##metric for imbalanaced data
Sensitivity <- TP/(TP+FN)
Specificity <- TN/(TN+FP)
Precision <-TP/(TP+FP)
Recall <- TP/(TP+FN)
F_1 <- 2*TP/(2*TP+FP+FN)
CSI <- TP/(TP+FP+FN)
Metric_AR_12h_Gibbs <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)

Result_AR_12h_Gibbs
Metric_AR_12h_Gibbs

###24 hours ahead---------------------------------
TC_data_24h_AR <- read.csv("TC_data_24h_AR_Matched.csv")
TC_data_24h_AR_NonID <- TC_data_24h_AR[,- (1:7)]

###First step LASSO search------------------------------
##Set up parameter and store results
i <- 1
VarChose_AR_24h <- vector()
Result_AR_24h <- vector()
Result_AR_24h_validation <- vector()
Metric_AR_24h <- vector()
Metric_AR_24h_validation <- vector()
##Cross-validation set up
folds <- createFolds(TC_data_24h_AR_NonID$y, k = 10, list = TRUE, returnTrain = FALSE)

while(i <= 10)
{
  ###make a partition
  index <- folds[[i]]
  TC_test <- TC_data_24h_AR_NonID[index,]
  TC_train <- TC_data_24h_AR_NonID[-index,]
  ##SMOTE method to overcome the inbalanced data
  TC_train$y <- as.factor(TC_train$y)
  TC_AR_24h_New <- SMOTE(y~., TC_train)
  ##Apply LASSO to do feature selection
  X_AR_24h <- as.matrix(TC_AR_24h_New[,2:150])
  Y_AR_24h <- as.matrix(TC_AR_24h_New[,1])
  ##LASSO
  LASSO.cv_AR_24h <- cv.glmnet(X_AR_24h, Y_AR_24h, family = c("binomial"))
  ##Results of LASSO
  coefs <- drop(predict(LASSO.cv_AR_24h, type = "coef"))
  ##Results of variable selection
  VarChose_AR_24h <- rbind(VarChose_AR_24h, as.numeric(coefs[-1]!=0))
  ##Prediction
  X_AR_24h <- as.matrix(TC_test[,2:150])
  Y_AR_24h <- as.matrix(TC_test[,1])
  Y_AR_24h_Pred <- predict(LASSO.cv_AR_24h, newx = X_AR_24h, type = c("class") )
  Y_AR_24h_Pred <- as.numeric(Y_AR_24h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_24h == 1 & Y_AR_24h_Pred == 1)
  FN <- sum(Y_AR_24h == 1 & Y_AR_24h_Pred == 0)
  FP <- sum(Y_AR_24h == 0 & Y_AR_24h_Pred == 1)
  TN <- sum(Y_AR_24h == 0 & Y_AR_24h_Pred == 0)
  Result_AR_24h_test <- cbind(TP,FN,FP,TN)
  Result_AR_24h <- rbind(Result_AR_24h, Result_AR_24h_test)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_24h_test <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_24h <- rbind(Metric_AR_24h, Metric_AR_24h_test)
  ##Validation on the train set
  X_AR_24h <- as.matrix(TC_train[,2:150])
  Y_AR_24h <- as.matrix(TC_train[,1])
  Y_AR_24h_Pred <- predict(LASSO.cv_AR_24h, newx = X_AR_24h, type = c("class") )
  Y_AR_24h_Pred <- as.numeric(Y_AR_24h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_24h == 1 & Y_AR_24h_Pred == 1)
  FN <- sum(Y_AR_24h == 1 & Y_AR_24h_Pred == 0)
  FP <- sum(Y_AR_24h == 0 & Y_AR_24h_Pred == 1)
  TN <- sum(Y_AR_24h == 0 & Y_AR_24h_Pred == 0)
  Result_AR_24h_train <- cbind(TP,FN,FP,TN)
  Result_AR_24h_validation <- rbind(Result_AR_24h_validation, Result_AR_24h_train)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_24h_train <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_24h_validation <- rbind(Metric_AR_24h_validation, Metric_AR_24h_train)
  ##Looper
  #print(i)
  i <- i+1
}
Result_AR_24h
Result_AR_24h_validation
Metric_AR_24h
Metric_AR_24h_validation

###Selected data-----------------------------------
s <- colSums(VarChose_AR_24h)/nrow(VarChose_AR_24h)
#barplot(s)
VarSel_AR_24h <- c(1,s) > occurance
TC_data_24h_AR_Select <- TC_data_24h_AR_NonID[VarSel_AR_24h]

###SMOTE sample
TC_data_24h_AR_Select$y <- as.factor(TC_data_24h_AR_Select$y)
TC_AR_24h_SMOTE <- SMOTE(y~., TC_data_24h_AR_Select)
#cor(TC_data_24h_SPO_Select)

###Fast Gibbs search------------------------------------------------------------------
z <- TC_AR_24h_SMOTE
p <- dim(z)[2]-1
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel_24h_AR <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel_24h_AR <- SelectedLModel_24h_AR[,colSums(SelectedLModel_24h_AR)!=1]

LModel_24h_AR_table <- apply(SelectedLModel_24h_AR,2,AICcselect)
table(LModel_24h_AR_table)

z_24h_AR <- z[,SelectedLModel_24h_AR[,which.min(LModel_24h_AR_table)]]
LModel_24h_AR <- glm(y~., data = z_24h_AR, family = distribution)
summary(LModel_24h_AR)
stopCluster(cl)


##Prediction
z1_24h_AR <- TC_data_24h_AR_Select[,SelectedLModel_24h_AR[,which.min(LModel_24h_AR_table)]]
Y_AR_24h_Pred <- predict(LModel_24h_AR, newdata = z1_24h_AR[,-1], type = "response")
Y_AR_24h_Pred <- as.numeric(Y_AR_24h_Pred > 0.5)
Y_AR_24h <- z1_24h_AR[,1]
##confusion matrix for binary classification
TP <- sum(Y_AR_24h == 1 & Y_AR_24h_Pred == 1)
FN <- sum(Y_AR_24h == 1 & Y_AR_24h_Pred == 0)
FP <- sum(Y_AR_24h == 0 & Y_AR_24h_Pred == 1)
TN <- sum(Y_AR_24h == 0 & Y_AR_24h_Pred == 0)
Result_AR_24h <- cbind(TP,FN,FP,TN)
##metric for imbalanaced data
Sensitivity <- TP/(TP+FN)
Specificity <- TN/(TN+FP)
Precision <-TP/(TP+FP)
Recall <- TP/(TP+FN)
F_1 <- 2*TP/(2*TP+FP+FN)
CSI <- TP/(TP+FP+FN)
Metric_AR_24h <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)

Result_AR_24h
Metric_AR_24h


###36 hours ahead---------------------------------
TC_data_36h_AR <- read.csv("TC_data_36h_AR_Matched.csv")
TC_data_36h_AR_NonID <- TC_data_36h_AR[,- (1:7)]

###First step LASSO search------------------------------
##Set up parameter and store results
i <- 1
VarChose_AR_36h <- vector()
Result_AR_36h <- vector()
Result_AR_36h_validation <- vector()
Metric_AR_36h <- vector()
Metric_AR_36h_validation <- vector()
##Cross-validation set up
folds <- createFolds(TC_data_36h_AR_NonID$y, k = 10, list = TRUE, returnTrain = FALSE)

while(i <= 10)
{
  ###make a partition
  index <- folds[[i]]
  TC_test <- TC_data_36h_AR_NonID[index,]
  TC_train <- TC_data_36h_AR_NonID[-index,]
  ##SMOTE method to overcome the inbalanced data
  TC_train$y <- as.factor(TC_train$y)
  TC_AR_36h_New <- SMOTE(y~., TC_train)
  ##Apply LASSO to do feature selection
  X_AR_36h <- as.matrix(TC_AR_36h_New[,2:150])
  Y_AR_36h <- as.matrix(TC_AR_36h_New[,1])
  ##LASSO
  LASSO.cv_AR_36h <- cv.glmnet(X_AR_36h, Y_AR_36h, family = c("binomial"))
  ##Results of LASSO
  coefs <- drop(predict(LASSO.cv_AR_36h, type = "coef"))
  ##Results of variable selection
  VarChose_AR_36h <- rbind(VarChose_AR_36h, as.numeric(coefs[-1]!=0))
  ##Prediction
  X_AR_36h <- as.matrix(TC_test[,2:150])
  Y_AR_36h <- as.matrix(TC_test[,1])
  Y_AR_36h_Pred <- predict(LASSO.cv_AR_36h, newx = X_AR_36h, type = c("class") )
  Y_AR_36h_Pred <- as.numeric(Y_AR_36h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_36h == 1 & Y_AR_36h_Pred == 1)
  FN <- sum(Y_AR_36h == 1 & Y_AR_36h_Pred == 0)
  FP <- sum(Y_AR_36h == 0 & Y_AR_36h_Pred == 1)
  TN <- sum(Y_AR_36h == 0 & Y_AR_36h_Pred == 0)
  Result_AR_36h_test <- cbind(TP,FN,FP,TN)
  Result_AR_36h <- rbind(Result_AR_36h, Result_AR_36h_test)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_36h_test <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_36h <- rbind(Metric_AR_36h, Metric_AR_36h_test)
  ##Validation on the train set
  X_AR_36h <- as.matrix(TC_train[,2:150])
  Y_AR_36h <- as.matrix(TC_train[,1])
  Y_AR_36h_Pred <- predict(LASSO.cv_AR_36h, newx = X_AR_36h, type = c("class") )
  Y_AR_36h_Pred <- as.numeric(Y_AR_36h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_36h == 1 & Y_AR_36h_Pred == 1)
  FN <- sum(Y_AR_36h == 1 & Y_AR_36h_Pred == 0)
  FP <- sum(Y_AR_36h == 0 & Y_AR_36h_Pred == 1)
  TN <- sum(Y_AR_36h == 0 & Y_AR_36h_Pred == 0)
  Result_AR_36h_train <- cbind(TP,FN,FP,TN)
  Result_AR_36h_validation <- rbind(Result_AR_36h_validation, Result_AR_36h_train)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_36h_train <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_36h_validation <- rbind(Metric_AR_36h_validation, Metric_AR_36h_train)
  ##Looper
  #print(i)
  i <- i+1
}
Result_AR_36h
Result_AR_36h_validation
Metric_AR_36h
Metric_AR_36h_validation

###Selected data-----------------------------------
s <- colSums(VarChose_AR_36h)/nrow(VarChose_AR_36h)
#barplot(s)
VarSel_AR_36h <- c(1,s) > occurance
TC_data_36h_AR_Select <- TC_data_36h_AR_NonID[VarSel_AR_36h]

###SMOTE sample
TC_data_36h_AR_Select$y <- as.factor(TC_data_36h_AR_Select$y)
TC_AR_36h_SMOTE <- SMOTE(y~., TC_data_36h_AR_Select)
#cor(TC_data_36h_SPO_Select)

###Fast Gibbs search------------------------------------------------------------------
z <- TC_AR_36h_SMOTE
p <- dim(z)[2]-1
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel_36h_AR <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel_36h_AR <- SelectedLModel_36h_AR[,colSums(SelectedLModel_36h_AR)!=1]

LModel_36h_AR_table <- apply(SelectedLModel_36h_AR,2,AICcselect)
table(LModel_36h_AR_table)

z_36h_AR <- z[,SelectedLModel_36h_AR[,which.min(LModel_36h_AR_table)]]
LModel_36h_AR <- glm(y~., data = z_36h_AR, family = distribution)
summary(LModel_36h_AR)
stopCluster(cl)


##Prediction
z1_36h_AR <- TC_data_36h_AR_Select[,SelectedLModel_36h_AR[,which.min(LModel_36h_AR_table)]]
Y_AR_36h_Pred <- predict(LModel_36h_AR, newdata = z1_36h_AR[,-1], type = "response")
Y_AR_36h_Pred <- as.numeric(Y_AR_36h_Pred > 0.5)
Y_AR_36h <- z1_36h_AR[,1]
##confusion matrix for binary classification
TP <- sum(Y_AR_36h == 1 & Y_AR_36h_Pred == 1)
FN <- sum(Y_AR_36h == 1 & Y_AR_36h_Pred == 0)
FP <- sum(Y_AR_36h == 0 & Y_AR_36h_Pred == 1)
TN <- sum(Y_AR_36h == 0 & Y_AR_36h_Pred == 0)
Result_AR_36h <- cbind(TP,FN,FP,TN)
##metric for imbalanaced data
Sensitivity <- TP/(TP+FN)
Specificity <- TN/(TN+FP)
Precision <-TP/(TP+FP)
Recall <- TP/(TP+FN)
F_1 <- 2*TP/(2*TP+FP+FN)
CSI <- TP/(TP+FP+FN)
Metric_AR_36h <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)

Result_AR_36h
Metric_AR_36h

###48 hours ahead---------------------------------
TC_data_48h_AR <- read.csv("TC_data_48h_AR_Matched.csv")
TC_data_48h_AR_NonID <- TC_data_48h_AR[,- (1:7)]

###First step LASSO search------------------------------
##Set up parameter and store results
i <- 1
VarChose_AR_48h <- vector()
Result_AR_48h <- vector()
Result_AR_48h_validation <- vector()
Metric_AR_48h <- vector()
Metric_AR_48h_validation <- vector()
##Cross-validation set up
folds <- createFolds(TC_data_48h_AR_NonID$y, k = 10, list = TRUE, returnTrain = FALSE)

while(i <= 10)
{
  ###make a partition
  index <- folds[[i]]
  TC_test <- TC_data_48h_AR_NonID[index,]
  TC_train <- TC_data_48h_AR_NonID[-index,]
  ##SMOTE method to overcome the inbalanced data
  TC_train$y <- as.factor(TC_train$y)
  TC_AR_48h_New <- SMOTE(y~., TC_train)
  ##Apply LASSO to do feature selection
  X_AR_48h <- as.matrix(TC_AR_48h_New[,2:150])
  Y_AR_48h <- as.matrix(TC_AR_48h_New[,1])
  ##LASSO
  LASSO.cv_AR_48h <- cv.glmnet(X_AR_48h, Y_AR_48h, family = c("binomial"))
  ##Results of LASSO
  coefs <- drop(predict(LASSO.cv_AR_48h, type = "coef"))
  ##Results of variable selection
  VarChose_AR_48h <- rbind(VarChose_AR_48h, as.numeric(coefs[-1]!=0))
  ##Prediction
  X_AR_48h <- as.matrix(TC_test[,2:150])
  Y_AR_48h <- as.matrix(TC_test[,1])
  Y_AR_48h_Pred <- predict(LASSO.cv_AR_48h, newx = X_AR_48h, type = c("class") )
  Y_AR_48h_Pred <- as.numeric(Y_AR_48h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_48h == 1 & Y_AR_48h_Pred == 1)
  FN <- sum(Y_AR_48h == 1 & Y_AR_48h_Pred == 0)
  FP <- sum(Y_AR_48h == 0 & Y_AR_48h_Pred == 1)
  TN <- sum(Y_AR_48h == 0 & Y_AR_48h_Pred == 0)
  Result_AR_48h_test <- cbind(TP,FN,FP,TN)
  Result_AR_48h <- rbind(Result_AR_48h, Result_AR_48h_test)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_48h_test <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_48h <- rbind(Metric_AR_48h, Metric_AR_48h_test)
  ##Validation on the train set
  X_AR_48h <- as.matrix(TC_train[,2:150])
  Y_AR_48h <- as.matrix(TC_train[,1])
  Y_AR_48h_Pred <- predict(LASSO.cv_AR_48h, newx = X_AR_48h, type = c("class") )
  Y_AR_48h_Pred <- as.numeric(Y_AR_48h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_48h == 1 & Y_AR_48h_Pred == 1)
  FN <- sum(Y_AR_48h == 1 & Y_AR_48h_Pred == 0)
  FP <- sum(Y_AR_48h == 0 & Y_AR_48h_Pred == 1)
  TN <- sum(Y_AR_48h == 0 & Y_AR_48h_Pred == 0)
  Result_AR_48h_train <- cbind(TP,FN,FP,TN)
  Result_AR_48h_validation <- rbind(Result_AR_48h_validation, Result_AR_48h_train)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_48h_train <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_48h_validation <- rbind(Metric_AR_48h_validation, Metric_AR_48h_train)
  ##Looper
  #print(i)
  i <- i+1
}
Result_AR_48h
Result_AR_48h_validation
Metric_AR_48h
Metric_AR_48h_validation

###Selected data-----------------------------------
s <- colSums(VarChose_AR_48h)/nrow(VarChose_AR_48h)
#barplot(s)
VarSel_AR_48h <- c(1,s) > occurance
TC_data_48h_AR_Select <- TC_data_48h_AR_NonID[VarSel_AR_48h]

###SMOTE sample
TC_data_48h_AR_Select$y <- as.factor(TC_data_48h_AR_Select$y)
TC_AR_48h_SMOTE <- SMOTE(y~., TC_data_48h_AR_Select)
#cor(TC_data_48h_SPO_Select)

###Fast Gibbs search------------------------------------------------------------------
z <- TC_AR_48h_SMOTE
p <- dim(z)[2]-1
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel_48h_AR <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel_48h_AR <- SelectedLModel_48h_AR[,colSums(SelectedLModel_48h_AR)!=1]

LModel_48h_AR_table <- apply(SelectedLModel_48h_AR,2,AICcselect)
table(LModel_48h_AR_table)

z_48h_AR <- z[,SelectedLModel_48h_AR[,which.min(LModel_48h_AR_table)]]
LModel_48h_AR <- glm(y~., data = z_48h_AR, family = distribution)
summary(LModel_48h_AR)
stopCluster(cl)


##Prediction
z1_48h_AR <- TC_data_48h_AR_Select[,SelectedLModel_48h_AR[,which.min(LModel_48h_AR_table)]]
Y_AR_48h_Pred <- predict(LModel_48h_AR, newdata = z1_48h_AR[,-1], type = "response")
Y_AR_48h_Pred <- as.numeric(Y_AR_48h_Pred > 0.5)
Y_AR_48h <- z1_48h_AR[,1]
##confusion matrix for binary classification
TP <- sum(Y_AR_48h == 1 & Y_AR_48h_Pred == 1)
FN <- sum(Y_AR_48h == 1 & Y_AR_48h_Pred == 0)
FP <- sum(Y_AR_48h == 0 & Y_AR_48h_Pred == 1)
TN <- sum(Y_AR_48h == 0 & Y_AR_48h_Pred == 0)
Result_AR_48h <- cbind(TP,FN,FP,TN)
##metric for imbalanaced data
Sensitivity <- TP/(TP+FN)
Specificity <- TN/(TN+FP)
Precision <-TP/(TP+FP)
Recall <- TP/(TP+FN)
F_1 <- 2*TP/(2*TP+FP+FN)
CSI <- TP/(TP+FP+FN)
Metric_AR_48h <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)

Result_AR_48h
Metric_AR_48h

###60 hours ahead---------------------------------
TC_data_60h_AR <- read.csv("TC_data_60h_AR_Matched.csv")
TC_data_60h_AR_NonID <- TC_data_60h_AR[,- (1:7)]

###First step LASSO search------------------------------
##Set up parameter and store results
i <- 1
VarChose_AR_60h <- vector()
Result_AR_60h <- vector()
Result_AR_60h_validation <- vector()
Metric_AR_60h <- vector()
Metric_AR_60h_validation <- vector()
##Cross-validation set up
folds <- createFolds(TC_data_60h_AR_NonID$y, k = 10, list = TRUE, returnTrain = FALSE)

while(i <= 10)
{
  ###make a partition
  index <- folds[[i]]
  TC_test <- TC_data_60h_AR_NonID[index,]
  TC_train <- TC_data_60h_AR_NonID[-index,]
  ##SMOTE method to overcome the inbalanced data
  TC_train$y <- as.factor(TC_train$y)
  TC_AR_60h_New <- SMOTE(y~., TC_train)
  ##Apply LASSO to do feature selection
  X_AR_60h <- as.matrix(TC_AR_60h_New[,2:150])
  Y_AR_60h <- as.matrix(TC_AR_60h_New[,1])
  ##LASSO
  LASSO.cv_AR_60h <- cv.glmnet(X_AR_60h, Y_AR_60h, family = c("binomial"))
  ##Results of LASSO
  coefs <- drop(predict(LASSO.cv_AR_60h, type = "coef"))
  ##Results of variable selection
  VarChose_AR_60h <- rbind(VarChose_AR_60h, as.numeric(coefs[-1]!=0))
  ##Prediction
  X_AR_60h <- as.matrix(TC_test[,2:150])
  Y_AR_60h <- as.matrix(TC_test[,1])
  Y_AR_60h_Pred <- predict(LASSO.cv_AR_60h, newx = X_AR_60h, type = c("class") )
  Y_AR_60h_Pred <- as.numeric(Y_AR_60h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_60h == 1 & Y_AR_60h_Pred == 1)
  FN <- sum(Y_AR_60h == 1 & Y_AR_60h_Pred == 0)
  FP <- sum(Y_AR_60h == 0 & Y_AR_60h_Pred == 1)
  TN <- sum(Y_AR_60h == 0 & Y_AR_60h_Pred == 0)
  Result_AR_60h_test <- cbind(TP,FN,FP,TN)
  Result_AR_60h <- rbind(Result_AR_60h, Result_AR_60h_test)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_60h_test <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_60h <- rbind(Metric_AR_60h, Metric_AR_60h_test)
  ##Validation on the train set
  X_AR_60h <- as.matrix(TC_train[,2:150])
  Y_AR_60h <- as.matrix(TC_train[,1])
  Y_AR_60h_Pred <- predict(LASSO.cv_AR_60h, newx = X_AR_60h, type = c("class") )
  Y_AR_60h_Pred <- as.numeric(Y_AR_60h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_60h == 1 & Y_AR_60h_Pred == 1)
  FN <- sum(Y_AR_60h == 1 & Y_AR_60h_Pred == 0)
  FP <- sum(Y_AR_60h == 0 & Y_AR_60h_Pred == 1)
  TN <- sum(Y_AR_60h == 0 & Y_AR_60h_Pred == 0)
  Result_AR_60h_train <- cbind(TP,FN,FP,TN)
  Result_AR_60h_validation <- rbind(Result_AR_60h_validation, Result_AR_60h_train)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_60h_train <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_60h_validation <- rbind(Metric_AR_60h_validation, Metric_AR_60h_train)
  ##Looper
  #print(i)
  i <- i+1
}
Result_AR_60h
Result_AR_60h_validation
Metric_AR_60h
Metric_AR_60h_validation

###Selected data-----------------------------------
s <- colSums(VarChose_AR_60h)/nrow(VarChose_AR_60h)
#barplot(s)
VarSel_AR_60h <- c(1,s) > occurance
TC_data_60h_AR_Select <- TC_data_60h_AR_NonID[VarSel_AR_60h]

###SMOTE sample
TC_data_60h_AR_Select$y <- as.factor(TC_data_60h_AR_Select$y)
TC_AR_60h_SMOTE <- SMOTE(y~., TC_data_60h_AR_Select)
#cor(TC_data_60h_SPO_Select)

###Fast Gibbs search------------------------------------------------------------------
z <- TC_AR_60h_SMOTE
p <- dim(z)[2]-1
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel_60h_AR <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel_60h_AR <- SelectedLModel_60h_AR[,colSums(SelectedLModel_60h_AR)!=1]

LModel_60h_AR_table <- apply(SelectedLModel_60h_AR,2,AICcselect)
table(LModel_60h_AR_table)

z_60h_AR <- z[,SelectedLModel_60h_AR[,which.min(LModel_60h_AR_table)]]
LModel_60h_AR <- glm(y~., data = z_60h_AR, family = distribution)
summary(LModel_60h_AR)
stopCluster(cl)


##Prediction
z1_60h_AR <- TC_data_60h_AR_Select[,SelectedLModel_60h_AR[,which.min(LModel_60h_AR_table)]]
Y_AR_60h_Pred <- predict(LModel_60h_AR, newdata = z1_60h_AR[,-1], type = "response")
Y_AR_60h_Pred <- as.numeric(Y_AR_60h_Pred > 0.5)
Y_AR_60h <- z1_60h_AR[,1]
##confusion matrix for binary classification
TP <- sum(Y_AR_60h == 1 & Y_AR_60h_Pred == 1)
FN <- sum(Y_AR_60h == 1 & Y_AR_60h_Pred == 0)
FP <- sum(Y_AR_60h == 0 & Y_AR_60h_Pred == 1)
TN <- sum(Y_AR_60h == 0 & Y_AR_60h_Pred == 0)
Result_AR_60h <- cbind(TP,FN,FP,TN)
##metric for imbalanaced data
Sensitivity <- TP/(TP+FN)
Specificity <- TN/(TN+FP)
Precision <-TP/(TP+FP)
Recall <- TP/(TP+FN)
F_1 <- 2*TP/(2*TP+FP+FN)
CSI <- TP/(TP+FP+FN)
Metric_AR_60h <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)

Result_AR_60h
Metric_AR_60h

###72 hours ahead---------------------------------
TC_data_72h_AR <- read.csv("TC_data_72h_AR_Matched.csv")
TC_data_72h_AR_NonID <- TC_data_72h_AR[,- (1:7)]

###First step LASSO search------------------------------
##Set up parameter and store results
i <- 1
VarChose_AR_72h <- vector()
Result_AR_72h <- vector()
Result_AR_72h_validation <- vector()
Metric_AR_72h <- vector()
Metric_AR_72h_validation <- vector()
##Cross-validation set up
folds <- createFolds(TC_data_72h_AR_NonID$y, k = 10, list = TRUE, returnTrain = FALSE)

while(i <= 10)
{
  ###make a partition
  index <- folds[[i]]
  TC_test <- TC_data_72h_AR_NonID[index,]
  TC_train <- TC_data_72h_AR_NonID[-index,]
  ##SMOTE method to overcome the inbalanced data
  TC_train$y <- as.factor(TC_train$y)
  TC_AR_72h_New <- SMOTE(y~., TC_train)
  ##Apply LASSO to do feature selection
  X_AR_72h <- as.matrix(TC_AR_72h_New[,2:150])
  Y_AR_72h <- as.matrix(TC_AR_72h_New[,1])
  ##LASSO
  LASSO.cv_AR_72h <- cv.glmnet(X_AR_72h, Y_AR_72h, family = c("binomial"))
  ##Results of LASSO
  coefs <- drop(predict(LASSO.cv_AR_72h, type = "coef"))
  ##Results of variable selection
  VarChose_AR_72h <- rbind(VarChose_AR_72h, as.numeric(coefs[-1]!=0))
  ##Prediction
  X_AR_72h <- as.matrix(TC_test[,2:150])
  Y_AR_72h <- as.matrix(TC_test[,1])
  Y_AR_72h_Pred <- predict(LASSO.cv_AR_72h, newx = X_AR_72h, type = c("class") )
  Y_AR_72h_Pred <- as.numeric(Y_AR_72h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_72h == 1 & Y_AR_72h_Pred == 1)
  FN <- sum(Y_AR_72h == 1 & Y_AR_72h_Pred == 0)
  FP <- sum(Y_AR_72h == 0 & Y_AR_72h_Pred == 1)
  TN <- sum(Y_AR_72h == 0 & Y_AR_72h_Pred == 0)
  Result_AR_72h_test <- cbind(TP,FN,FP,TN)
  Result_AR_72h <- rbind(Result_AR_72h, Result_AR_72h_test)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_72h_test <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_72h <- rbind(Metric_AR_72h, Metric_AR_72h_test)
  ##Validation on the train set
  X_AR_72h <- as.matrix(TC_train[,2:150])
  Y_AR_72h <- as.matrix(TC_train[,1])
  Y_AR_72h_Pred <- predict(LASSO.cv_AR_72h, newx = X_AR_72h, type = c("class") )
  Y_AR_72h_Pred <- as.numeric(Y_AR_72h_Pred)
  ##confusion matrix for binary classification
  TP <- sum(Y_AR_72h == 1 & Y_AR_72h_Pred == 1)
  FN <- sum(Y_AR_72h == 1 & Y_AR_72h_Pred == 0)
  FP <- sum(Y_AR_72h == 0 & Y_AR_72h_Pred == 1)
  TN <- sum(Y_AR_72h == 0 & Y_AR_72h_Pred == 0)
  Result_AR_72h_train <- cbind(TP,FN,FP,TN)
  Result_AR_72h_validation <- rbind(Result_AR_72h_validation, Result_AR_72h_train)
  ##metric for imbalanaced data
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(TN+FP)
  Precision <-TP/(TP+FP)
  Recall <- TP/(TP+FN)
  F_1 <- 2*TP/(2*TP+FP+FN)
  CSI <- TP/(TP+FP+FN)
  Metric_AR_72h_train <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)
  Metric_AR_72h_validation <- rbind(Metric_AR_72h_validation, Metric_AR_72h_train)
  ##Looper
  #print(i)
  i <- i+1
}
Result_AR_72h
Result_AR_72h_validation
Metric_AR_72h
Metric_AR_72h_validation

###Selected data-----------------------------------
s <- colSums(VarChose_AR_72h)/nrow(VarChose_AR_72h)
#barplot(s)
VarSel_AR_72h <- c(1,s) > occurance
TC_data_72h_AR_Select <- TC_data_72h_AR_NonID[VarSel_AR_72h]

###SMOTE sample
TC_data_72h_AR_Select$y <- as.factor(TC_data_72h_AR_Select$y)
TC_AR_72h_SMOTE <- SMOTE(y~., TC_data_72h_AR_Select)
#cor(TC_data_72h_SPO_Select)

###Fast Gibbs search------------------------------------------------------------------
z <- TC_AR_72h_SMOTE
p <- dim(z)[2]-1
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel_72h_AR <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel_72h_AR <- SelectedLModel_72h_AR[,colSums(SelectedLModel_72h_AR)!=1]

LModel_72h_AR_table <- apply(SelectedLModel_72h_AR,2,AICcselect)
table(LModel_72h_AR_table)

z_72h_AR <- z[,SelectedLModel_72h_AR[,which.min(LModel_72h_AR_table)]]
LModel_72h_AR <- glm(y~., data = z_72h_AR, family = distribution)
summary(LModel_72h_AR)
stopCluster(cl)


##Prediction
z1_72h_AR <- TC_data_72h_AR_Select[,SelectedLModel_72h_AR[,which.min(LModel_72h_AR_table)]]
Y_AR_72h_Pred <- predict(LModel_72h_AR, newdata = z1_72h_AR[,-1], type = "response")
Y_AR_72h_Pred <- as.numeric(Y_AR_72h_Pred > 0.5)
Y_AR_72h <- z1_72h_AR[,1]
##confusion matrix for binary classification
TP <- sum(Y_AR_72h == 1 & Y_AR_72h_Pred == 1)
FN <- sum(Y_AR_72h == 1 & Y_AR_72h_Pred == 0)
FP <- sum(Y_AR_72h == 0 & Y_AR_72h_Pred == 1)
TN <- sum(Y_AR_72h == 0 & Y_AR_72h_Pred == 0)
Result_AR_72h <- cbind(TP,FN,FP,TN)
##metric for imbalanaced data
Sensitivity <- TP/(TP+FN)
Specificity <- TN/(TN+FP)
Precision <-TP/(TP+FP)
Recall <- TP/(TP+FN)
F_1 <- 2*TP/(2*TP+FP+FN)
CSI <- TP/(TP+FP+FN)
Metric_AR_72h <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)

Result_AR_72h
Metric_AR_72h

###LOOCV -----------------------------------------------

cv.glm2 <- function(data, glmfit, K=n){
  n <- nrow(data)
  #    glm.f <- formula(glmfit)
  glm.y <- glmfit$y
  #cost.0 <- cost(glm.y, fitted(glmfit))
  seq_len <- 1:n
  #CV <- 0
  Hindcast <- vector()
  Call <- glmfit$call
  for(i in 1:n) {
    j.out <- seq_len == i
    j.in <- seq_len != i
    ## we want data from here but formula from the parent.
    Call$data <- SMOTE(y~., data[j.in, , drop=FALSE])
    d.glm <- eval.parent(Call)
    #p.alpha <- 1/n
    Hindcast[i] <- predict(d.glm, data[j.out, , drop=FALSE], type = "response")
    #cost.i <- cost(glm.y[j.out], Hindcast[i])
    #CV <- CV + p.alpha * cost.i
    #cost.0 <- cost.0 - p.alpha *
      #cost(glm.y, predict(d.glm, data, type = "response"))
  }
  list(K = K,
       #delta = as.numeric(c(CV, CV + cost.0)),  # drop any names
       Hindcast = Hindcast)
}

##12h ahead-----------------------------------
LOOCV_AR_12h <- cv.glm2(z1_12h_AR, LModel_12h_AR)
Hindcast_Y_AR_12h <- as.numeric(LOOCV_AR_12h$Hindcast > 0.5)
##confusion matrix for binary classification
TP <- sum(Y_AR_12h == 1 & Hindcast_Y_AR_12h == 1)
FN <- sum(Y_AR_12h == 1 & Hindcast_Y_AR_12h == 0)
FP <- sum(Y_AR_12h == 0 & Hindcast_Y_AR_12h == 1)
TN <- sum(Y_AR_12h == 0 & Hindcast_Y_AR_12h == 0)
Result_AR_12h_Hindcast <- cbind(TP,FN,FP,TN)
##metric for imbalanaced data
Sensitivity <- TP/(TP+FN)
Specificity <- TN/(TN+FP)
Precision <-TP/(TP+FP)
Recall <- TP/(TP+FN)
F_1 <- 2*TP/(2*TP+FP+FN)
CSI <- TP/(TP+FP+FN)
Metric_AR_12h_Hindcast <- cbind(Sensitivity,Specificity,Precision,Recall,F_1,CSI)

Result_AR_12h_Hindcast
Metric_AR_12h_Hindcast

par(pty="s")
roc(Y_AR_12h, Y_AR_12h_Pred_Prob, plot = TRUE, legacy.axes = TRUE )
roc(Y_AR_12h, LOOCV_AR_12h$Hindcast, plot = TRUE, legacy.axes = TRUE )

###SVM
SVM_AR_12h <- svm(y~., data = z1_12h_AR )
summary(SVM_AR_12h)

x1 = attributes(predict(SVM_AR_12h, z1_12h_AR, 
                                      decision.values = TRUE))$decision.values
x2 <- predict(SVM_AR_12h, z1_12h_AR, 
              decision.values = TRUE)
roc(Y_AR_12h, as.numeric(x2), plot = TRUE, legacy.axes = TRUE )
roc(Y_AR_12h, x1, plot = TRUE, legacy.axes = TRUE )

tune_out = tune(svm, y~., data = z1_12h_AR , kernel = "radial",
                ranges = list(cost = c(0.1,1,10,100,1000), gamma = c(0.5,1,2,3,4)))
bestmod = tune_out$best.model
summary(bestmod)
