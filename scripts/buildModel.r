######################################################################################
#                   Libraries
######################################################################################

# function to install missing packages
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

# load packages
set.seed(1)
usePackage("e1071")
usePackage("randomForest")
usePackage("ROCR") 
usePackage("party")
usePackage("mlbench")
usePackage("caret")
usePackage("class")
usePackage("pROC")
usePackage("ROCR")

######################################################################################
#                   Features and load
######################################################################################
# seleciting attributes
features <- c("Class", "ORF_Length", "ORF_proportion", 'GCGG', 'TTTT', 'TCG', 'AAAA', 'ACG', 'TTGT', 'TAT', 'TAC',
              'GTT', 'GTG', 'AGT', 'CCGA', 'TACC', 'CGTG', 'CGCT', 'TACG',
              'TTAG', 'CGTA', 'ACCG', 'CCGT', 'CGGT', 'CGAC', 'CGCA', 'GCGT',
              'GTAG', 'CGTT', 'CGAA', 'GCGA', 'CGAT', 'TAGT')

# Data loading and pre-processing
train <- read.csv('temp/train.csv', sep=";", header = TRUE)
test  <- read.csv('temp/test.csv',  sep=";", header = TRUE)

# seleciting attributes
train_filtered <- train[features]
test_filtered  <- test[features]

row.has.na <- apply(train_filtered, 1, function(x){
  any(is.na(x))
  })
train_filtered <- train_filtered[!row.has.na,]

row.has.na <- apply(test_filtered, 1, function(x){
  any(is.na(x))
})
test_filtered <- test_filtered[!row.has.na,]

accuracy <- 0

print("Training and Test set loaded with success")

######################################################################################
#                   SVM
######################################################################################

print("Building svm model")

# tune svm parameters
svm_tune <- tune(svm, Class ~ ., data = train_filtered, kernel = "radial", ranges = list(cost=10^(-1:2), gamma=c(.5,1,2)),
                 tunecontrol = tune.control(sampling = "cross"))

svm_cost <- svm_tune$best.parameters$cost
svm_gamma <- svm_tune$best.parameters$gamma

print("SVM tuned.")
print("Cost:")
print(svm_cost)
print("Gamma:")
print(svm_gamma)


# building SVM model for train.csv
svm_model <- svm(Class ~ ., data = train_filtered, type = "C-classification", 
                 cost = svm_cost, gamma = svm_gamma)
# summary(svm_model) 

print("SVM model built.")

# predicts test set
prediction_svm <- predict(svm_model, newdata = test_filtered)
confusion_svm  <- table(pred = prediction_svm, true = test_filtered$Class)
confusion_svm2  <- confusionMatrix(prediction_svm, test_filtered$Class )

# evaluation metrics
# confusion matrix
accuracy_svm = sum(prediction_svm == test_filtered$Class)/length(test_filtered$Class)
precision = confusion_svm[1,1]/sum(confusion_svm[,1])
recall = confusion_svm[1,1]/sum(confusion_svm[1,])
f = 2 * (precision * recall) / (precision + recall)
cat(paste("Accuracy:\t" , format(accuracy_svm, digits=2) , "\n",sep=" "))
cat(paste("Precision:\t", format(precision, digits=2), "\n",sep=" "))
cat(paste("Recall:\t\t" , format(recall, digits=2)   , "\n",sep=" "))
cat(paste("F-measure:\t", format(f, digits=2)        , "\n",sep=" "))

if(accuracy_svm > accuracy){
  accuracy <- accuracy_svm
}


######################################################################################
#                   RF
######################################################################################
print("Building RF model")

# tunning RF
# Create model with default paramters
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(ncol(train_filtered))
tunegrid <- expand.grid(.mtry=mtry)
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
#control <- trainControl(method="repeatedcv", number=5, repeats=2, search="grid")

set.seed(seed)
#tunegrid <- expand.grid(.mtry=c(1:2))
tunegrid <- expand.grid(.mtry=c(1:15))

rf_gridsearch <- train(Class ~ ., data = train_filtered, method="rf",
                       metric = metric, tuneGrid=tunegrid, trControl=control)

mtry_tune <- rf_gridsearch$bestTune$mtry

print("RF model tuned.")
print("mtry:")
print(mtry_tune)
print("RF model built.")

png(filename="result/rf_grid.png")
plot(rf_gridsearch)
dev.off()

# rf model
rf_model <- randomForest(Class ~ ., data = train_filtered, 
                         mtry = mtry_tune, importance=TRUE, ntree = 1500)
# print(rf_model) 

print("Feature importance:\n")
# print(importance(rf_model, type = 2)) 
png(filename="result/rf_features.png")
varImpPlot(rf_model,type=2)
dev.off()

# predicts test set
prediction_rf <- predict(rf_gridsearch, newdata = test_filtered)
confusion_rf  <- table(pred = prediction_rf, true = test_filtered$Class)
confusion_rf2  <- confusionMatrix(prediction_rf, test_filtered$Class )

# evaluation metrics
# confusion matrix
print(confusion_rf2)
accuracy_rf = sum(prediction_rf == test_filtered$Class)/length(test_filtered$Class)
precision = confusion_rf[1,1]/sum(confusion_rf[,1])
recall = confusion_rf[1,1]/sum(confusion_rf[1,])
f = 2 * (precision * recall) / (precision + recall)
cat(paste("Accuracy:\t" , format(accuracy_rf, digits=2) , "\n",sep=" "))
cat(paste("Precision:\t", format(precision, digits=2), "\n",sep=" "))
cat(paste("Recall:\t\t" , format(recall, digits=2)   , "\n",sep=" "))
cat(paste("F-measure:\t", format(f, digits=2)        , "\n",sep=" "))


if(accuracy_rf > accuracy){
  accuracy <- accuracy_rf
}


######################################################################################
#                  C-tree
######################################################################################
print("Building Ctree model")

# tunning ctree
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(ncol(train_filtered))
tunegrid <- expand.grid(.mtry=mtry)
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
set.seed(seed)
tunegrid <- expand.grid(.mincriterion = .95, 
                        .maxdepth = as.integer(seq(5, 10, 2)))
ct_gridsearch <- train(Class ~ ., data = train_filtered, method="ctree2",
                       metric = metric, tuneGrid=tunegrid, trControl=control,
                       controls = ctree_control(minbucket = 10))
# print(ct_gridsearch)

png(filename="result/ctree_tuning.png")
plot(ct_gridsearch)
dev.off()


maxdepth <- ct_gridsearch$bestTune$maxdepth
mincriterion <- ct_gridsearch$bestTune$mincriterion

print("Ctree tuned")
print("Max Depth:")
print(maxdepth)
print("Min Criterion")
print(mincriterion)


print("Ctree built")


# ctree model
ct_model <- ctree(Class ~ ., data = train_filtered,
                  controls = ctree_control(maxdepth = maxdepth, 
                                           mincriterion = mincriterion))
png(filename="result/ctree_features.png")
plot(ct_model, type="simple")
dev.off()


# predicts test set
prediction_ct <- predict(ct_gridsearch, newdata = test_filtered)
confusion_ct  <- table(pred = prediction_ct, true = test_filtered$Class)
confusion_ct2  <- confusionMatrix(prediction_ct, test_filtered$Class )

# evaluation metrics
# confusion matrix
print(confusion_ct2)
accuracy_ct = sum(prediction_ct == test_filtered$Class)/length(test_filtered$Class)
precision = confusion_ct[1,1]/sum(confusion_ct[,1])
recall = confusion_ct[1,1]/sum(confusion_ct[1,])
f = 2 * (precision * recall) / (precision + recall)
cat(paste("Accuracy:\t" , format(accuracy_ct, digits=2) , "\n",sep=" "))
cat(paste("Precision:\t", format(precision, digits=2), "\n",sep=" "))
cat(paste("Recall:\t\t" , format(recall, digits=2)   , "\n",sep=" "))
cat(paste("F-measure:\t", format(f, digits=2)        , "\n",sep=" "))


if(accuracy_ct > accuracy){
  accuracy <- accuracy_ct
}


######################################################################################
#                  KNN
######################################################################################
print("Building KNN model")

# tunning knn
seed <- 7
metric <- "Accuracy"
set.seed(seed)
mtry <- sqrt(ncol(train_filtered))
tunegrid <- expand.grid(.mtry=mtry)
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="grid")
set.seed(seed)
knn_gridsearch <- train(Class ~ ., data = train_filtered, method="knn",
                       metric = metric, trControl=control,
                       preProcess = c("center","scale"), tuneLength = 20)
print(knn_gridsearch)
png(filename="result/knn_tuning.png")
plot(knn_gridsearch)
dev.off()

k <- knn_gridsearch$bestTune$k

print("Knn model tuned")
print("K:")
print(k)

# predicts test set
prediction_knn <- predict(knn_gridsearch, newdata = test_filtered)
confusion_knn  <- table(pred = prediction_knn, true = test_filtered$Class)
confusion_knn2  <- confusionMatrix(prediction_knn, test_filtered$Class )

# confusion matrix
print(confusion_knn2)

accuracy_knn = sum(prediction_knn == test_filtered$Class)/length(test_filtered$Class)
precision = confusion_knn[1,1]/sum(confusion_knn[,1])
recall = confusion_knn[1,1]/sum(confusion_knn[1,])
f = 2 * (precision * recall) / (precision + recall)
cat(paste("Accuracy:\t" , format(accuracy_knn, digits=2) , "\n",sep=" "))
cat(paste("Precision:\t", format(precision, digits=2), "\n",sep=" "))
cat(paste("Recall:\t\t" , format(recall, digits=2)   , "\n",sep=" "))
cat(paste("F-measure:\t", format(f, digits=2)        , "\n",sep=" "))

if(accuracy_knn > accuracy){
  accuracy <- accuracy_knn
}


######################################################################################
#                  INTERSECTION
######################################################################################
print("Building consensus ensemble model")

tam             <- length(prediction_svm)
total_counter   <- 0
prediction_all  <- 1:tam

# checks for prediction intersections
for (index in 1:tam){
  if(prediction_svm[index] == prediction_rf[index] &&
     prediction_svm[index] == prediction_ct[index] &&
     prediction_svm[index] == prediction_knn[index]){
    prediction_all[index] <- prediction_svm[index]
    total_counter <- total_counter + 1
  }else{
    prediction_all[index] <- NA
  }
}


# evaluation metrics
confusion_all  <- table(pred = prediction_all, true = test_filtered$Class)
print(confusion_all)
accuracy = (confusion_all[1,1] + confusion_all[2,2])/total_counter
precision = confusion_all[1,1]/sum(confusion_all[,1])
recall = confusion_all[1,1]/sum(confusion_all[1,])
f = 2 * (precision * recall) / (precision + recall)
cat(paste("Accuracy:\t" , format(accuracy, digits=2) , "\n",sep=" "))
cat(paste("Precision:\t", format(precision, digits=2), "\n",sep=" "))
cat(paste("Recall:\t\t" , format(recall, digits=2)   , "\n",sep=" "))
cat(paste("F-measure:\t", format(f, digits=2)        , "\n",sep=" "))



######################################################################################
#                  MAJORITY
######################################################################################
print("Building Majority ensemble model")

tam             <- length(prediction_svm)
total_counter   <- 0
prediction_all2  <- 1:tam
linc <- NA
pct <- NA

# class labels
for (index in 1:tam){
  voting_linc <- 0
  voting_pct <- 0
  
  if(prediction_svm[index] == 'linc'){
    linc <- prediction_svm[index]
  }else if(prediction_knn[index] == 'pct'){
    pct <- prediction_svm[index]
  }
  
  if(!is.na(linc) && !is.na(pct)){
    print(pct)
    print(linc)
    
    break;
  }
}

# checks for prediction intersections
for (index in 1:tam){
  voting_linc <- 0
  voting_pct <- 0
  
  if(prediction_svm[index] == 'linc'){
    voting_linc <- voting_linc + 1
    
    if( accuracy ==  accuracy_svm){
      voting_linc <- voting_linc + 0.1
    }
  }
  if(prediction_rf[index] == 'linc'){
    voting_linc <- voting_linc + 1
    
    if( accuracy ==  accuracy_rf){
      voting_linc <- voting_linc + 0.1
    }
  }
  if(prediction_ct[index] == 'linc'){
    voting_linc <- voting_linc + 1
    
    if( accuracy ==  accuracy_ct){
      voting_linc <- voting_linc + 0.1
    }
  }
  if(prediction_knn[index] == 'linc'){
    voting_linc <- voting_linc + 1
    
    if( accuracy ==  accuracy_knn){
      voting_linc <- voting_linc + 0.1
    }
  }
  
  if(prediction_svm[index] == 'pct'){
    pct <- prediction_svm[index]
    voting_pct <- voting_pct + 1
    
    if( accuracy ==  accuracy_svm){
      voting_pct <- voting_pct + 0.1
    }
  }
  if(prediction_rf[index] == 'pct'){
    voting_pct <- voting_pct + 1
    
    if( accuracy ==  accuracy_rf){
      voting_pct <- voting_pct + 0.1
    }
  }
  if(prediction_ct[index] == 'pct'){
    voting_pct <- voting_pct + 1
    
    if( accuracy ==  accuracy_ct){
      voting_pct <- voting_pct + 0.1
    }
  }
  if(prediction_knn[index] == 'pct'){
    voting_pct <- voting_pct + 1
    
    if( accuracy ==  accuracy_knn){
      voting_pct <- voting_pct + 0.1
    }
  }
  
  if(voting_pct > 2){
    prediction_all2[index] <- pct
    total_counter <- total_counter + 1
  }else if(voting_linc > 2){
    prediction_all2[index] <- linc
    total_counter <- total_counter + 1
  }else{
    prediction_all2[index] <- NA
  }
}

# evaluation metrics
confusion_all2  <- table(pred = prediction_all2, true = test_filtered$Class)
print(confusion_all2)
accuracy = (confusion_all2[1,1] + confusion_all2[2,2])/total_counter
precision = confusion_all2[1,1]/sum(confusion_all2[,1])
recall = confusion_all2[1,1]/sum(confusion_all2[1,])
f = 2 * (precision * recall) / (precision + recall)
cat(paste("Accuracy:\t" , format(accuracy, digits=2) , "\n",sep=" "))
cat(paste("Precision:\t", format(precision, digits=2), "\n",sep=" "))
cat(paste("Recall:\t\t" , format(recall, digits=2)   , "\n",sep=" "))
cat(paste("F-measure:\t", format(f, digits=2)        , "\n",sep=" "))


######################################################################################
#               ROC AUC
######################################################################################
pred_svm        <- prediction(as.numeric(prediction_svm), as.numeric(test_filtered$Class))
performance_svm <- performance(pred_svm, "tpr","fpr")
auc_svm         <- performance(pred_svm, "auc")@y.values[[1]]
pred_rf         <- prediction(as.numeric(prediction_rf), as.numeric(test_filtered$Class))
performance_rf  <- performance(pred_rf, "tpr","fpr")
auc_rf          <- performance(pred_rf, "auc")@y.values[[1]]
pred_ct         <- prediction(as.numeric(prediction_ct), as.numeric(test_filtered$Class))
performance_ct  <- performance(pred_ct, "tpr","fpr")
auc_ct          <- performance(pred_ct, "auc")@y.values[[1]]
pred_knn        <- prediction(as.numeric(prediction_knn), as.numeric(test_filtered$Class))
performance_knn <- performance(pred_knn, "tpr","fpr")
auc_knn         <- performance(pred_knn, "auc")@y.values[[1]]
pred_all        <- prediction(as.numeric(prediction_all), as.numeric(test_filtered$Class))
performance_all <- performance(pred_all, "tpr","fpr")
auc_all         <- performance(pred_all, "auc")@y.values[[1]]
pred_all2        <- prediction(as.numeric(prediction_all2), as.numeric(test_filtered$Class))
performance_all2 <- performance(pred_all2, "tpr","fpr")
auc_all2         <- performance(pred_all2, "auc")@y.values[[1]]

png(filename="result/roc.png")
plot(performance_svm, col='blue', main="ROC curve comparing classification performance")
legend(0.6, 0.6, c('svm', 'rf', 'ctree','knn', 'majority', 'consensus'), 2:7, fill = c('blue','red','gray','green','lightblue','purple'))
plot(performance_rf, col='red', add=TRUE,borders = "black")
plot(performance_ct, col='gray', add=TRUE,borders = "black")
plot(performance_knn, col='green', add=TRUE,borders = "black")
plot(performance_all2, col='lightblue', add=TRUE,borders = "black")
plot(performance_all, col='purple', add=TRUE,borders = "black")
dev.off()

print("AUC model:")
print(auc_all)


######################################################################################
#               Print predicted lincs
######################################################################################
write.row <- function(x, file = filepath, append = TRUE, quote = TRUE, sep=",", 
                      eol = "\n", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE, 
                      qmethod = c("escape", "double"), fileEncoding = ""){
  write.table(as.matrix(t(x)), file, append , quote, sep, 
              eol, na, dec, row.names, col.names, 
              qmethod, fileEncoding)
}

# print predicted lincs
tam             <- length(test$Class)
total_counter   <- 0

write.row("id", file="result/consensus_lincs.csv", append=FALSE, sep=";")

# checks for prediction intersections
for (index in 1:tam){
  if(!is.na(prediction_all[index])){
    # type = if(prediction_all[index] == 1) 'linc' else 'pct'
    if(prediction_all[index] == 1){
      if(!is.null(test$id[index])){
        write.row(toString(test$id[index]), file="result/consensus_lincs.csv", sep=";")
      }
    }
  }
}


# print predicted lincs
tam             <- length(test$Class)
total_counter   <- 0

write.row("id", file="result/majority_lincs.csv", append=FALSE, sep=";")

# checks for prediction intersections
for (index in 1:tam){
  if(!is.na(prediction_all2[index])){
    # type = if(prediction_all[index] == 1) 'linc' else 'pct'
    if(prediction_all2[index] == 1){
      if(!is.null(test$id[index])){
        write.row(toString(test$id[index]), file="result/majority_lincs.csv", sep=";")
      }
    }
  }
}


