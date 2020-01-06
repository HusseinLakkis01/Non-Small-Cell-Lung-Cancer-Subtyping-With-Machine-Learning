# load the required libraries
library(plyr)
library(dyplr)
library(Biobase)
library(BootMRMR)
library(mRMRe)
library(genefilter)
library(globaltest)
library(GSEABase)
library(annotate)
library(org.Hs.eg.db)
library(praznik)
library(factoextra)
library(FactoMineR)
library(caret)

# set the working directory to the one where the data is stored
setwd("/Users/husseinalilakkis/Desktop/CAPSTONE_DATA/")

#*****************************************************************************************************************
# Data Reading and Setting
#*****************************************************************************************************************
# read the meta data foe both datasets (sample description)
AD_pheno <- read.csv(file = "Adeno_metadata.csv", header = TRUE)
SC_pheno <- read.csv(file = "Squamous_metadata.csv", header = TRUE)

# set a dataframe from the csv files
AD_pheno = as.data.frame(AD_pheno)
SC_pheno = as.data.frame(SC_pheno)

#check some data
head(SC_pheno)

# adding the type of cancer to the 2 metadata frames, if sample is 2 then it is normal and append this coloumn 
# else it is cancer
AD_pheno$Status = ifelse(AD_pheno$Sam_Tissue == 2, "NORMAL", "ADLC")
SC_pheno$Status = ifelse(SC_pheno$Sam_Tissue == 2, "NORMAL", "SCLC")

# check the count to make sure it is correct using plyr method "count"
count(AD_pheno, "Status")
count(SC_pheno, "Status")

# read the final dataset that includes the expression dat and the label
#(shortcut step not to repeat the previous steps)
data = read.csv(file = "finalData.csv", header = TRUE)
data = as.data.frame(data)

# set the class vector or target vector
class <- data$STATUS

# read the gene expression data from the expression text files
AD_EXPR <- read.table(file = "Adeno_expression.txt", header = TRUE)
SC_EXPR <- read.table(file = "Squamous_expression.txt", header = TRUE)

# turn the read data into matrices 
AD_EXPR = as.matrix(AD_EXPR)
SC_EXPR = as.matrix(SC_EXPR)

# transpose the matrices
AD_EXPR = t(AD_EXPR)
SC_EXPR = t(SC_EXPR)

# instead of having numbers AS row names, use the names of samples
rownames(AD_EXPR) = AD_pheno$Sam_Patient
rownames(SC_EXPR) = SC_pheno$Sam_Patient

# transform the matrices into 2 data frames
AD_EXPR = as.data.frame(AD_EXPR)
SC_EXPR = as.data.frame(SC_EXPR)

# append the status column to each dataframe...
AD_EXPR$STATUS = AD_pheno$Status
SC_EXPR$STATUS = SC_pheno$Status

# transpore the dataset so that we have the same format as regular expression sets
# that is probes as rows and samples as columns
AD_EXPR = t(AD_EXPR)
SC_EXPR= t(SC_EXPR)

# Merge the 2 dataframes by rows to get the complete dataset
result = rbind(AD_EXPR,SC_EXPR )
my.data = as.data.frame(result)

# save the results so that we can start from here
write.csv(result, file = "finalData.csv")

#*****************************************************************************************************************
# Filtering Step
#*****************************************************************************************************************
# MRMR feauture selection AND select top 2000 genes, threads = 0 means its parallelized
selection.result = MRMR(my.data, targetvector, k = 2000, threads = 0)

# SAVE the selected probes as a csv file as this step is computational expensive
write.csv(selection.result, file = "features.csv")

# read the features
features = as.data.frame(read.csv(file = "features.csv"))

# select the features from the whole dataset
ready_data = my.data[,c(as.factor(features$features))]

# Takes the entrezIDs of the probes of the expression data
geneSymbols = as.character(features$features)

# connect to the database org.Hs.eg.db and fetch the SYMBOLS
annotations = select(org.Hs.eg.db, geneSymbols, "SYMBOL", "ENTREZID")

# Removes the annotations for those genes for which it didn't find the EntrezID
annotations = annotations[!is.na(annotations$ENTREZID),]
annotations

# Set colnames of the dataframe to the gene symbols instead of numbers
colnames(ready_data) = annotations$SYMBOL

# if a gene symbol is not fetched, drop the probe
ready_data <- ready_data[!is.na(names(ready_data))]

# append the label vector to the dataframe
ready_data$class = class

# check if the naming is done
colnames(ready_data)

#******************************************************************************************************************
# Correlation checking
##*****************************************************************************************************************
# set a seed to ensure the results are repeatable
set.seed(7)

# load the data
dim(ready_data)

# calculate correlation matrix, drop the label out of te calculation
correlationMatrix <- cor(ready_data[ ,1:2000])

# summarize the correlation matrix
print(correlationMatrix)

# find attributes that are highly corrected (ideally >0.75)
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)

# print indexes of highly correlated attributes
print(highlyCorrelated)

# print the names of highly correlated genes
print(colnames(ready_data[,highlyCorrelated]))


#*****************************************************************************************************************
# Random Forests Step on the 2000 genes
#*****************************************************************************************************************
# import these two libraries to multithread
library(parallel)
library(doParallel)

# convention to leave 1 core for OS
cluster <- makeCluster(detectCores() - 1) 

# register a cluster
registerDoParallel(cluster)

# grid search for the optimal number of parameters at each split
mtry <- sqrt(ncol(ready_data))
tunegrid <- expand.grid(.mtry=mtry)

# ensure results are repeatable
set.seed(7)

# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=3, allowParallel = TRUE)

# train the model
model <- train(class~., 
                    data=ready_data, 
                    method='rf', 
                    metric='Accuracy', importance = TRUE, 
                    tuneGrid=tunegrid, 
                    trControl=control)

# terminate the reserved cluster
stopCluster(cluster)
registerDoSEQ()
print(model)

# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)

importance

# Open a PDF for plotting; units are inches by default
pdf("/Users/husseinalilakkis/Desktop/CAPSTONE_DATA/importance_var_RF.pdf", width=80, height=1000)

# plot importance of variables
plot(importance)

# Close the PDF file's associated graphics device (necessary to finalize the output)
dev.off()

#************************************************************************************************************
# Recursive Feature Elimination, not in the report, takes so much time
#************************************************************************************************************
# define the control using a random forest selection function

control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(ready_data[,1:2000], ready_data[,2001], sizes=c(1:20), rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))

#************************************************************************************************************
# Principle Component Analysis
#************************************************************************************************************
# set a split
## 60% of the sample size
smp_size <- floor(0.6 * nrow(ready_data))

## set the seed to make your partition reproducible
set.seed(123)
#create a list of training indices
train_ind <- sample(seq_len(nrow(ready_data)), size = smp_size)

#PCA with scaling on the first 2000 feature which are genes
# exclude labels
prin_comp <- prcomp(ready_data[,1:2000], scale = TRUE)

# check standard deviation
std_dev <- prin_comp$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
sum(prop_varex[1:40])

dim(ready_data)

# create a new dataframe using the pc as variables 
data <- data.frame(class = ready_data$class, prin_comp$x)

# split the data
train.data <- data[train_ind ,]
test.data <- data[-train_ind,]

#*****************************************************************************************************************
# SVM on PC
#*****************************************************************************************************************

# import the library used for SVM
library(e1071)

# lopping to find the optimal number of pc to include
for(i in seq(2, 100, by = 1)){
  
  # subset columns ie. iterate over PC 
  train.data <- data[train_ind ,1:i]
  test.data <- data[-train_ind,1:i]
  
  svm = svm(as.factor(class) ~ ., data = train.data, method = "C-classificaton", 
            kernel = "linear", 
            cost = 1, cross = 10)
  
  list[[i]] <- svm$tot.accuracy
  
}
# Plot the bar chart.
plot(list,type = "o",col = "black", 
     xlab = "Number of PC", ylab = "Accuracy", 
     main = "Variation of Accuracy")

train.data <- data[train_ind ,1:7]
test.data <- data[-train_ind,1:7]

# train the model on the already benchmarked parameters
svm = svm(as.factor(class) ~ ., data = train.data, 
          method = "C-classificaton", 
          kernel = "linear", cost = 1,
          cross = 10)


# get a summary from the training
summary(svm)

# do the prediction
svm.predict = predict(svm, test.data)
svm.predict = as.data.frame(svm.predict)
svm.predict$svm.predict

# set a confusion matrix and get the statistics
confusionMatrix(as.factor(svm.predict$svm.predict),as.factor(test.data$class), mode = "prec_recall")

 #*****************************************************************************************************************
# RF on PC
 #*****************************************************************************************************************
library(randomForest)

#Create control function for training with 10 folds and keep 3 folds for training. search method is grid.
# mtry is the number of predictors at each split, it's usually the sqrt
mtry <- sqrt(ncol(train.data))
tunegrid <- expand.grid(.mtry=mtry)

# set the control
control <- trainControl(method='repeatedcv', 
                         number=10, 
                         repeats=3,
                         allowParallel = TRUE)
 
# allow parallel and register A CLUSTER
cluster <- makeCluster(detectCores() - 1) 
registerDoParallel(cluster)

# train the forest on 7 pc
rf <- train(class~., 
            data=train.data, 
            method='rf', 
            metric='Accuracy', importance = TRUE, 
            trControl=control,
            tunegrid = tune.control())

# plot the accuracy vs mtry 
plot(rf)

rf$results

# do the predictions 
rf.predict = predict(rf, test.data)
rf.predict = as.data.frame(rf.predict)
rf.predict$rf.predict
#generate a confusion matrix
confusionMatrix(as.factor(rf.predict$rf.predict),as.factor(test.data$class), mode = "prec_recall")

#*****************************************************************************************************************
# gbt on PC
#*****************************************************************************************************************
## 10-fold CV training control
fitControl <- trainControl(method = "cv",
                           number = 10,
                           allowParallel = TRUE)

# set the tuning grid 
gbmGrid <-  expand.grid(interaction.depth = c(1:9), 
                        n.trees = (1:30)*50, 
                        shrinkage = 0.1,
                        n.minobsinnode = 20)
# register a cluster
cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
registerDoParallel(cluster)

# set seed so that it is reproducible
set.seed(825)

# train the model and do grid search
gbmFit1 <- train(class ~ ., data = train.data, 
                 method = "gbm", 
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = FALSE,
                 tuneGrid = gbmGrid)


stopCluster(cluster)
registerDoSEQ()

# plot the results of the grid search in a caret theme
trellis.par.set(caretTheme())
plot(gbmFit1) 
trellis.par.set(caretTheme())
plot(gbmFit1, metric = "Accuracy", plotType = "level",
     scales = list(x = list(rot = 90)))

# get the best tune and the best model
gbmFit1$bestTune
sorted = sort(gbmFit1$results[,5])

sorted

# do the predictions and generate a confusion matrix
gbm.predict = predict(gbmFit1, test.data)
gbm.predict = as.data.frame(gbm.predict)
gbm.predict$gbm.predict

confusionMatrix(as.factor(gbm.predict$gbm.predict),as.factor(test.data$class), mode = "prec_recall")

#************************************************************************************************************
#*****************************************************************************************************************
# rda on PC, THIS IS EXTRA BUT I WONT INCLUDE IT IN THE REPORT
#*****************************************************************************************************************
set.seed(825)

fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10,
  allowParallel = TRUE,
  classProbs = TRUE)
rdaFit <- train(class ~ ., data = train.data, 
                method = "rda", 
                trControl = fitControl, 
                tuneLength = 4,
               metric = "ROC")


rdaFit$finalModel      

rda.predict = predict(rdaFit, test.data)
rda.predict = as.data.frame(rda.predict)
rda.predict$rda.predict

confusionMatrix(as.factor(rda.predict$rda.predict),as.factor(test.data$class), mode = "prec_recall")

#************************************************************************************************************
#*****************************************************************************************************************
# KNN on PC
#*****************************************************************************************************************
set.seed(400)

# set the training control
ctrl <- trainControl(method="repeatedcv",repeats = 10, number = 10) 
#,classProbs=TRUE,summaryFunction = twoClassSummary)
# do grid search for the best number of neighbors
knnFit <- train(class ~ ., data = train.data, method = "knn", 
                trControl = ctrl, tuneLength = 20)

# plot accuracy vs num of neighbors
plot(knnFit)

knnFit$bestTune

# do the predictions
knn.predict = predict(knnFit, test.data)
knn.predict = as.data.frame(knn.predict)
knn.predict$knn.predict

#generate a confusion matrix
conf = confusionMatrix(as.factor(knn.predict$knn.predict),as.factor(test.data$class), mode = "prec_recall")
conf

#************************************************************************************************************
#*****************************************************************************************************************
# NB on PC
#*****************************************************************************************************************

set.seed(400)
# set the controls
ctrl <- trainControl(method="repeatedcv",repeats = 5, allowParallel = TRUE) 

# set y as label (class) and predictors as x
x <- train.data[,2:7]
y <- train.data[,1]

# set the grid, use parametric and non parametric methods (kernels)
search_grid <- expand.grid(usekernel = c(TRUE, FALSE),
                           fL = 0:5,
                           adjust = seq(0, 5, by = 1))

# train the naive bayes model
nbFit <- train(x = x,y = y,method = "nb",
               trControl = ctrl,
               tuneGrid = search_grid)
plot(nbFit)

# do the the predictions on the test set
nb.predict = predict(nbFit, test.data)
nb.predict = as.data.frame(nb.predict)
nb.predict$nb.predict

#generate a confusion matrix
conf = confusionMatrix(as.factor(nb.predict$nb.predict),as.factor(test.data$class), mode = "prec_recall")
conf

#************************************************************************************************************
#*****************************************************************************************************************
# NNs on PC
#*****************************************************************************************************************
# Train on entire training set.
# training <- data

# set the controls
numFolds <- trainControl(method = 'repeatedcv', number = 10, classProbs = TRUE, verboseIter = TRUE)

# set the grid for search with size being the number of hidden layers, and choose best decay
nnetGrid <-  expand.grid(size = seq(from = 1, to = 10, by = 1),
                         decay = seq(from = 0.1, to = 0.7, by = 0.1))

# train the neural net, use ROC as metric
nnFit <- train(class ~ . -class, 
              data = train.data,
              method = 'nnet', 
              preProcess = c('center', 'scale'), 
              metric = "ROC",
              trControl = numFolds, 
              tuneGrid = nnetGrid)

plot(nnFit)
nnFit$bestTune

# do the predictions on the test data
nn.predict = predict(nnFit, test.data)
nn.predict = as.data.frame(nn.predict)
nn.predict$nn.predict

#generate a confusion matrix

conf = confusionMatrix(as.factor(nn.predict$nn.predict),as.factor(test.data$class), mode = "prec_recall")
conf

# check the probaility of the predictions
probs <- predict(nnFit, newdata=test.data, type='prob')
plot(probs)


