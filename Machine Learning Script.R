# Load the ML packages and set seed (for the random number generator.
#This ensures that the sequence of random numbers generated will be the same each time you run the code.
#This is crucial for reproducibility in statistical analysis and simulations, allowing others (or yourself at a later time) to obtain the same results.
library(caret)
library(DALEX)
library(pROC)
library(dplyr)
library(tidyr)
# Step 1 -----> Working with merged_data.csv which was named in the enviroment(log_transformed file) and the file was uploaded via environment
barcode<-merged_data$barcode
genders<-merged_data$gender
# Step 2
SDs = apply(merged_data, 2, sd)
# This line calculates the standard deviation for each column (each gene) in trans_data and stores the results in the vector SDs.
# A higher standard deviation indicates greater variability in the gene expression levels.
topPredicts = order(SDs, decreasing = T)[1:3000]
#After calculating the standard deviation, this lines ensures that the top 2000 genes with the highest SD are selected, and in descending order.
log_data_transposed = merged_data[, topPredicts]
log_transformed_data <- cbind(barcode = barcode, log_data_transposed)
#The trans_data is now reduced to only include the top 2000 genes that have the highest variability,
# which are often more informative for subsequent analyses like classification
# Merge the data
carcinoma_merged_data <- merge(log_transformed_data, carcinoma_meta, by = "barcode") #merging of both main and meta data
dim(carcinoma_merged_data)
View(carcinoma_merged_data)
log_transformed_data <- merge(metadata.tcga.coad,log_transformed_data,by= "barcode")
# To remove near zero variation
all_zero <- preProcess(log_transformed_data, method = "nzv", uniqueCut = 15)
log_transformed_data <- predict(all_zero, log_transformed_data)
dim(log_transformed_data)
# center. Centering helps to stabilize numerical computations
#and can be particularly important for algorithms sensitive to the scale of the data (like PCA, KNN, etc.).
all_center <- preProcess(log_transformed_data, method = "center")
carcinoma_merged_data <- predict(all_center, log_transformed_data)
dim(carcinoma_merged_data)
# to remove highly correlated values.
#Reduce Multicollinearity: High correlations between features can lead to issues in model training,
#such as inflated standard errors and difficulties in interpreting coefficients.
#Removing redundant features helps create a more stable model.
# Improve Model Performance: By reducing the number of features,
# you can improve the efficiency and performance of certain machine learning algorithms, especially those sensitive to multicollinearity.
all_corr <-preProcess(carcinoma_merged_data, method = "corr", cutoff = 0.5)
carcinoma_merged_data <- predict(all_corr, carcinoma_merged_data)
dim(carcinoma_merged_data)
# Splitting into training and Data Sets
# Assuming your data is already loaded in the variable carcinoma_merged_data
data <- carcinoma_merged_data
# One-hot encoding the gender (convert gender into 0 and 1)
data$gender_encoded <- ifelse(data$gender == "female", 0, 1)
# Encode tumor stage using label encoding (you can use one-hot encoding if preferred)
data$tumor_stage_encoded <- as.numeric(factor(data$tumor_stage))
# Select only the gene expression columns (assuming columns starting with 'ENSG' are gene expression data)
gene_expression_columns <- grep("^ENSG", colnames(data), value = TRUE)
X <- data[, gene_expression_columns]
y <- data$gender_encoded
# Split the data into training and testing sets
set.seed(123)  # For reproducibility
trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]
# Train a KNN model
knn_model <- train(X_train, as.factor(y_train), method = "knn", tuneLength = 5)
# Predict on the test set
y_pred <- predict(knn_model, X_test)
# Evaluate the model using confusion matrix
conf_matrix <- confusionMatrix(y_pred, as.factor(y_test))
print(conf_matrix)
# Calculate accuracy
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(paste("Accuracy: ", accuracy))
# Now perform permutation importance
set.seed(123)  # For reproducibility
importance <- varImp(knn_model, scale = FALSE)
print(importance)
