# Load the ML packages and set seed (for the random number generator. 
#This ensures that the sequence of random numbers generated will be the same each time you run the code. 
#This is crucial for reproducibility in statistical analysis and simulations, allowing others (or yourself at a later time) to obtain the same results.)

library(caret)
library(DALEX)
library(pROC)
library(DALEXtra)

library(dplyr) 
library(tidyr)

set.seed(123)

# Load the main and meta data
carcinoma_data <- read.csv(file = "raw_expression_data_tcga-coad.csv", header = TRUE)
carcinoma_meta <- read.csv(file = "metadata tcga-coad0.csv", header = TRUE)




# Preprocessing the main data


rownames(carcinoma_data) <- carcinoma_data$X # Making column "X" the row names, rather than numerical "1", "2", "3"...
carcinoma_data$X <- NULL # Removing duplicate columns as row names


boxplot(carcinoma_data, col = "lightblue")
carcinoma_data <- log2(carcinoma_data + 1) # A log transformation to normalize the data
boxplot(carcinoma_data, col = "lightblue")


colnames(carcinoma_data) <- gsub("\\.", "-", colnames(carcinoma_data)) # Editing the format of the main data so it becomes similar to the rownames of the meta data

# Transpose the main data
carcinoma_data <- data.frame(t(carcinoma_data))

SDs = apply(carcinoma_data, 2, sd)
# This line calculates the standard deviation for each column (each gene) in carcinoma_data and stores the results in the vector SDs. 
# A higher standard deviation indicates greater variability in the gene expression levels.

topPredicts = order(SDs, decreasing = T)[1:3000]
#After calculating the standard deviation, this lines ensures that the top 3000 genes with the highest SD are selected, and in descending order.


carcinoma_data = carcinoma_data[, topPredicts]
#The trans_data is now reduced to only include the top 3000 genes that have the highest variability, which are often more informative for subsequent analyses like classification



# Preprocessing the meta data
anyNA(carcinoma_meta)
sum(is.na(carcinoma_meta))
carcinoma_meta <- carcinoma_meta %>% drop_na()
anyNA(carcinoma_meta)



rownames(carcinoma_meta) <- carcinoma_meta$barcode #changed the rownames from numerical "1, 2, 3, 4, 5...." to the Barcode

carcinoma_meta$barcode <- NULL #to remove the row name duplicate


# Merge both data
carcinoma_merged_data <- merge(carcinoma_data, carcinoma_meta, by = "row.names") #merging of both main and meta data

dim(carcinoma_merged_data)
View(carcinoma_merged_data)
rownames(carcinoma_merged_data) <- carcinoma_merged_data$Row.names # to ensure the row names are the samples
carcinoma_merged_data$Row.names <- NULL #to remove duplicate columns of row names


# Step 4: Further preprocessing steps

# To remove near zero variation
all_zero <- preProcess(carcinoma_merged_data, method = "nzv", uniqueCut = 15)
carcinoma_merged_data <- predict(all_zero, carcinoma_merged_data)
dim(carcinoma_merged_data)

# center. Centering helps to stabilize numerical computations and can be particularly important for algorithms sensitive to the scale of the data (like PCA, KNN, etc.).
all_center <- preProcess(carcinoma_merged_data, method = "center")
carcinoma_merged_data <- predict(all_center, carcinoma_merged_data)
dim(carcinoma_merged_data)

# to remove highly correlated values.
#Reduce Multicollinearity: High correlations between features can lead to issues in model training, such as inflated standard errors and difficulties in interpreting coefficients.
#Removing redundant features helps create a more stable model. Improve Model Performance: By reducing the number of features,
# you can improve the efficiency and performance of certain machine learning algorithms, especially those sensitive to multicollinearity.
all_corr <-preProcess(carcinoma_merged_data, method = "corr", cutoff = 0.5)
carcinoma_merged_data <- predict(all_corr, carcinoma_merged_data)
dim(carcinoma_merged_data)


# Step 5: preparing the data for the ML model
# Assuming your data is already loaded in the variable carcinoma_merged_data
data <- carcinoma_merged_data

# One-hot encoding the gender (convert gender into 0 and 1)
data$gender <- ifelse(data$gender == "female", 0, 1)


# Select only the gene expression columns (assuming columns starting with 'ENSG' are gene expression data)
gene_expression_columns <- grep("^ENSG", colnames(data), value = TRUE)

X <- data[, gene_expression_columns]
y <- data$gender

# Step 6: Split the data into training and testing sets
set.seed(123)  # For reproducibility
trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]

# Step 7: Train a KNN model
knn_model <- train(X_train, as.factor(y_train), method = "knn", tuneLength = 5)

# Predict on the test set
y_pred <- predict(knn_model, X_test)
x_pred <- predict(knn_model, X_train)
# Evaluate the model using confusion matrix
conf_matrix <- confusionMatrix(y_pred, as.factor(y_test))
conf_matrix1 <- confusionMatrix(x_pred, as.factor(y_train))
print(conf_matrix)
print(conf_matrix1)
# Step 8: Calculate accuracy
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
print(paste("Accuracy: ", accuracy))
# Now perform permutation importance
set.seed(123)  # For reproducibility
importance <- varImp(knn_model, scale = FALSE)
print(importance)
