######################################
#### HackBio internship - stage 3 #### 
######################################

# Biomarkers analysis #### 

## load necessary libraries 

library("TCGAbiolinks")
library(SummarizedExperiment)
library(edgeR)
library(gplots)
library(ggplot2)
library(biomaRt)

## pick any cancer type/subtype and download from TCGA 

# get project information 
getProjectSummary("TCGA-COAD")

# download the dataset
tcga_coad <- GDCquery(project = "TCGA-COAD",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification")
GDCdownload(tcga_coad) 
coad_data <- GDCprepare(tcga_coad) 
head(coad_data) 
View(coad_data)

# explore metadata information
coad_data$barcode

coad_data$race 
table(coad_data$race) 

coad_data$tumor_descriptor 
table(coad_data$tumor_descriptor) 

coad_data$ajcc_pathologic_stage 
table(coad_data$ajcc_pathologic_stage)

coad_data$ajcc_pathologic_m 
table(coad_data$ajcc_pathologic_m)

coad_data$gender
table(coad_data$gender)

# create a simple metadata 
metadata_df <- data.frame("barcode" = coad_data$barcode,
                          "race" = coad_data$race,
                          'tumor_type' = coad_data$tumor_descriptor,
                          'stage' = coad_data$ajcc_pathologic_stage,
                          'metastasis_status' = coad_data$ajcc_pathologic_m,
                          'gender' = coad_data$gender)
View(metadata_df)

# subset the metadata female vs male 
coad_rawdata <- assays(coad_data) 
dim(coad_rawdata$unstranded) 
View(coad_rawdata$unstranded)

# downsize the data to 20 vs 20  
samples <- c(subset(metadata_df, gender == "female")$barcode[c(1:20)],
             subset(metadata_df, gender == "male")$barcode[c(1:20)])
final_df <- coad_rawdata$unstranded[ , c(samples)]
dim(final_df)
View(final_df)

# clean and preprocess the data, handling missing values, normalizing gene expression data ####
table(is.na(final_df)) # no NAs
norm_data <- TCGAanalyze_Normalization(tabDF = final_df, geneInfo = geneInfoHT, method = "geneLength")

filt_data <- TCGAanalyze_Filtering(tabDF = norm_data,
                                   method = "quantile",
                                   qnt.cut = 0.25)


## differential expression analysis 
DEA <- TCGAanalyze_DEA(mat1 = filt_data[ , c(samples)[1:20]], 
                       mat2 = filt_data[ , c(samples)[21:40]],
                       Cond1type = "female",
                       Cond2type = "male",
                       pipeline = "edgeR",
                       fdr.cut = 0.05,
                       logFC.cut = 1)

DEA.Level <- 
  TCGAanalyze_LevelTab(DEA, "female", "male",
                       filt_data[ , c(samples)[1:20]],
                       filt_data[ , c(samples)[21:40]])
View(DEA.Level)

# visualization of top DEGs with a heatmap color coded based on the samples (female = red, male = blue)
heat.DEGs <- filt_data[rownames(DEA.Level), ]

# color code
gender <- c(rep("female", 20), rep("male", 20))
ccodes <- c()
for (i in gender) {
  if ( i == "female") {
    ccodes <- c(ccodes, "red") 
  } else {
    ccodes <- c(ccodes, "blue") 
  }
}

# heatmap
heatmap.2(x = as.matrix(heat.DEGs),
          col = hcl.colors(10, palette = 'Blue-Red 2'),
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.9, cexCol = 0.7,
          main = "Heatmap of Females vs Males",
          na.color = 'black',
          ColSideColors = ccodes,
          margins = c(11,10))
legend("left",                       
       legend = c("Female", "Male"),     
       col = c("red", "blue"),           
       lty = 1,                          
       lwd = 4,                         
       cex = 0.7,
       xpd = TRUE,
       inset = c(-0.2, 2))                      

# volcano plot
ggplot(DEA.Level, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(FDR < 0.05 & abs(logFC) > 1,
                                ifelse(logFC > 1, "Upregulated", "Downregulated"), "Not significant")), size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 FDR", title = "Volcano Plot of Females vs Males", color = "Gene Regulation") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top")



## functional enrichment analysis 

# selection of up- and down-regulated genes from the DEA 
upreg.genes <- rownames(subset(DEA.Level, logFC > 1 & FDR < 0.05))
dnreg.genes <- rownames(subset(DEA.Level, logFC < -1 & FDR < 0.05))

# convert ensemble IDs to gene IDs
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg.genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = upreg.genes,
                     mart = mart)$hgnc_symbol

dnreg.genes <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                     filters = 'ensembl_gene_id',
                     values = dnreg.genes,
                     mart = mart)$hgnc_symbol

# EA for up- and down-regulated genes 
up.EA <- TCGAanalyze_EAcomplete(TFname = "Upregulated", upreg.genes) 
dn.EA <- TCGAanalyze_EAcomplete(TFname = "Downregulated", dnreg.genes)

# barplot for enriched pathways in up-regulated genes 
TCGAvisualize_EAbarplot(tf = rownames(up.EA$ResBP), 
                        GOBPTab = up.EA$ResBP, 
                        GOCCTab = up.EA$ResCC,
                        GOMFTab = up.EA$ResMF,
                        PathTab = up.EA$ResPat, 
                        nRGTab = upreg.genes, 
                        nBar = 5, 
                        text.size = 2,
                        fig.width = 30,
                        fig.height = 15)

# barplot for enriched pathways in down-regulated genes
TCGAvisualize_EAbarplot(tf = rownames(dn.EA$ResBP), 
                        GOBPTab = dn.EA$ResBP, 
                        GOCCTab = dn.EA$ResCC,
                        GOMFTab = dn.EA$ResMF,
                        PathTab = dn.EA$ResPat, 
                        nRGTab = dnreg.genes, 
                        nBar = 5, 
                        text.size = 2, 
                        fig.width = 30,
                        fig.height = 15)



# Machine learning model development #### 

## load necessary libraries and set seed (for random number generation)

#This ensures that the sequence of random numbers generated will be the same each time you run the code. 
#This is crucial for reproducibility in statistical analysis and simulations, allowing others (or yourself at a later time) to obtain the same results.)

library(caret)
library(DALEX)
library(pROC)
library(DALEXtra)
library(dplyr) 
library(tidyr)

set.seed(123)


## Preprocessing the data
carcinoma_data <- final_df

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
carcinoma_meta <-  metadata_df

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


# Preparing the data for the ML model
# Assuming your data is already loaded in the variable carcinoma_merged_data
data <- carcinoma_merged_data

# One-hot encoding the gender (convert gender into 0 and 1)
data$gender <- ifelse(data$gender == "female", 0, 1)


# Select only the gene expression columns (assuming columns starting with 'ENSG' are gene expression data)
gene_expression_columns <- grep("^ENSG", colnames(data), value = TRUE)

X <- data[, gene_expression_columns]
y <- data$gender

# Split the data into training and testing sets
set.seed(123)  
trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[trainIndex, ]
X_test <- X[-trainIndex, ]
y_train <- y[trainIndex]
y_test <- y[-trainIndex]


## Train a KNN model
knn_model <- train(X_train, as.factor(y_train), method = "knn", tuneLength = 5)

# Predict on the test set
y_pred <- predict(knn_model, X_test)
x_pred <- predict(knn_model, X_train)

# Evaluate the model using confusion matrix
conf_matrix <- confusionMatrix(y_pred, as.factor(y_test))
conf_matrix1 <- confusionMatrix(x_pred, as.factor(y_train))
print(conf_matrix)
print(conf_matrix1)


# Now perform permutation importance
set.seed(123) 
importance <- varImp(knn_model, scale = FALSE)
print(importance)
