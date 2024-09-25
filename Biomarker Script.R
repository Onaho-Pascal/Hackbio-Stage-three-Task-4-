######################################
#### HackBio internship - stage 3 #### 
######################################

# load necessary libraries

library("TCGAbiolinks")
library(SummarizedExperiment)
library(edgeR)
library(gplots)
library(ggplot2)
library(biomaRt)

# pick any cancer type/subtype and download from TCGA ####

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

# differential expression analysis ####
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


# functional enrichment analysis ####

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
