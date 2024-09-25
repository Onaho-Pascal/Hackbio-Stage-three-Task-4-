# Integrating Machine Learning and Differential Expression Analysis for Biomarker Discovery in Colon Adenocarcinoma 🧬💻🤖
 **Authors:** Manal Agbada, Ojiaku Confidence Chisom, Abdulrahman Walid Elbagoury, Rahma Nabil Sallam, Pascal Onaho, Hagar Haitham, Elazab, Mercy Francis, Ariyo Adesokan
 ## Introduction 
Colon adenocarcinoma (COAD) is a prevalent and deadly cancer, often complicated by late diagnosis and treatment resistance [1]. The Cancer Genome Atlas (TCGA) has significantly advanced COAD research by providing extensive omics and clinical data, enhancing our understanding of COAD’s molecular landscape and revealing potential biomarkers essential for improving prognosis and developing personalized treatment strategies [2].
## Dataset Description
For this analysis, we utilized the COAD dataset from TCGA, which includes omics and clinical information. RNA-seq data was downloaded in RStudio using the ‘GDCdownload’ and ‘GDCprepare’ functions from the ‘TCGAbiolinks’ package in Bioconductor. We created a simplified metadata dataset by selecting specific features and reduced the number of samples by selecting 20 female and 20 male samples. Raw expression counts were normalized for gene length and filtered to remove lowly expressed genes. 
## Methodology and Results
### Biomarker
Using the ‘TCGAbiolinks’ library and the ‘edgeR’ package, we performed differential expression analysis (DEA) to identify differentially expressed genes (DEGs) between the two groups. Genes were considered differentially expressed if they had an adjusted p-value < 0.05 and an absolute log fold change > 1. This analysis identified 752 upregulated and 371 downregulated genes. Enrichment analysis (EA) was then performed using the ‘TCGAbiolinks’ library, focusing on Gene Ontology (GO) terms, with results visualized in a bar plot.

![volcanoplot](https://github.com/user-attachments/assets/b8cf0aa3-79d6-4f6e-9c34-704e09dddaef)

![heatmap](https://github.com/user-attachments/assets/cbb8a389-3b1d-4685-a436-396c02bbf272)

### Machine Learning
The K-Nearest Neighbor (KNN) machine learning model was employed to identify potential biomarkers and important genes for predicting gender. In bioinformatics, KNN classifies biological data, such as gene expression profiles and DNA sequences, by comparing an unknown sample to its nearest neighbors in a labeled dataset. The KNN model achieved an overall accuracy of 66.67% in classifying gender, with a high sensitivity of 83.33% for correctly identifying females. However, its specificity for males was only 50%, reflecting a tendency to misclassify males as females, as indicated by a positive predictive value (PPV) of 62.5% for females and a negative predictive value (NPV) of 75% for males. The Kappa statistic of 0.3333 shows moderate agreement between predictions and actual values, while McNemar’s test indicates no significant difference in error types across genders. The model identified key genes associated with gender classification, including RPS4Y1, APOH, CYP1A1, GSTM1, and PF4V1. 

![Screenshot 2024-09-24 at 23 23 51](https://github.com/user-attachments/assets/e667d3a6-00cb-46bb-9945-ccc303ce325f)


## Conclusion
Despite reasonable performance, improvements are needed to enhance male classification. The model was trained using a merged dataset of genetic and meta data, which underwent normalization and preprocessing, before being split into training (70%) and test (30%) data for model development.

## References
1. Hu G, Yao H, Wei Z, Li L, Yu Z, Li J, Luo X, Guo Z. A bioinformatics approach to identify a disulfidptosis-related gene signature for prognostic implication in colon adenocarcinoma. Sci Rep. 2023 Jul 31;13(1):12403.
2. The Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. Nature 487, 330–337 (2012).
3. Agarwal A, Singh K, Kant S, Bahadur RP. A comparative analysis of machine learning classifiers for predicting protein-binding nucleotides in RNA sequences. Comput Struct Biotechnol J. 2022 Jun 17;20:3195-3207..

