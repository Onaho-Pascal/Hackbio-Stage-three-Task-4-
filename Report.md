# Integrating Machine Learning and Differential Expression Analysis for Biomarker Discovery in Colon Adenocarcinoma 🧬💻🤖

 **Authors:** Manal Agdada (@Manal), Ojiaku Confidence Chisom (@Areta), Abdulrahman Walid Elbagoury (@Willeau), Rahma Nabil Sallam (@rahmanabil2002), Pascal Onaho (@PascalOnaho), Hagar Haitham Elazab (@HBONH33), Mercy Francis (@Mercylee), Ariyo Adesokan (@Adesokan_ariyo1)

**R script:** https://github.com/Onaho-Pascal/Hackbio-Stage-three-Task-4-/blob/main/Code/stage3_script.R

**Figures:** https://github.com/Onaho-Pascal/Hackbio-Stage-three-Task-4-/tree/main/Figures

 ## Introduction 
Colon adenocarcinoma (COAD) is a prevalent and deadly cancer, often complicated by late diagnosis and treatment resistance [1]. The Cancer Genome Atlas (TCGA) has significantly advanced COAD research by providing extensive omics and clinical data, enhancing our understanding of COAD’s molecular landscape and revealing potential biomarkers essential for improving prognosis and developing personalized treatment strategies [2].

## Dataset Description and preprocessing steps
For this analysis, we utilized the COAD dataset from TCGA, which includes omics and clinical information. RNA-seq data was downloaded in RStudio using the `GDCdownload` and `GDCprepare` functions from the `TCGAbiolinks` package in Bioconductor. We created a simplified metadata dataset by selecting specific features and reduced the number of samples by selecting 20 female and 20 male samples. Raw expression counts were normalized for gene length and filtered to remove lowly expressed genes. 

## Methodology and Results

### Biomarker discovery
Using the `TCGAbiolinks` library and the `edgeR` package, we performed differential expression analysis (DEA) to identify differentially expressed genes (DEGs) between the two groups. Genes were considered differentially expressed if they had an adjusted p-value < 0.05 and an absolute log fold change > 1. This analysis identified 752 upregulated and 371 downregulated genes. Enrichment analysis (EA) was then performed using the `TCGAbiolinks` library, focusing on Gene Ontology (GO) terms, with results visualized in a bar plot.

![heatmap](https://github.com/Onaho-Pascal/Hackbio-Stage-three-Task-4-/blob/main/Figures/Heatmap_DEGs_female_vs_male.png)
![Picture1](https://github.com/user-attachments/assets/fec828ed-1d55-4032-ba42-f187681b25af)

![Volcanoplot](https://github.com/user-attachments/assets/d93d2999-870e-45a9-a645-904ca8dabf90)
![Picture2](https://github.com/user-attachments/assets/ec9d42ec-6837-4683-a67d-dced5d3fb6cf)

![EA_up](https://github.com/Onaho-Pascal/Hackbio-Stage-three-Task-4-/blob/main/Figures/EA_upregulated_female_vs_male.pdf)
![Picture3](https://github.com/user-attachments/assets/7ce12ea5-b393-49fd-9012-33121a13f652)

![EA_dn](https://github.com/Onaho-Pascal/Hackbio-Stage-three-Task-4-/blob/main/Figures/EA_downregulated_female_vs_male.pdf)
![Picture4](https://github.com/user-attachments/assets/e0f12ec4-ade0-4588-bcf2-0b6a174e3c76)

### Machine Learning
The K-Nearest Neighbor (KNN) machine learning model was employed to identify potential biomarkers and important genes for predicting gender. In bioinformatics, KNN classifies biological data, such as gene expression profiles and DNA sequences, by comparing an unknown sample to its nearest neighbors in a labeled dataset [3]. The KNN model achieved an overall accuracy of 66.67% in classifying gender, with a high sensitivity of 83.33% for correctly identifying females. However, its specificity for males was only 50%, reflecting a tendency to misclassify males as females, as indicated by a positive predictive value (PPV) of 62.5% for females and a negative predictive value (NPV) of 75% for males. The Kappa statistics of 0.3333 shows moderate agreement between predictions and actual values, while McNemar’s test indicates no significant difference in error types across genders. The model identified key genes associated with gender classification, including RPS4Y1, APOH, CYP1A1, GSTM1, and PF4V1. 

## Conclusion
Despite reasonable performance, improvements are needed to enhance male classification. The model was trained using a merged dataset of genetic and meta data, which underwent normalization and preprocessing, before being split into training (70%) and test (30%) data for model development.

## References
1. Hu G, Yao H, Wei Z, Li L, Yu Z, Li J, Luo X, Guo Z. A bioinformatics approach to identify a disulfidptosis-related gene signature for prognostic implication in colon adenocarcinoma. Sci Rep. 2023 Jul 31;13(1):12403.
2. The Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. Nature 487, 330–337 (2012).
3. Agarwal A, Singh K, Kant S, Bahadur RP. A comparative analysis of machine learning classifiers for predicting protein-binding nucleotides in RNA sequences. Comput Struct Biotechnol J. 2022 Jun 17;20:3195-3207..

