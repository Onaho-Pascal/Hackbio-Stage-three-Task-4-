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

![Heatmap](https://github.com/user-attachments/assets/b82af36c-ca2b-43fa-b59a-17462a42ae06)
![Caption1](https://github.com/user-attachments/assets/fec828ed-1d55-4032-ba42-f187681b25af)

![Volcano_plot_DEGs_female_vs_male](https://github.com/user-attachments/assets/3fdcea58-bf21-41c5-a639-39dee41332fd)
![Caption2](https://github.com/user-attachments/assets/ec9d42ec-6837-4683-a67d-dced5d3fb6cf)

![EA_upregulated_female_vs_male](https://github.com/user-attachments/assets/cce6c98f-4fa1-4c06-a079-cdf62ef35da5)
![Caption3](https://github.com/user-attachments/assets/487e9150-07d0-4c9f-be02-f20b1f885e79)

![EA_downregulated_female_vs_male](https://github.com/user-attachments/assets/78463bfd-bea6-4647-ade2-fb4d47cad95a)
![Caption4](https://github.com/user-attachments/assets/e0f12ec4-ade0-4588-bcf2-0b6a174e3c76)

### Machine Learning
The K-Nearest Neighbor (KNN) model, which is an ML model that classifies biological data, such as gene expression profiles and DNA sequences, by comparing an unknown sample to its nearest neighbors in a labeled dataset [3], was applied to identify potential biomarkers and important genes for predicting gender where 0 represents female and 1 represents male. However, the model showed poor performance on the test data, achieving only 27.27% accuracy. It correctly classified 2 females but misclassified 6 males as females, with low sensitivity (25%) and specificity (33.33%). A negative Kappa value (-0.2941) indicated poor agreement between predictions and actual values. The model's high training accuracy of 82% suggests overfitting, as it struggled to generalize to new data, likely due to the small sample size of 40, with a high misclassification, particularly of males as females.  Additionally, key genes identified include RPS4Y1, APOH, CYP1A1, GSTM1, and PF4V1.

## Conclusion
This project successfully integrated machine learning and differential expression analysis to explore potential biomarkers in colon adenocarcinoma. While differential expression analysis identified key upregulated and downregulated genes, the K-Nearest Neighbor model's poor performance on the test data, likely due to a small sample size and potential overfitting, highlights the need for larger datasets and further optimization to improve predictive accuracy. Future research should focus on refining models and exploring additional biomarkers to enhance the detection and treatment of COAD.
## References
1. Hu G, Yao H, Wei Z, Li L, Yu Z, Li J, Luo X, Guo Z. A bioinformatics approach to identify a disulfidptosis-related gene signature for prognostic implication in colon adenocarcinoma. Sci Rep. 2023 Jul 31;13(1):12403.
2. The Cancer Genome Atlas Network. Comprehensive molecular characterization of human colon and rectal cancer. Nature 487, 330–337 (2012).
3. Agarwal A, Singh K, Kant S, Bahadur RP. A comparative analysis of machine learning classifiers for predicting protein-binding nucleotides in RNA sequences. Comput Struct Biotechnol J. 2022 Jun 17;20:3195-3207..

