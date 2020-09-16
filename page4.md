## Step 4 searching in Cancer databases

To evaluate the expression pattern of *ACE2* in cancer cell lines, we used CCLE - Cancer Cell Line Encyclopedia. The aims of the database were to characterize of a large panel of human cancer models, to develop integrated computational analyses that link distinct pharmacologic vulnerabilities to genomic patterns and to translate cell line integrative genomics integrative genomics into cancer patient stratification.
When search by *ACE2*, the database screen show the Entrez ID (59272), the description (angiotensin I converting enzyme 2), and several external Links. Below, the interactive plot shows the distribution by lineage (RNA-seq expression) and the option to download the data (figure 5)       
<img src= "./images/ccle.PNG">

Other options are: scatter plots of the query gene, with the option to compare with other differentially expressed genes across the available datasets, mutation data (users with previously login only), Fusions/Translocations data, and methylation - CpG viewer.

In order to compare the cell line information with clinical data, we chose the TCGA data portal (via GDC data portal, <https://portal.gdc.cancer.gov/>).
Here we can been search for a specific gene, Cases by Major Primary Site, or for a specific barcode (e.g., TCGA-56-7731, primary site: Bronchus and lung).
The search by *ACE2* showed the Gene Summary, External References,  Cancer Distribution of Single Nucleotide Variants (SNV) and Copy Number Variation (CNV), functional impact of the variants with the option to trade the gene transcripts, and the most frequent somatic mutations (figure 6). The user can be download all of these results in tables or figures.
<img src= "./images/tcga1.PNG">
<img src= "./images/tcga2.PNG">
<img src= "./images/tcga3.PNG">
<img src= "./images/tcga4.PNG">
In the exploration menu, we select "TCGA-LUAD" (LUNG ADENOMAS AND ADENOCARCINOMAS) and obtain 37 cases, with clinical and survival data, mutations and oncogrid (figure 7) 
<img src= "./images/tcga-survivalace2.PNG">
<img src= "./images/tcga-survival.PNG">
<img src= "./images/tcga-luad-oncogrid.PNG">
