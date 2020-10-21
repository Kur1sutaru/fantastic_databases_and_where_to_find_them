## Step 4 searching in Cancer databases

To evaluate the expression pattern of *ACE2* in cancer cell lines, we used CCLE - Cancer Cell Line Encyclopedia. The aims of the database are to characterize a large panel of human cancer models, to develop integrated computational analyses that link distinct pharmacologic vulnerabilities to genomic patterns and to translate cell line integrative genomics into cancer patient stratification.
When search by *ACE2*, the database screen shows the Entrez ID (59272), the description (angiotensin I converting enzyme 2), and several external links. Below, the interactive plot shows the distribution by lineage (RNA-seq expression) and the option to download the data.  

<img src= "./images/ccle.PNG">

Other options are: scatter plots of the query gene, with the option to compare with other differentially expressed genes across the available datasets, mutation data (users with previously login only), fusion and translocation data, and methylation - CpG viewer.

To compare the cell line information with clinical data, we chose the TCGA data portal (via the GDC data portal, https://portal.gdc.cancer.gov/). Here we can search for a specific gene, cases by major primary site, or a specific barcode (e.g., TCGA-56-7731, primary site: Bronchus and lung). The search by ACE2 showed the Gene Summary, External References, Cancer Distribution of Single Nucleotide Variants (SNV), and Copy Number Variation (CNV), the functional impact of the variants with the option to visualize the gene transcripts, and the most frequent somatic mutations (figure 6). The user can download all of these results in tables or figures. In the exploration menu, we selected "TCGA-LUAD" (LUNG ADENOMAS AND ADENOCARCINOMAS) and obtained 37 cases, with clinical and survival data, mutations, and oncogrid.

<img src= "./images/tcga1.PNG">
<img src= "./images/tcga2.PNG">
<img src= "./images/tcga3.PNG">
<img src= "./images/tcga4.PNG">

In the exploration menu, we select "TCGA-LUAD" (LUNG ADENOMAS AND ADENOCARCINOMAS) and obtain 37 cases, with clinical and survival data, mutations and oncogrid. The pins represents the mutations across the *ACE2* gene.

<img src= "./images/tcga-survivalace2.PNG">
<img src= "./images/tcga-survival.PNG">
<img src= "./images/tcga-luad-oncogrid.PNG">

[INTRODUCTION](./index.md) [Previous - Step3](./page3.md) [Next - Step5](./page5.md)
