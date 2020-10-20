## Step 9 miRNA, lncRNA and circRNA Databases

We used one miRNA database, one lncRNA and circRNA database and one viral miRNA database in order to improve our results.

1) miRNA - miRDB - miRDB is an online database for miRNA target prediction and functional annotations. All the targets in miRDB were predicted by a bioinformatics tool, MirTarget, which was developed by analyzing thousands of miRNA-target interactions from high-throughput sequencing experiments. Common features associated with miRNA binding and target downregulation have been identified and used to predict miRNA targets with machine learning methods. miRDB hosts predicted miRNA targets in five species: human, mouse, rat, dog and chicken. Users may also provide their own sequences for custom target prediction using the updated prediction algorithm. In addition, through combined computational analyses and literature mining, functionally active miRNAs in humans and mice were identified. These miRNAs, as well as associated functional annotations, are presented in the FuncMir Collection in miRDB. As a recent update, miRDB presents the expression profiles of hundreds of cell lines and the user may limit their search for miRNA targets that are expressed in a cell line of interest. The user can choose one of the following search options: Search by miRNA name OR Search by gene target
The initial search returns 50 miRNA candidates to bind the *ACE2* gene:

<img src= "./images/mirdb1.PNG">

Select the first result, the user can visualize the miRNA target and description:
<img src= "./images/mirdb2.PNG">

The user may acess the miRNA information and related references:
<img src= "./images/mirdb3.PNG">

2) lncRNA and circRNA - exoRBase is a repository of circular RNA (circRNA), long non-coding RNA (lncRNA) and messenger RNA (mRNA) derived from RNA-seq data analyses of human blood exosomes. Experimental validations from published literature are also included. exoRBase features the integration and visualization of RNA expression profiles based on normalized RNA-seq data spanning both normal individuals and patients with different diseases.
The user can visualize the data in three different plots: line chart, heatmap and box plot.

<img src= "./images/exorbase1.png">
<img src= "./images/exorbase-heatmap.png">
<img src= "./images/chart (1).png">

3) Viral miRNA -  VirMir: a database containing predicted viral miRNA candidate hairpins. 
The microRNA hairpin discovery pipeline was also applied to discover viral encoded microRNAs. All viral genomes were obtain from NCBI. The classification of virus is based on the taxonomy table of NCBI (Jun, 2006). Totally, the genomes of 2266 viruses were analyzed. Users may query the putative miRNA hairpins of a specific viral species by the hierarchical menu or by search function using the GenBank Identifier or RefSeq accession number(e.g. 30844336 or NC_003663 for Cowpox virus). In addition, users can also search for the putative target genes of a particular viral miRNA hairpins by a RNAhybrid service link. The 3'-UTR regions of human, mouse, rat, zebrafish, arabidopsis and rice genes are available for search.
Results:
<img src= "./images/virmir-result.PNG">

Exploring Human coronavirus OC43 - first result - 22979
<img src= "./images/virmir-result1.PNG">

Target prediction with RNAhybrid
For each putative or known miRNA, user may query the target sites in the 3'UTR regions of genes of human, mouse, rat, zebrafish, rice or Arabidopsis. From the human.rna.fna file downloaded from NCBI, we fetched the sequences of human 3'UTRs based on the CDS (coding sequence) positions acquired from the human.rna.gbff file. We adopted the same strategy for the mouse, rat, zebrafish, rice and Arabidopsis 3'UTRs.
In our case, the result send by e-mail was empty, which means there's no matching sequence returned from RNAhybrid.
When operating RNAhybrid, the pipeline first calculates the optimal free energy of miRNA at which it binds to a perfectly complementary site, then the miRNA/mRNA duplex mfe (minimum free energy). An alignment whose RNA duplex mfe is more than the selected percentage of its correspondent optimal free energy is regarded as a positive alignment. The alignment result will be sent to the users via their e-mail addresses.


[INTRODUCTION](./index.md)        [Previous - Step 8](./page8.md)     [Next - Step 10](./page10.md)
