## Step 7 methylation databases
To analyze the methylation pattern of *ACE2* gene, we select BECon: A tool for interpreting DNA methylation findings from blood in the context of brain, and DiseaseMeth version 2.0 - The human disease methylation database.

The aim of BECon (Blood-Brain Epigenetic Concordance) is to allow for improved interpretation of surrogate methylation results by looking at the relationship human blood and brain methylation.
<img src= "./images/BECON.PNG">
Spearman correlation values of methylation between blood and the listed brain region. In these columns colors are used to show the percentile of all correlations each value reaches (see legend). Darker grey and orange represent lower or higher correlations, repetitively. The cell composition metrics are the delta betas between data unadjusted for cell composition and data adjusted for cell composition. So how much the beta values change on average at a CpG with cell composition adjustment.

To acess the methylation pattern of *ACE2* in a cancer context, we used DiseaseMeth. The Search Engine was redesigned into three query entrance: GeneSearch (query DNA methylation level for disease associated genes by disease type), DiseaseSearch (query DNA methylation level across diseases for a specified gene/genomic position), AdvanceSearch (query disease methylome or DNA methylation level for a specified gene/genomic position).

* Disease option is used to select the interested diseases. We select all the cancer types and diseases in order to compare the results.
* Gene name/Transcript ID/Genomic Interval are used to set interested genomic regions. Users can input multiple items. We select *ACE2* only.
* Technology Experimental Platform option is used to select the datatype of methylome. We select all the methylation platforms.
* The option of controls from the same tissue/cellline is used the select all controls from differential researches but the same tissue/cellline. The option of controls from the same research is to select controls in the same research with cases. In some researches, maybe there is no controls matched with cases.
* Method of Differential Analysis: select the method for case-control differential analysis. There are 4 methods including t-test, minfi, samr, edgeR. If t-test is selected, the differential analysis will be applied through Student's t test (we used the default t test).
* Significant p-value is used to set the significant p-value for the analysis. The significant results will be marked with red color in the results panel. The options were 0.01 and 0.05 (we select 0.05).
* Absolute methylation difference: this option is used to filter the differentially methylated region together with the p-value.
* The result will display the outcome of differential analysis, disease-gene association and methlation profile of the given region in the selected diseases. Differential analysis will be applied in the TSS or genomic region by t-test between cases and controls. Survival analysis and functional annotation can also be applied. The Gene-Gene Relationship Analysis and Disease-Disease Relationship Analysis will display the Pearson correlation coefficient between genes or diseases.

Differential methylated results - case vs control : Here we provide 3 examples of the query results:
<img src= "./images/dismeth-lung.PNG">
<img src= "./images/dismeth-lusc.PNG">
<img src= "./images/dismeth-neuroblastoma.PNG">

Disease gene association:
<img src= "./images/dis-meth-dise-geneass.PNG">

Survival analysis: at this moment, the survival analysis was fail
<img src= "./images/dishmethsurvival.PNG">
<img src= "./images/deuruim.PNG">

**Note** - the error is also for the survival analysis module.

[INTRODUCTION](./index.md)
[Previous - Step 6](./page6.md) 
[Next - Step 8](./page8.md)
