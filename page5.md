## Step 5 Comparative databases

We used TISSUES database to compare the *ACE2* gene expression across human and animal model data.
<img src= "./images/TISSUES-MOUSE.PNG">
<img src= "./images/TISSUES-RAT.PNG">
<img src= "./images/TISSUES-HUMAN.PNG">
<img src= "./images/TISSUES-PIG.PNG">

TISSUES also provides experimental curated data, text mining interactions, and orthologs and paralogs. The confidence of each association is signified by stars, where ★★★★★ is the highest confidence and ★☆☆☆☆ is the lowest. 

For the next step, we search a list of genes which interacts with *ACE2* using the GeneMania web <https://genemania.org/>. The search return the following gene list: *ACE2, GK, 
ACE, AGT, CTD-2501B8.1, TMEM27, AGTR2, REN, SULT1E1, RP11-668G10.2, GK2, BDKRB2, EIF4G2, ENPEP, PTP4A3, WEE1, GK5, GHRL, JUND, TUB, SERPIND1, FGGY, CMA1, AQP8* .
We use these gene list with input into TopCluster web server <https://toppcluster.cchmc.org/>. The next step show a feature to explore the data (e.g., gene ontology, drugs, Pubmed, related diseases) and statistical options, such as Correction methods (we choose Bonferroni) and p-value cutoff (0.05). Besides, we may choose start the analysis in the
 Functional Enrichment module or Interaction Network. The database also provide the genes found in the gene list. In our case, 3 genes was not found in the database. We select the top score interactions for run the next step of the analysis. 
 
Note: the quality of the image generated is unclear.

<img src= "./images/data.png">

[INTRODUCTION](./index.md) [Previous - Step 4](./page4.md) [Next - Step 6](./page6.md)
