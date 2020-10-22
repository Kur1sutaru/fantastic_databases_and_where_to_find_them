## Step 3 Searching of *ACE2* Gene expression

To evaluate the gene expression of the ACE2 gene, we chose two databases: ESCAPE and the GTEx portal. In ESCAPE, there are the options: RNAi screens, protein lists, differentially expressed genes after knock-down or over-expression experiments, and target genes for transcription factors and histone modifications determined by ChIP-seq. The input can be  a list of genes for overlap with gene lists from the ESCAPE database. Users can input lists with Entrez gene symbols and submit to start the enrichment analysis. 

<img src= "./images/escape-output.PNG"> 

Lists from the ESCAPE database are visualized as shown below. On the left there is the input list. On the right are the enrichment results and p-value. The figure in the middle represents the experiments leading to the results in a color coded manner. The circles indicate the presence of input genes and brightness indicates local similarity among the lists.

Additionally, at the bottom of the gene list, we found the EnrichR option, which provides gene terms related to Transcription Factors, Pathways, Ontologies, Disease/Drugs, Cell types, and other resources. Figure 3 exemplifies the result from EnrichR to the gene input generated in ESCAPE.


<img src= "./images/escape-enrichr.PNG"> 

The Genotype-Tissue Expression (GTEx) project is an ongoing effort to build a comprehensive public resource to study tissue-specific gene expression and regulation. Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals, primarily for molecular assays including WGS, WES, and RNA-Seq. Remaining samples are available from the GTEx Biobank. The GTEx Portal provides open access to data including gene expression, QTLs, and histology images.

The initial screen shows several options for input data:
* **Browse**  By gene ID - Browse and search all data by gene
              By variant or rs ID - Browse and search all data by variant
              
* **By Tissue** - 	Browse and search all data by tissue
* **Histology Image Viewer**	- Browse and search GTEx histology images
* **Expression** 
              Multi-Gene Query - Browse and search expression by gene and tissue
              Top 50 Expressed Genes - Visualize the top 50 expressed genes in each tissue
              Transcript Browser	Visualize transcript expression and isoform structures
* **QTL**
              Locus Browser	Visualize QTLs by gene in the Locus Browser
              Locus Browser - Visualize QTLs by variant in the Locus Browser VC (Variant Centric)
              IGV Browser	 - Visualize tissue-specific eQTLs and coverage data in the IGV Browser
              eQTL Dashboard	- Batch query eQTLs by gene and tissue
              eQTL Calculator	- Test your own eQTLs
* **eGTEx**
              H3K27ac	Browse - H3K27ac ChIP-seq data in IGV Browser
* **Biobank**
              Access Biospecimens	Search and request available GTEx biospecimens

In the option "Browse By gene ID", we can visualize the expression of *ACE2* across more than 30 tissues.

<img src= "./images/gene-exp-plot.svg"> 

The expression of *ACE2* gene was enriched in several gastric tissues, in testis and heart, suggesting a important role in the regulation of cardiovascular and gastrointestinal function, as well as in the fertility.

We also visualize the Exon Expression of *ACE2*:

<img src= "./images/gtex-exon-expression.PNG"> 
<img src= "./images/gtex-gene-model-exon.PNG"> 

In the **Expression** - Multi-Gene Query - the input is a gene list, we chose compare the expression of the *ACE* and *ACE2* genes:
<img src= "./images/gtex-multiquerygene.PNG"> 

[INTRODUCTION](./index.md)        [Previous - Step2](./page2.md)             [Next - Step4](./page4.md)
