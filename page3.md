## Step 3 searching of ACE2 Gene expression

To evaluate the gene expression of *ACE2* gene, we choose two databases: ESCAPE and GTEx portal.
In ESCAPE, we have the options RNAi screens, protein lists from IP-MS pull-downs, genes differentially expressed after knock-down or over-expression, and target genes for transcription factors and histone modifications as determined by ChIP-seq. The input can be gene-lists for overlap with gene lists from the ESCAPE database. On the left, users can cut and paste lists of Entrez gene symbols and then press Submit to perform the enrichment analysis. In the middle, most of the lists from the ESCAPE database are visualized on a canvas. Each square represents a list. The color indicates the experiment type, and the brightness indicates the level of local similarity among the lists. The enriched terms appear as circles on top of the colored squares representing the gene lists from the ESCAPE database on the canvas: the brighter the circle the more significant the overlap with the input list. The results are also available in a table with the associated p-values on the right (Figure 2).

<img src= "./images/escape-output.PNG"> 

Additionally, in the bottom of gene list, we found the EnrichR option, which provide gene terms related to Transcription, Pathways, Ontologies, Disease/Drugs, Cell types, and other resources. Figure 3 exemplifies a Pathway result of the gene input generated in ESCAPE:

<img src= "./images/escape-enrichr.PNG"> 

The Genotype-Tissue Expression (GTEx) project is an ongoing effort to build a comprehensive public resource to study tissue-specific gene expression and regulation. Samples were collected from 54 non-diseased tissue sites across nearly 1000 individuals, primarily for molecular assays including WGS, WES, and RNA-Seq. Remaining samples are available from the GTEx Biobank. The GTEx Portal provides open access to data including gene expression, QTLs, and histology images.

The initial screen show several options for input data:
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

In the option "Browse By gene ID", we can visualize the expression of *ACE2* across more than 30 tissues (Figure 4) 

<img src= "./images/gene-exp-plot.svg"> 

We also visualize the Exon Expression of *ACE2*:

<img src= "./images/gtex-exon-expression.PNG"> 
<img src= "./images/gtex-gene-model-exon.PNG"> 

In the **Expression** - Multi-Gene Query - the input as a gene list, we chose compare the expression of the *ACE* and *ACE2* genes:
<img src= "./images/gtex-multiquerygene.PNG"> 

[INTRODUCTION](./index.md)        [Previous - Step2](./page2.md)             [Next - Step4](./page4.md)
