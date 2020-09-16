## Step 8 Regulatory Databases
Here we used the Interferome db and TTSMI (Triplex Target DNA Site Mapping and Integration).
The Interferome db - users guide has the link broken (PDF not found). The parameters are: 
* Interferon Type: (ALL) subtype: (ALL)
* Treatment Concentration: Any Treatment Time: Any 
* In Vivo / In Vitro: All Species: ALL
* System:	(Cell or tissues): ALL
* Sample Types:	Any (Normal AND Abnormal)
* Fold Change option - Up: 2.0   Down: 2.0  (The Fold Change value must be greater or equal to 1)
* Gene Symbol List: *ACE, ACE2* 

The output show several options: Search Conditions, Gene Summary, Experiment Data, Ontology Analysis, TF Analysis, Chromosome, IFN Type, and Basal Expression
The Experiment data option:

<img src= "./images/interferome-gene-experiment.PNG">

Transcription Factor provided by TRANSFAC database

<img src= "./images/interferome-tf.PNG">

The location of transcription factors on each of the genes from the region spaning -1500bp to + 500 bp from the start site. Each coloured box represents a specific transcription factor, the key for which is provided below the graphic. The user can move the cursor over the transciption factor's coloured box to reveal the TRANSFAC* predicted match between the transcription factor and its predicted binding site.

Interferon type: The Venn diagram shows the number of genes regulated by one or more IFN type (Type I, II or III). It should be noted that there are far more datasets for genes regulated by type I than for types II or III, this imbalance introduces the risk of false nagatives and bias for the under-represented types II and III; caution is therefore encouraged in interpreting low or negative results from these types.

<img src= "./images/interferome-ifntype.PNG">

TTSMI (TTS Mapping and Integration) database allows the user to facilitate in (i) searching of the TTS using several search terms including genomic location, gene ID, NCBI RefSeq ID, TTS ID, gene symbol and gene description keywords, (ii) interactive filtering of the TTS co-localized with several other gene regulatory signals, (iii) exploring of dynamic combination of structural and functional annotations of specific TTS, and (iv) viewing of the TTS simultaneously with diverse annotation tracks via UCSC genome browser link. The user can quickly search and filter the TTS candidates from our largest collection of unique TTS location. In terms of finding therapeutic targets, the specificity of TTS is the most important factor to control off-target effects. 

Search by *ACE2* with all default parameters, the website was not return any results.
<img src= "./images/ttsmi.PNG">
