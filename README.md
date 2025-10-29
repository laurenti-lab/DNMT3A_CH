# DNMT3A_CH

This repository contains the supporting code for the paper: 

**"DNMT3A-R882 mutations intrinsically drive dysfunctional neutropoiesis from human haematopoietic stem cells"**

The files are organised into five folders, corresponding to the five analyses present in the paper.

1. scRNA-seq analysis of *in-vitro* myeloid differentiation
* The scripts are numbered following the analysis pipeline:
- 00 preprocessing
- 01 QC and data transformation
- 02 clustering and label transfer
- 03 differential gene expression

* Author: G.  Mantica

2. bulk RNA-seq analysis of *in-vivo* myeloid differentiation
* The scripts are numbered following the analysis pipeline
* Author: G. Mantica

3. DNMT3A genotyping of single-cell derived colonies
* The scripts are organised in folders, the input is included in each folder
* Author: A. Vedi

4. phenotyping analysis of single-cell derived colonies
* The scripts are numbered following the analysis pipeline 
* Author: G. Mantica

5. multiparameter flow cytometry analysis of single-cell derived colonies with Cytotree
* Author: A. Krzywon

6. Shannon diversity and global methylation analysis on phylogenetic tree
* Author: G. Mantica

The data is currently hosted at:

1. scRNA-seq: GEO GSE298242 

2. bulk RNA-seq: GEO GSE306029

3. variant calling results: supplementary table S2 https://www.biorxiv.org/content/10.1101/2025.10.08.679225v1.full

4. phenotyping input (output from FlowJo analysis): supplementary table S2 https://www.biorxiv.org/content/10.1101/2025.10.08.679225v1.full

6. WGS files and phylogenetic tree: EGAS00001004280,EGAD00001015750 ; Mendeley data: DOI:10.17632/rcv6tkvbfy.1
