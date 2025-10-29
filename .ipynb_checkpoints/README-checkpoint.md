# DNMT3A_CH

This repository contains the supporting code for the paper: 

**"DNMT3A-R882 mutations intrinsically drive dysfunctional neutropoiesis from human haematopoietic stem cells"**

The preprint of paper is available at: https://www.biorxiv.org/content/10.1101/2025.10.08.679225v1.full

The codes are organised into folders that correspond to the analyses present in the paper.

1. scRNA-seq analysis of *in-vitro* myeloid differentiation
    * Author: G.  Mantica
    * The scripts are numbered following the analysis pipeline:
        - 00 preprocessing
        - 01 QC and data transformation
        - 02 clustering, label transfer, and pseudotime
        - 03 differential gene expression, GSVEA, GSVA
        - 04 Cell Chat

2. bulk RNA-seq analysis of *in-vivo* myeloid differentiation
    * Author: G. Mantica
    * The scripts are numbered following the analysis pipeline

3. DNMT3A genotyping of single-cell derived colonies
    * Author: A. Vedi
    * The scripts are organised in folders, the input is included in each folder

4. phenotyping analysis of single-cell derived colonies
    * Author: G. Mantica
    * The scripts are numbered following the analysis pipeline 

5. multiparameter flow cytometry analysis of single-cell derived colonies with Cytotree
    * Author: A. Krzywon

6. Shannon diversity and global methylation analysis on phylogenetic trees
    * Author: G. Mantica

The data is currently hosted at:

1. scRNA-seq: GEO GSE298242 

2. bulk RNA-seq: GEO GSE306029

3. variant calling results for genotyping: supplementary table S2 https://www.biorxiv.org/content/10.1101/2025.10.08.679225v1.full

4. phenotyping input (output from FlowJo analysis): supplementary table S2 https://www.biorxiv.org/content/10.1101/2025.10.08.679225v1.full

6. WGS files and phylogenetic tree: EGAS00001004280,EGAD00001015750 ; Mendeley data: DOI:10.17632/rcv6tkvbfy.1