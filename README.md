# DNMT3A_CH

This repository contains the supporting code for the paper: 

**"DNMT3A-R882 mutations intrinsically drive dysfunctional neutropoiesis from human haematopoietic stem cells"**

The preprint of paper is available at: https://www.biorxiv.org/content/10.1101/2025.10.08.679225v1.full

The **code** is organised into folders that correspond to the respective analyses present in the paper.

1. scRNA-seq analysis of *in-vitro* myeloid differentiation
    * Author: G.  Mantica
    * The scripts are numbered following the analysis pipeline:
        - 00 preprocessing
        - 01 QC and data transformation
        - 02 clustering, label transfer, and pseudotime
        - 03 differential gene expression, GSEA, GSVA
        - 04 Cell Chat
        - 05 Annotated Heatmap of DEGs in mature neutrophil cluster 

2. bulk RNA-seq analysis of *in-vivo* myeloid differentiation
    * Author: G. Mantica
    * The scripts are numbered following the analysis pipeline:
        - 00 preprocessing
        - 01 QC, data transformation, PCA
        - 02 differential gene expression, GSEA

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

We acknowledge Hugo P. Bastos for the overall supervision of the bioinformatic work

The source code present in this repository is under the MIT license whereas the datasets (from genotyping of single-cell derived colonies) are under the CC-BY-4.0 license.

The **datasets** associated with this work are available at:

1. scRNA-seq: GEO [GSE298242](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298242) 

2. bulk RNA-seq: GEO [GSE306029](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE306029)

3. variant calling results for genotyping: [supplementary table S2](https://www.biorxiv.org/content/biorxiv/early/2025/10/08/2025.10.08.679225/DC1/embed/media-1.xlsx?download=true)

4. phenotyping input (output from FlowJo analysis): [supplementary table S2](https://www.biorxiv.org/content/biorxiv/early/2025/10/08/2025.10.08.679225/DC1/embed/media-1.xlsx?download=true)

6. WGS files: [EGAS00001004280](https://ega-archive.org/studies/EGAS00001004280), [EGAD00001015750](https://ega-archive.org/datasets/EGAD00001015750) ; and phylogenetic tree: [DOI:10.17632/rcv6tkvbfy.1](https://data.mendeley.com/datasets/rcv6tkvbfy/1)
