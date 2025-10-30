# %% [markdown]
# # 01.00_QC_normalisation_PCA
# Author:GM <br>

# %%
suppressPackageStartupMessages({library(tidyverse)
                                library(RColorBrewer)
                                library(edgeR)
                                library(gridExtra)
                                library(PCAtools)
                                })

# %% [markdown]
# R version: R version 4.3.2 (2023-10-31) 
#
#   Package      Version
#   <chr>        <chr>  
# 1 edgeR        4.0.16 
# 2 gridExtra    2.3    
# 3 PCAtools     2.14.0 
# 4 RColorBrewer 1.1-3  
# 5 tidyverse    2.0.0

# %%
basedir <- paste(laurenti, user, project, sep='/')
basename <- '_CHIP30_'
print(basedir)

# %%
outdir<-paste(basedir, 'output', sep='/') 
figdir<- paste(outdir, 'figures', '01.02/', sep='/')

# %% [markdown]
# ### import data

# %%
setwd('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/output/tables')

# %%
v2 <- read.csv('202503_CHIP30_countmatrix_v2_clean.txt')

# %% [markdown]
# ### select only coding genes and remove MT and X Y genes

# %% [markdown]
# #### 0. extract gene lengths

# %%
gene_length_v2 <- v2[, c(1,2)]
row.names(v2)<-v2[,1]
v2<-v2[,-c(1,2)]

# %%
colnames(v2) <- sapply(colnames(v2), function(x) gsub("_.*", "", x))
colnames(v2) <- sapply(colnames(v2), function(x) gsub("HUE7565", "", x))

# %% [markdown]
# #### 1. gene info to get coding genes

# %%
gene_info <- read.delim(paste0(basedir,'/ref/geneInfo.tab'), header = FALSE, stringsAsFactors = FALSE)

# %%
gene_info <-gene_info[-1,]

# %%
colnames(gene_info)<-c('ENS_id', 'Geneid', 'type')

# %%
setdiff(v2$Geneid, gene_info$Geneid) #all the genes are present

# %%
coding<- gene_info %>% filter(type=='protein_coding')

# %% [markdown]
# #### 2. gtf file to get chromosome location

# %%
gtf <- read.table(paste0(basedir,'/ref/gencode.v19.annotation.gtf'), sep = "\t", header = FALSE, comment.char = "#", stringsAsFactors = FALSE)

# %%
colnames(gtf)<-c('chromosome', 'V2', 'V3','start','end','V6','V7','V8','V9')

# %%
gtf$gene_name <- sapply(gtf$V9, function(x) gsub('"', '', sub(".*gene_name ([^;]+);.*", "\\1", x)))

# %%
setdiff(coding$Geneid, gtf$gene_name) #all the genes appear in the gtf file

# %%
gene_chr <- gtf %>% select(chromosome, gene_name)

# %%
gene_chr_somatic <- gene_chr %>% filter(!(chromosome %in% c('chrX','chrY','chrM')))

# %% [markdown]
# #### 3. filter for coding and somatic genes only

# %%
v2_filt <- v2[row.names(v2) %in% coding$Geneid,]
v2_filt <- v2_filt[row.names(v2_filt) %in% gene_chr_somatic$gene_name,]

# %%
v2_filt %>% dim()

# %% [markdown]
# ### clean up meta for exports

# %%
meta<-read.csv(paste0(basedir,'/meta/RNASeq_Laurenti_Sample table.csv'))

# %%
meta$SAMPLE_ID <- sapply(meta$SAMPLE_ID, function(x) gsub("HUE7565", "", x))

# %%
col_order <- meta$SAMPLE_ID

# %%
v2_filt <- v2_filt[,col_order]

# %%
meta<-write_csv(meta, paste0(basedir,'/meta/', prefix, 'meta_cleaned.csv'))

# %%
paste0(basedir,'/meta/', prefix, 'meta_cleaned.csv')

# %% [markdown]
# #### 4. export raw matrix only coding and somatic genes

# %%
v2_filt_exp <- cbind(gene = rownames(v2_filt), v2_filt)

# %%
# write_csv(v2_filt_exp, paste0(outdir, '/tables/01.02/', prefix, '_v2_rawcounts_somaticcodinggenes.csv'))

# %% [markdown]
# ### raw counts QC overview

# %% [markdown]
# #### 1. lowly expressed genes filtering

# %%
v2.bin.10= ifelse(v2_filt>=10, yes=1,no=0)

# %%
v2.bin.3= ifelse(v2_filt>=3, yes=1,no=0)

# %%
v2.bin.10 %>% rowMeans %>% density() %>% plot(main='v2', xlab='rowmeans thresh10') 
abline(v = 0.2, col = "red", lty = 2, lwd = 4)

# %%
v2.bin.3 %>% rowMeans %>% density() %>% plot(main='v2', xlab='rowmeans thresh3') 
abline(v = 0.2, col = "red", lty = 2, lwd = 4)

# %%
dim(v2.bin.10)
sum(rowMeans(v2.bin.10) >= 0.20)

# %% [markdown]
# ##### 1. filter for genes exp >10 in 80% of the samples

# %%
keep.v2 <- rowMeans(v2_filt >= 10) >= 0.20
v2.mtx <- v2_filt[keep.v2, ]

# %% [markdown]
# #### 2. sample qc

# %%
plot(colSums(v2.mtx) + 1, log = "y", main='v2', pch=19)
text(x = seq_along(totals), y = totals + 1, labels = colnames(v2.mtx), pos = 3, cex = 0.9)

# %% [markdown]
# #### 3. export filtered matrix

# %%
v2.mtx_exp <- cbind(gene = rownames(v2.mtx), v2.mtx)

# %%
# write_csv(v2.mtx_exp, paste0(outdir, '/tables/01.02/', prefix, '_v2_rawcounts_somaticcodinggenes_genefiltered.csv'))

# %% [markdown]
# ### logTPM normalisation

# %%
y2 <-DGEList(v2.mtx)

# %%
meta$group <- paste(meta$GENOTYPE, meta$CELL_TYPE, meta$TREATMENT, sep = "_")

# %%
group<-meta[meta$SAMPLE_ID %in% rownames(y2$samples),] %>% pull(group)
mouse <- meta[meta$SAMPLE_ID %in% rownames(y2$samples),] %>% pull(MOUSE_ID)
genotype <- meta[meta$SAMPLE_ID %in% rownames(y2$samples),] %>% pull(GENOTYPE)
VAF <- meta[meta$SAMPLE_ID %in% rownames(y2$samples),] %>% pull(VAF)
cell_type<-meta[meta$SAMPLE_ID %in% rownames(y2$samples),] %>% pull(CELL_TYPE)
treatment<-meta[meta$SAMPLE_ID %in% rownames(y2$samples),] %>% pull(TREATMENT)

# %%
y2$samples$group <- as_factor(group) #not sure I will keep this group for DGE
y2$samples$mouse <- as_factor(mouse)
y2$samples$genotype <-  as_factor(genotype)
y2$samples$cell_type <- as_factor(cell_type)
y2$samples$treatment <- as_factor(treatment)
y2$samples$VAF <- as_factor(VAF) #to change if we ever get the continuous values

# %%
print("Dimensions v2 without subsetting:") 
print(dim(y2))
print("")
keep <- filterByExpr(y2) 
print("Dimensions v2 after subsetting:") 
print(dim(y2))

# %%
y2 <- calcNormFactors(y2)

# %%
rownames(gene_length_v2)<- gene_length_v2$Geneid
gene_length_v2 <- gene_length_v2[rownames(y2$counts), ]

# %%
v2.rpkm <- rpkm(y2, gene.length = gene_length_v2$Length) 

# %%
# Convert RPKM to TPM:

v2.TPM <- ((v2.rpkm) / colSums(v2.rpkm)[col(v2.rpkm)] )*10^6

# %%
v2.TPM.mtx <- cbind(gene = rownames(v2.TPM), v2.TPM)

# %%
v2.TPM.mtx<-as_tibble(v2.TPM.mtx)

# %%
write_csv(v2.TPM.mtx, paste0(outdir, '/tables/01.02/', prefix, '_v2.TPM.mtx.csv'))

# %% [markdown]
# #####  log tx TPM

# %%
v2.log_mtx <- log2(v2.TPM+1)

# %%
v2.log_mtx %>% head()

# %%
v2.log_mtx <- cbind(gene = rownames(v2.log_mtx), v2.log_mtx)

# %%
v2.log_mtx<-as_tibble(v2.log_mtx)

# %%
write_csv(v2.log_mtx, paste0(outdir,'/tables/01.02/',prefix,'_v2.log2norm.TPM.mtx.csv'))

# %% [markdown]
# ### TMM normalisation and PCA

# %%
#adjust meta to project on PCA
meta<-meta[-7]
rownames(meta)<-meta$SAMPLE_ID
meta$v2.libsize <- colSums(v2.mtx)
meta$v3.libsize <- colSums(v3.mtx)

# %%
meta.var <- colnames(meta)
meta.var <- meta.var[-c(1)]

# %%
dge_TMM.v2<- DGEList(counts=v2.mtx)
dge_TMM.v2<-calcNormFactors(dge_TMM.v2, method = 'TMM')
edger.df.v2 <-cpm(dge_TMM.v2, log=TRUE, prior.count=4)

# %%
var_genes.v2 <- apply(edger.df.v2 , 1, var)
select_var.v2 <- names(sort(var_genes.v2, decreasing=TRUE))[1:1000]
topvar.genes.v2<-list(select_var.v2)

# %%
topvar_data.v2 <- edger.df.v2[select_var.v2,]

# %%
p.v2 <- pca(topvar_data.v2, metadata = meta, center = T, scale = T)
elbow.v2 <-findElbowPoint(p.v2$variance)
elbow.v2 

# %%
perc_var_3pc.v2<-sum(p.v2$variance[1:3])
perc_var_3pc.v2

# %%
fig(6,7)
p<-screeplot(p.v2, components = getComponents(p.v2, components = c(1:15)), vline = c(elbow.v2))
p
#ggsave(paste0(figdir,prefix,'_screeplot_v2.pdf'),p, w=6, h=7)

# %% [markdown]
# ####  PCA boxplot

# %%
PCA<-p.v2$rotated[,1:2]

# %%
PCA$SAMPLE_ID <- row.names(PCA)

# %%
PCA<-left_join(PCA, meta, by='SAMPLE_ID')

# %%
PCA$CELL_TYPE <- factor(PCA$CELL_TYPE, levels=c('PRE', 'NEU'))
PCA$GENOTYPE <- factor(PCA$GENOTYPE, levels=c('WT', 'R882H'))

# %%
fig(5,6)
p5 = ggplot(PCA %>% filter(TREATMENT=='CTRL'), aes(x= GENOTYPE, y=PC1, fill= GENOTYPE)) + 
        geom_boxplot(outlier.shape = NA) + 
        geom_jitter(height = 0, width = 0.2, alpha = 0.5 , size=2.4, color='black') + 
        ylab('PCA1 values') + 
        xlab('') + 
        theme_classic(base_size=20) +
        theme(legend.position='none', axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_fill_brewer(palette="Set1") +
        facet_grid(. ~ CELL_TYPE, scales = "free_x", space = "free_x") +
        theme(strip.background = element_blank(),
                            strip.text = element_text(size = 16, face = "bold"),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           legend.position='none')
        ylim(-45,28)

p5

# %%
ggsave(paste0(figdir,prefix,'_PCA1_v2_boxplot.pdf'),p5, w=5, h=6)

# %%
PC1_wt_neu<- PCA %>% filter(TREATMENT=='CTRL', CELL_TYPE=='NEU', GENOTYPE=='WT') %>% pull(PC1)
PC1_mut_neu<- PCA %>% filter(TREATMENT=='CTRL', CELL_TYPE=='NEU', GENOTYPE=='R882H') %>% pull(PC1)


# %%
t.test(PC1_wt_neu,PC1_mut_neu)

# %%
wilcox.test(PC1_wt_neu,PC1_mut_neu)

# %% [markdown]
# ####  PCA plots

# %%
for (i in 1:length(meta.var)){
pairplot.p.v2<-pairsplot(p.v2, components = getComponents(p.v2, seq_len(8)), colby = meta.var[i], pointSize = 1.7, plotaxes=F, trianglelabSize=8)
legend<-pairsplot(p.v2, components = getComponents(p.v2, seq_len(2)), colby = meta.var[i], pointSize = 1.7, plotaxes=F, trianglelabSize=6, legendPosition = "left", legendLabSize = 10)
pdf(paste0(figdir,prefix,'_PCA_v2_', meta.var[i], ".pdf"), height = 14, width = 12)
print(pairplot.p.v2)
print(legend)
dev.off()
}

# %%
#plots <- list()
for (i in 1:length(meta.var[-7])){
p<-biplot(p.v2,
    colby = meta.var[i],
    hline = 0, vline = 0,
    legendPosition = 'right')

assign(paste0('p', i), p)
}

# %%
fig(16,20)
pca.v2 <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 3)
pca.v2
#ggsave(paste0(figdir,prefix,'_PCA12_v2.pdf'),pca.v2, w=16, h=20)

# %%
fig(7,6)
PCA_v2_merged <- biplot(p.v2,
       lab=NA,
       colby = 'GENOTYPE',
       shape = 'CELL_TYPE',
       legendPosition = "right" )
PCA_v2_merged

# %%
ggsave(paste0(figdir,prefix,'_PCA12_v2_merged.pdf'),PCA_v2_merged, w=7, h=6)

