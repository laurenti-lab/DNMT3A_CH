# %% [markdown]
# # B_03.03_GSVA on external dataset of cardiovascular disease
# Author:GM <br>
# Data from Pekayvaz, Losert, Knottenberg et al : https://www.nature.com/articles/s41591-024-02953-4

# %%
suppressPackageStartupMessages({library(tidyverse)
                                library(RColorBrewer)
                                library(edgeR)
                                library(GSVA)
                                library(gridExtra)
                                })

# %% [markdown]
# R version: R version 4.3.2 (2023-10-31) 
#   Package      Version
#   <chr>        <chr>  
# 1 edgeR        4.0.16 
# 2 gridExtra    2.3    
# 3 GSVA         1.50.5 
# 4 RColorBrewer 1.1-3  
# 5 tidyverse    2.0.0 

# %%
outdir<-paste(basedir, 'output', sep='/')
figdir<- paste(outdir, 'figures', '06.02/', sep='/')

# %%
options(repr.matrix.max.cols=20) #Inf

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %%
boxpl <- function(data, parameter){
                p <-ggplot(data, aes(x = classification, y = !!sym(parameter), fill = classification))+
                      geom_boxplot(outliers = FALSE) + 
                      geom_dotplot(binaxis = "y", stackdir = "center",
                                    dotsize = 0.4,
                                    binwidth = 0.015,
                                    position = position_jitter(0.3), alpha = 0.5, fill= 'black') +
                      stat_summary(fun = median, geom = "crossbar", width = 0.8, fatten = 0.8, color = "black") +
                      theme_classic(base_size = 16) +
                      labs(x = "", y = parameter) +
                      guides(fill = "none") +
                      scale_fill_brewer(palette = 'Paired') +
                      facet_grid(. ~ measurement, scales = "free_x", space = "free_x") +
                      theme(strip.background = element_blank(),
                            strip.text = element_text(size = 16, face = "bold"),
                            axis.text.x = element_text(angle = 45, hjust = 1),
                            legend.position = "none")
    
                plot_name= paste0('p_',parameter)
                assign(plot_name, p, .GlobalEnv) 
    
                return(p)
                           }

# %% [markdown]
# ### import

# %% [markdown]
# #### data

# %%
counts_raw <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/output/tables/06.01/20250105_neut_primeseq_counts.csv', 
                 stringsAsFactors = FALSE)


# %%
rownames(counts_raw) <- counts_raw$sample_id

# %%
counts_raw <- counts_raw[,-1]

# %%
count_matrix <-as.matrix(t(counts_raw))

# %% [markdown]
# #### BED file hg19

# %%
coding <- read_tsv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/references/gencode.v39.protein_coding.txt')

# %%
coding <- coding[!coding$chr %in% c("chrM", "chrY"),]  #remove chrM and chrY

# %% [markdown]
# #### metadata

# %%
meta_raw <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/references/Pekayvaz2024/Merged_Sample_Meta_Data_CRPclean.csv')

# %%
meta_raw <- meta_raw[,-1]

# %%
meta <- meta_raw %>% select(sample_id, age, sex, classification, group, CRP, Troponin, CK, Subject, measurement)

# %% [markdown]
# #### gene sets

# %% [markdown]
# ##### hallmarks

# %%
hallmark_genesets_path <- "~/giovanna/references/h.all.v2023.2.Hs.symbols.gmt"

# %%
hallmarks <- fgsea::gmtPathways(hallmark_genesets_path)

# %% [markdown]
# ##### DNMT3A Neutrophils (Cluster 6) DGE

# %%
col_names=c('gene', 'logFC', 'logCPM', 'F', 'PValue', 'FDR')

# %%
FDR_thres = 0.05
logFC_thres = 1

# %%
cl6_model1_MUT_vs_WT<-read_tsv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/output/tables/edgeR/20241127_6_model1_untreated_groupMUT-groupWT_DGE_results.txt',
                               skip=1,
                              col_names = col_names)


# %%
cl6_m1_UP <- cl6_model1_MUT_vs_WT %>% filter(FDR < FDR_thres, logFC > logFC_thres) %>% pull(gene)

# %%
hallmarks$NEUT_DNMT3A_UP<- cl6_m1_UP

# %% [markdown]
# ### QC and filtering

# %%
counts_raw %>% rowMeans %>% plot(pch=19)

# %%
counts_raw.bin= ifelse(counts_raw!=0, yes=1,no=0)

# %%
col_means <- colMeans(counts_raw.bin)

data <- data.frame(value = col_means)

ggplot(data, aes(x = "", y = value)) +
  geom_violin(fill = "coral2") +
  geom_hline(yintercept = 0.2, color = "blue", linetype = "dashed", size = 1.5) + 
  theme_classic()

# %%
ggplot(data, aes(x = value)) +  # Use the column means as x
  geom_histogram(binwidth = 0.005, fill = "coral2", color = "black", alpha = 0.7) +  # Histogram
  geom_vline(xintercept = 0.20, color = "blue", linetype = "dashed", size = 1.5) +  # Vertical line at 0.15
  theme_classic() + 
  xlab("Column Means") +  # X-axis label
  ylab("Frequency")  # Y-axis label

# %%
gene_sample_pct <- counts_raw.bin %>% colMeans
gene_sample_pct <- gene_sample_pct[gene_sample_pct>0.2]

# %%
counts_filt1<- counts_raw[,colnames(counts_raw) %in% names(gene_sample_pct)]

# %%
counts_filt1.bin= ifelse(counts_filt1 !=0, yes=1,no=0)

counts_filt1.bin %>% colMeans %>% plot()
abline(h = 0.8, col = "blue", lty = 2, lwd = 2)

# %%
counts_filt1.bin %>% rowMeans %>% plot()
abline(h = 0.25, col = "blue", lty = 2, lwd = 2)

# %%
sample_gene_pct <- counts_filt1.bin %>% rowMeans
sample_gene_pct <- sample_gene_pct[sample_gene_pct>0.25]

# %%
counts_filt2 <- counts_filt1[rownames(counts_filt1) %in% names(sample_gene_pct),]

# %%
counts_filt2_mx <- as.matrix(t(counts_filt2))

# %%
counts_filt2_mx <- counts_filt2_mx[rownames(counts_filt2_mx) %in% coding$gene_name, ]

# %% [markdown]
# ### data transformation

# %% [markdown]
# #### TPM normalisation

# %%
y <- DGEList(counts_filt2_mx)

# %%
print("Dimensions without subsetting:")
print(dim(y))
print("")

# %%
y <- calcNormFactors(y) #library size correction

# %%
coding_filt <- coding[!duplicated(coding$gene_name), ]

# %%
coding_filt <- coding_filt[coding_filt$gene_name %in% rownames(counts_filt2_mx),] %>% arrange(match(gene_name, rownames(counts_filt2_mx)))

# %%
matr.expr.norm.rpkm <- rpkm(y, gene.length = coding_filt$end - coding_filt$start) 

# %%
# Convert RPKM to TPM:
matr.expr.norm.TPM <- ((matr.expr.norm.rpkm) / colSums(matr.expr.norm.rpkm)[col(matr.expr.norm.rpkm)] )*10^6

# %%
matr.expr.norm.TPM %>% head()

# %%
mtx <- cbind(gene = rownames(matr.expr.norm.TPM), matr.expr.norm.TPM)

# %%
mtx<-as_tibble(mtx)

# %%
#write_csv(mtx, paste0(outdir, '/tables/06.02/', prefix, '_NeutPrimeSeq_normalisedTPM_filter2.csv'))

# %% [markdown]
# #### log tx TPM

# %%
log_mtx <- log2(matr.expr.norm.TPM+1)

# %%
log_mtx %>% head()

# %%
log_mtx <- cbind(gene = rownames(log_mtx), log_mtx)

# %%
log_mtx %>% head(n=2)

# %%
log_mtx<-as_tibble(log_mtx)

# %%
#write_csv(log_mtx, paste0(outdir,'/tables/06.02/',prefix,'_NeutPrimeSeq_normalisedTPM_log2normalised_filter2.csv'))

# %% [markdown]
# ##### stopping point, can reimport matrix

# %%
logTPM <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/output/tables/06.02/20250110_NeutPrimeSeq_normalisedTPM_log2normalised_filter2.csv', check.names = FALSE)

# %%
logTPM%>%head()

# %%
rownames(logTPM)<-logTPM$gene

# %%
logTPM<- logTPM[,-1]

# %%
logTPM<-as.matrix(logTPM)

# %% [markdown]
# ### GSVA

# %%
gsvaPar <- gsvaParam(logTPM, hallmarks)

# %%
gsvaPar

# %%
gsva_values <- gsva(gsvaPar, verbose=TRUE)

# %%
gsva_val <- t(gsva_values)

# %% [markdown]
# #### merge with meta

# %%
gsva_val <- cbind(sample_id = rownames(gsva_val), gsva_val)
gsva_val <- as_tibble(gsva_val)
gsva_val[] <- lapply(gsva_val, as.double)
gsva<-left_join(gsva_val, meta, by = 'sample_id')

# %% [markdown]
# ##### No enrichment of hallmark inflammatory response gene set

# %%
fig(3,7)
p2<- ggplot(gsva %>% filter(measurement=='TP0'),aes(x = group, y = HALLMARK_INFLAMMATORY_RESPONSE, fill=group))+
                      geom_boxplot(outliers = FALSE) + 
                       geom_dotplot(binaxis = "y", stackdir = "center",
                                    dotsize = 0.5,
                                    binwidth = 0.015,
                                    position = position_jitter(0.1), alpha = 0.5, fill='black') +
                      stat_summary(fun = median, geom = "crossbar", width = 0.8, fatten = 0.8, color = "black") +
                      theme_classic(base_size = 16) +
                      #labs(x = "", y = parameter) +
                      #guides(fill = "none") +
                      scale_fill_brewer(palette='Paired') +
                      facet_grid(. ~ measurement, scales = "free_x", space = "free_x") +
                      #scale_y_log10() +
                      theme(strip.background = element_blank(),
                            strip.text = element_text(size = 16, face = "bold"),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           legend.position='none')

p2

# %%
l_mod <- lm(HALLMARK_INFLAMMATORY_RESPONSE ~ group + age, data = gsva %>% filter(measurement=='TP0'))
res <- summary(l_mod)
res

# %% [markdown]
# ##### No enrichment of hallmark inflammatory response gene set

# %%
fig(3,7)
p2<- ggplot(gsva %>% filter(measurement=='TP0'),aes(x = group, y = HALLMARK_INTERFERON_GAMMA_RESPONSE, fill=group))+
                      geom_boxplot(outliers = FALSE) + 
                       geom_dotplot(binaxis = "y", stackdir = "center",
                                    dotsize = 0.5,
                                    binwidth = 0.015,
                                    position = position_jitter(0.1), alpha = 0.5, fill='black') +
                      stat_summary(fun = median, geom = "crossbar", width = 0.8, fatten = 0.8, color = "black") +
                      theme_classic(base_size = 16) +
                      #labs(x = "", y = parameter) +
                      #guides(fill = "none") +
                      scale_fill_brewer(palette='Paired') +
                      facet_grid(. ~ measurement, scales = "free_x", space = "free_x") +
                      #scale_y_log10() +
                      theme(strip.background = element_blank(),
                            strip.text = element_text(size = 16, face = "bold"),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           legend.position='none')

p2

# %%
l_mod <- lm(HALLMARK_INTERFERON_GAMMA_RESPONSE ~ group + age, data = gsva %>% filter(measurement=='TP0'))
res <- summary(l_mod)
res

# %% [markdown]
# ##### Enrichment of DNMT3A NEUT gene set

# %%
fig(3,7)
p2<- ggplot(gsva %>% filter(measurement=='TP0'),aes(x = group, y = NEUT_DNMT3A_UP, fill=group))+
                      geom_boxplot(outliers = FALSE) + 
                       geom_dotplot(binaxis = "y", stackdir = "center",
                                    dotsize = 0.5,
                                    binwidth = 0.015,
                                    position = position_jitter(0.1), alpha = 0.5, fill='black') +
                      stat_summary(fun = median, geom = "crossbar", width = 0.8, fatten = 0.8, color = "black") +
                      theme_classic(base_size = 16) +
                      #labs(x = "", y = parameter) +
                      #guides(fill = "none") +
                      scale_fill_brewer(palette='Paired') +
                      facet_grid(. ~ measurement, scales = "free_x", space = "free_x") +
                      #scale_y_log10() +
                      theme(strip.background = element_blank(),
                            strip.text = element_text(size = 16, face = "bold"),
                           axis.text.x = element_text(angle = 45, hjust = 1),
                           legend.position='none')

p2

# %%
l_mod <- lm(NEUT_DNMT3A_UP ~ group + age, data = gsva %>% filter(measurement=='TP0'))
res <- summary(l_mod)
res
