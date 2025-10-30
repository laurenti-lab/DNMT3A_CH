# %% [markdown]
# # 02.00_EdgeRDGE_GSEA
# Author:GM <br>

# %%
suppressPackageStartupMessages({library(tidyverse)
                                library(RColorBrewer)
                                library(edgeR)
                                library(fgsea)
                                library(gridExtra)
                                })

# %% [markdown]
# R version: R version 4.3.2 (2023-10-31) 
#   Package      Version
#   <chr>        <chr>  
# 1 edgeR        4.0.16 
# 2 gridExtra    2.3    
# 3 RColorBrewer 1.1-3  
# 4 tidyverse    2.0.0 
# 5 fgsea        1.28.0 

# %%
basedir <- paste(laurenti, user, project, sep='/')
basename <- '_CHIP30_'
outdir<-paste(basedir, 'output', 'tables', '01.03', sep='/') 
figdir<- paste(outdir, 'figures', '01.03/', sep='/')

# %%
options(repr.matrix.max.cols=28) #Inf

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %% [markdown]
# ## import 

# %% [markdown]
# ### data

# %%
setwd(paste0(basedir, '/output/tables/01.02/'))

# %%
v2 <- read.csv('20250411_v2_rawcounts_somaticcodinggenes_genefiltered.csv')

# %%
rownames(v2)<- v2$gene
v2<-v2[,-1]

# %% [markdown]
# ### metadata

# %%
meta <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/meta/20250411meta_cleaned.csv')

# %% [markdown]
# ### gene lists

# %%
hallmark_genesets_path <- "~/giovanna/references/h.all.v2023.2.Hs.symbols.gmt"
reactome_genesets_path <- "~/giovanna/references/c2.cp.reactome.v2023.2.Hs.symbols.gmt"

# %%
hallmarks <- fgsea::gmtPathways(hallmark_genesets_path)
reactome <- fgsea::gmtPathways(reactome_genesets_path)

# %%
length(hallmarks)

# %% [markdown]
# #### import genes from Hackert et al
# get to their "core inflammatory programme" as performed in the paper 
# Also export it because it's a bit ridicoulous that I re-compute this every time

# %%
hackert_allgenes <- read_delim('~/giovanna/references/CHIP27/Hackert2023_neut_inflamm.txt', delim="\t")

# %%
hackert_allgenes %>% head(n=50)

# %%
hackert_allgenes <- hackert_allgenes %>% arrange(fisher_adjusted)

# %%
hackert_allgenes %>% head(n=50)

# %% [markdown]
# the table was already ranked

# %%
hackert_top500<- hackert_allgenes %>% slice_head(n = 500)

# %%
hackert_top500 %>% dim()

# %%
hackert_selected <- hackert_top500  %>%
                    filter(abs(mean_lfc) >=0.5)


# %%
hackert_selected_up <- hackert_selected %>%
                       filter(mean_lfc > 0)

# %%
hackert_selected_down <- hackert_selected %>%
                       filter(mean_lfc < 0)


# %%
write.table(hackert_selected_up$symbol,file = "~/giovanna/references/CHIP27/Hackert2023_neut_inflamm_GM_UP.txt", 
            quote = FALSE, row.names = FALSE,col.names = FALSE)
write.table(hackert_selected_down$symbol,file = "~/giovanna/references/CHIP27/Hackert2023_neut_inflamm_GM_DOWN.txt", 
            quote = FALSE, row.names = FALSE,col.names = FALSE)

# %%
hallmarks$HACKERT_NEUTS_INFLAMMATION_UP <- hackert_selected_up$symbol 

# %% [markdown]
# #### import genes from Van Galen et al 
#

# %%
AML_VanGalen_monocytes <- c(
  "ACAP2", "ACTG1", "ADA2", "ADAR", "AHNAK", "AHR", "AIF1", "ALDH2", "ALDOA", 
  "ANKRD28", "ANXA1", "AOAH", "APLP2", "APOBEC3A", "APOL6", "ARHGAP30", 
  "ARID4A", "ARPC5", "ATF4", "ATP2B1-AS1", "ATP5G2", "ATP8B4", "AZU1", "BCL6", 
  "BIN2", "BTG1", "C5orf24", "C6orf48", "CALR", "CAPZB", "CAT", "CCNL1", 
  "CCPG1", "CD163", "CD300E", "CD36", "CD4", "CD47", "CD74", "CDC42SE1", 
  "CDKN2D", "CFD", "CHD2", "CIRBP", "CITED2", "CLEC12A", "CLU", "CNBP", 
  "CNN2", "CNOT2", "COX4I1", "CPVL", "CR1", "CREBRF", "CRIP1", "CST3", 
  "CTSB", "CTSD", "CTSG", "CTSS", "CX3CR1", "CXCL8", "CXXC5", "CYBA", "DBET", 
  "DDIT3", "DDX17", "DDX5", "DMXL2", "DNASE2", "DOK2", "DSE", "DUSP1", 
  "DUSP22", "DUSP6", "EEF1A1", "EEF1B2", "EEF1G", "EEF2", "EFHD2", "EGR1", 
  "EIF1", "EIF3E", "EIF3F", "EIF3K", "EIF4A2", "ELANE", "ELF1", "EMB", 
  "EMP3", "EPSTI1", "EVI2A", "F13A1", "FAM46A", "FAU", "FCGR2A", "FCGR3A", 
  "FCGRT", "FGD2", "FGL2", "FGR", "FKBP5", "FLT3", "FOS", "FTH1P3", "FTL", 
  "FUT4", "FYB1", "GAPDH", "GAS5", "GBP1", "GBP2", "GBP4", "GGA1", "GGNBP2", 
  "GIMAP4", "GIMAP7", "GPX1", "GRB2", "GRN", "GSTK1", "GTF3A", "HCK", 
  "HCLS1", "HINT1", "HLA-A", "HLA-B", "HLA-C", "HLA-DPA1", "HLA-DPB1", 
  "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DRB6", "HLA-E", "HLA-F", "HMGN3", 
  "HMOX1", "HOXA9", "HSD17B11", "HSPA5", "IER2", "IFI16", "IFI30", "IFI44L", 
  "IFI6", "IFITM3", "IL10RA", "IRF2BP2", "ISG15", "IST1", "ITGA4", "ITGAX", 
  "ITGB2", "ITM2B", "JAML", "JUN", "JUNB", "KLF2", "KLF6", "KLHL24", 
  "LGALS2", "LGALS3", "LMNA", "LOC100506675", "LOC100996385", "LPAR6", 
  "LRP1", "LRRC75A-AS1", "LST1", "LUC7L3", "LY6E", "LYST", "MAFB", "MALAT1", 
  "MARCH1", "MARCKS", "MBNL1", "METTL7A", "MFSD1", "MGST1", "MRC1", 
  "MRPL33", "MS4A10", "MS4A6A", "MS4A7", "MSRB3", "MTSS1", "MX1", "MX2", 
  "N4BP2L2", "NCF1", "NCF2", "NCOA4", "NEAT1", "NFKBIA", "NFKBIZ", 
  "NOTCH2NL", "NPC2", "NPM1", "NUDT16", "NUFIP2", "OAS1", "PABPC1", "PAK2", 
  "PARP14", "PARP9", "PCNP", "PDK4", "PFDN5", "PFN1", "PLAC8", "PLD3", 
  "PLEK", "PLEKHO1", "PLIN2", "PNISR", "PNRC1", "POU2F2", "PPP1R15A", 
  "PRR13", "PSAP", "PTEN", "PTP4A2", "PTPN6", "PYCARD", "RAB12", "RAB13", 
  "RACK1", "RGS18", "RHOB", "RICTOR", "RIOK3", "RIPOR2", "RN7SK", "RNA5-8S", 
  "RNASET2", "RNF213", "RNU4-2", "RNU4ATAC", "RNU5F-1", "RSRP1", "RUNX3", 
  "SCIMP", "SCPEP1", "SERP1", "SERPINA1", "SESN3", "SF1", "SFPQ", 
  "SH3BGRL3", "SLC1A3", "SMAP2", "SMIM3", "SNORD79", "SP110", "SQSTM1", 
  "SREK1", "SRGN", "SRRM2", "SSH2", "STAT1", "TAGLN2", "TBK1", "TLR4", 
  "TMEM107", "TMEM70", "TNFAIP2", "TNFRSF1B", "TNFSF10", "TPT1", "TRIM38", 
  "TXNIP", "TYMP", "TYROBP", "UBA52", "UBB", "UBC", "UBE2J1", "UQCRB", 
  "VAMP8", "VMP1", "VSIR", "VTRNA1-1", "WARS", "WDR26", "XAF1", "YPEL2", 
  "ZBTB7A", "ZCCHC6", "ZFAS1", "ZFC3H1", "ZFP36", "ZFP36L1", "ZFP36L2", 
  "ZFR"
)


# %%
hallmarks$VAN_GALEN_AML_MONO <- AML_VanGalen_monocytes

# %% [markdown]
# #### import genes from Alshetaiwi 2020
# already filtered for abs(logFC)>1 and adjusted P value < 0.05

# %%
Alsh_GMDCs <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/references/Alshetaiwi2020/20250204_aay6017_table_s3_GMDSCs_versus_neutrophils_filtered_human.csv')
Alsh_MMDCs <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/references/Alshetaiwi2020/20250204_aay6017_table_s4_MMDSCs_versus_monocytes_filtered_human.csv')
Alsh_combo_MDCs<- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/references/Alshetaiwi2020/20250204_aay6017_table_s5_Combined_MDSC_signaturegenelist_filtered_human.csv')

# %%
Alsh_GMDCs %>% filter(avg_logFC>1, p_val_adj<0.05) %>% pull(external_gene_name) %>% unique() %>% length() #removed duplicates (they are due to ortholog matching)

# %%
Alsh_GMDCs %>% filter(avg_logFC<(-1), p_val_adj<0.05) %>% unique() %>% dim() # all unique

# %%
hallmarks$ALSH_GMDCs_UP <- Alsh_GMDCs %>% filter(avg_logFC>1, p_val_adj<0.05) %>% unique() %>% pull(external_gene_name)

# %%
hallmarks$ALSH_GMDCs_DOWN <- Alsh_GMDCs %>% filter(avg_logFC<(-1), p_val_adj<0.05) %>% unique() %>% pull(external_gene_name)

# %%
length(hallmarks)

# %% [markdown]
# #### import genes from Kirchberger 2024

# %%
Kirch <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/ref/Kirchberger2024.csv')

# %%
Kirch %>% head()

# %%
M1<-Kirch %>% filter(module_pan=='M1_pan') %>% dplyr::select(gene_symbol, module_pan)

# %%
M2<-Kirch %>% filter(module_pan=='M2_pan') %>% dplyr::select(gene_symbol, module_pan)

# %%
M3<-Kirch %>% filter(module_pan=='M3_pan') %>% dplyr::select(gene_symbol, module_pan)

# %% [markdown]
# ##### zebrafish to human with biomaRt
# useful link: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html#selecting-an-ensembl-biomart-database-and-dataset

# %%
ensembl <- useEnsembl(biomart = "genes")

# %%
human<- useMart("ensembl", dataset = "hsapiens_gene_ensembl") #using hg38

# %%
listAttributes(human) %>% filter(page == 'homologs')

# %%
listAttributes(human) %>% 
  filter(stringr::str_detect(name, "drerio_homolog_"))

# %%
listFilters(human)%>% 
  filter(stringr::str_detect(name, "drerio"))

# %%
attributes<- c('external_gene_name', 'drerio_homolog_associated_gene_name', "drerio_homolog_orthology_type")

# %%
is.atomic(attributes) 

# %%
orth.fish<-  getBM(attributes, filters="with_drerio_homolog",
                    values=TRUE, mart = human, uniqueRows=TRUE)

# %%
colnames(orth.fish)[2] <- "gene_symbol"

# %% [markdown]
# ###### M1

# %%
M1 <- left_join(M1, orth.fish, by= 'gene_symbol')

# %%
M1[is.na(M1$external_gene_name),]

# %%
#manually adding the gene annotations that have no ortholog for no believable reason (from zfin.org)
M1 <- M1 %>%
  mutate(external_gene_name = case_when(
    gene_symbol == "eif3ea"   ~ "EIF3E",
    gene_symbol == "zgc:171772" ~ "RPL37A", 
    gene_symbol == "mki67" ~ "MKI67",
    gene_symbol == "si:ch73-308l14.2" ~ "TCOF1",
    gene_symbol == "si:rp71-45k5.4" ~ "PSMA2",
    gene_symbol == "tardbpl" ~ "TARDBP",
    gene_symbol == "insig1" ~ "INSIG1",
    TRUE ~ external_gene_name 
  ))

# %%
M1_list<-unique(M1$external_gene_name)

# %% [markdown]
# ###### M2

# %%
M2 <- left_join(M2, orth.fish, by= 'gene_symbol')

# %%
M2[is.na(M2$external_gene_name),]

# %%
#manually adding the gene annotations that have no ortholog for no believable reason (from zfin.org)
M2 <- M2 %>%
  mutate(external_gene_name = case_when(
    gene_symbol == "mpc2"   ~ "MPC2",
    gene_symbol == "syne1b" ~ "SYNE1", 
    gene_symbol == "flna" ~ "FLNA",
    gene_symbol == "emc7" ~ "EMC7",
    gene_symbol == "actr1" ~ "ACTR1B",
    gene_symbol == "si:ch211-69g19.2" ~ "CDCA3",
    gene_symbol == "zgc:56525" ~ "GOLM1",
    gene_symbol == "npc2" ~ "NPC2",
    gene_symbol == "mkrn1" ~ "MKRN1",
    TRUE ~ external_gene_name 
  ))

# %%
M2_list<-unique(M2$external_gene_name)

# %% [markdown]
# ###### M3

# %%
M3 <- left_join(M3, orth.fish, by= 'gene_symbol')

# %%
M3[is.na(M3$external_gene_name),]

# %%
#manually adding the gene annotations that have no ortholog for no believable reason (from zfin.org)
M3 <- M3 %>%
  mutate(external_gene_name = case_when(
    gene_symbol == "slc9a3r1a"   ~ "NHERF1",
    gene_symbol == "txn" ~ "TXN", 
    gene_symbol == "zgc:77650" ~ "ARF",
    TRUE ~ external_gene_name 
  ))

# %%
M3_list<-unique(M3$external_gene_name)

# %%
hallmarks$KIRK_M1 <- M1_list
hallmarks$KIRK_M2 <- M2_list
hallmarks$KIRK_M3 <- M3_list

# %% [markdown]
# #### import genes from Montaldo et al 2022

# %%
Montado <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/ref/Montado2022_suppltab29.csv')

# %%
hallmarks$MONTADO_MATURE_BM <- Montado %>% filter(Neutrophil.maturation.stage=='Mature BM neutrophils')  %>% pull(Genes)

# %% [markdown]
# ## DGE

# %%
cells <- unique(meta$CELL_TYPE)
treatments <- unique(meta$TREATMENT)
model <- 'model1'
v <- 'v2'

# %%
for (cell in cells){
    for (treatment in treatments){
        sel_meta <- meta %>% filter(CELL_TYPE == cell, TREATMENT== treatment)
        sel_id<- sel_meta %>% pull(SAMPLE_ID)
        counts.mtx <-v2[,sel_id]
        
        y<-DGEList(counts.mtx)

        y$samples$group<- as_factor(sel_meta$GENOTYPE) #mouse is redundant here

        print("Dimensions before subsetting:")
        print(dim(y))
        keep <- filterByExpr(y)
        y <- y[keep, , keep.lib.sizes=FALSE]
        print("Dimensions after subsetting:")
        print(dim(y))
        
        y <- calcNormFactors(y)

        design <- model.matrix(~ 0 + group, data = y$samples)
        
        y <- estimateDisp(y, design = design)
        fit <- glmQLFit(y, design)

        plotMDS(y, col=ifelse(y$samples$group == "R882H", "red", "blue"))
        title(paste0(cell, '_', treatment))
        plotBCV(y)
        title(paste0(cell, '_', treatment))
        
        setcontrast='groupR882H-groupWT'
        myContrast <- makeContrasts(setcontrast, levels = y$design)
        qlf <- glmQLFTest(fit, contrast=myContrast)

        results <- as.data.frame(topTags(qlf,sort.by = "PValue", adjust.method="BH", n=dim(qlf$table)[1]))
        write.table(results, file = paste0(outdir,"/", prefix, '_', v, '_', model, '_', cell, "_", treatment, "_", setcontrast, '_DGE_results.txt'), 
                            row.names = TRUE, quote = FALSE, sep='\t')
        }
    }

# %% [markdown]
# ##### Import output from code above

# %%
model_date <- '20250413'
model<- 'model1'

# %%
cells <- c('NEU', 'PRE')
treatments <- c('CTRL', 'IFN')
contrast<-'groupR882H-groupWT'

# %%
for (cell in cells) {
    for (treatment in treatments) {
        data <- read.delim(paste0(datadir, '/tables/01.03/', model_date, '_v2_', model, "_", cell, "_", treatment, "_", contrast, '_DGE_results.txt'),
                           header = TRUE, stringsAsFactors = FALSE)
        data$treatment <- treatment
        data$cell <- cell
        data$model<- model
        data$coeff<- contrast
        data<-rownames_to_column(data, var = 'gene')
        assign(paste(model, cell, treatment, sep='_'), data)
        }
    }

# %% [markdown]
# ## GSEA

# %%
for (cell in cells) {
    for (treatment in treatments) {
        title<-paste(model, cell, treatment, sep='_')
        data<-get(title)
        data$fcsign <- sign(data$logFC)
        data$logP <- -log10(data$PValue)
        data$metric <- data$logP/data$fcsign
        data <- data[order(data$metric, decreasing= TRUE), ]
        gene_vector <- data$metric
        names(gene_vector) <- data$gene
        assign(paste('table', 'ranked', model, cell, treatment, sep='_'), data) #generated the table just to check that it was all ok but then the vector is enough
        assign(paste('ranked', model, cell, treatment, sep='_'), gene_vector)
        }
    }

# %%
gsea_hall<-list()

for (cell in cells) {
    for (treatment in treatments) {
        title<-paste('ranked', model, cell, treatment, sep='_')
        ranked_list<-get(title)
        
        fgseaRes.hallmarks <- (fgsea(pathways=hallmarks, stats=ranked_list, minSize=15, maxSize=400, nperm=10000) %>% as.data.frame)[1:5]
        fgseaRes.hallmarks.to.plot <- fgseaRes.hallmarks[fgseaRes.hallmarks$padj<0.05,]
        fgseaRes.hallmarks.to.plot$pathway <- gsub("HALLMARK_","", fgseaRes.hallmarks.to.plot$pathway)
        fgseaRes.hallmarks.to.plot$cell <- cell
        fgseaRes.hallmarks.to.plot$treatment <- treatment
        group<- paste(model, cell, treatment, sep='_')
        
        assign(paste('gsea_hall_', model, cell, treatment, sep='_'), fgseaRes.hallmarks.to.plot)
        gsea_hall[[group]] <- fgseaRes.hallmarks.to.plot
      }
}    

gsea_hallmarks <- bind_rows(gsea_hall, .id="group")


# %% [markdown]
# ### bubble plot

# %%
gsea_hallmarks$cell <- factor(gsea_hallmarks$cell, levels=c('PRE', 'NEU'))

# %%
fig(8.5,11)
p0<- ggplot(gsea_hallmarks, aes(cell, pathway, size=-log10(padj), col=NES)) + 
    geom_point() + 
    scale_color_gradient2(low='deepskyblue4', mid = 'gray90', high='firebrick3', name='NES') +
    xlab('') + ylab('') +
    geom_vline(xintercept = seq(1, 1, 1) + 0.5, alpha = 0.2)+
    #geom_vline(xintercept = seq(1, 4, 1) + 0.5, color = "grey20", alpha = 1, linewidth=0.5)+
    ggtitle(paste0('DNMT3A MUT vs WT') ) + 
    facet_wrap(~treatment, nrow=1)+
    theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=0.5),
         axis.text=element_text(size=15, color='black'),
         text=element_text(size=15),
            panel.border = element_rect(color='black', linewidth=1, fill = NA),
            panel.background = element_blank(),
            axis.line = element_blank(),
          legend.position = 'right',
          strip.background = element_blank(),
          strip.text = element_text(size = 15, face = "bold", color = "black"),
          legend.title = element_text(vjust=1),
        plot.margin=margin(1, 3, 1, 1, "lines")
)
p0

# %%
ggsave(paste0(figdir,prefix,'_model1_GSVEA_allsignificant.pdf'), p0, w=8.5, h=11)

# %%
sel_path<- c('E2F_TARGETS','MITOTIC_SPINDLE','G2M_CHECKPOINT','PI3K_AKT_MTOR_SIGNALING',
                       'KRAS_SIGNALING_UP','APOPTOSIS', 'P53_PATHWAY', 
                       'INFLAMMATORY_RESPONSE', 'HACKERT_NEUTS_INFLAMMATION_UP',
                       'TNFA_SIGNALING_VIA_NFKB',  'INTERFERON_ALPHA_RESPONSE', 
                       'INTERFERON_GAMMA_RESPONSE', 'IL6_JAK_STAT3_SIGNALING', 'TGF_BETA_SIGNALING'
                        ) %>% rev()

# %%
fig(8,7)
p1 <- ggplot(
  gsea_hallmarks %>%
    filter(pathway %in% sel_path) %>%
    mutate(pathway = factor(pathway, levels = sel_path)),
  aes(cell, pathway, size = -log10(padj), col = NES)
) +
  geom_point() +
  scale_color_gradient2(low = 'deepskyblue4', mid = 'gray90', high = 'firebrick3', name = 'NES') +
  xlab('') + ylab('') +
  geom_vline(xintercept = seq(1, 1, 1) + 0.5, alpha = 0.2) +
  ggtitle(paste0('DNMT3A MUT vs WT')) +
  facet_wrap(~treatment, nrow = 1) +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
    axis.text = element_text(size = 15, color = 'black'),
    text = element_text(size = 15),
    panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
    panel.background = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 15, face = "bold", color = "black"),
    legend.title = element_text(vjust = 1),
    plot.margin = margin(1, 3, 1, 1, "lines")
  )
p1


# %%
ggsave(paste0(figdir,prefix,'_model1_GSVEA_14significant_selected.pdf'), p1, w=8, h=7.5)

# %% [markdown]
# ### GSEA enrichment plot Neut cluster

# %%
gsea_plot<- fgsea(pathways=hallmarks, stats=ranked_model1_NEU_CTRL, minSize=15, maxSize=400, nperm=10000) #%>% as.data.frame

# %%
fig(8,8)
plot_pathway<-gsea_plot$pathway[54]

pd <- plotEnrichmentData(
    pathway = hallmarks[[plot_pathway]],
    stats = ranked_model1_NEU_CTRL
)

pval <- gsea_plot %>% 
  filter(pathway == plot_pathway) %>% 
  pull(padj) %>% round(6)
NES<- gsea_plot %>% 
  filter(pathway == plot_pathway) %>% 
  pull(NES) %>% round(3)

p2<- with(pd,
     ggplot(data=curve) +
         geom_line(aes(x=rank, y=ES), color="chartreuse3", size=2) + #chartreuse2, size1.5
         geom_ribbon(data=stats,
                     mapping=aes(x=rank, ymin=0,
                                 ymax=stat/maxAbsStat*(spreadES/4)),
                     fill="grey") +
         geom_segment(data=ticks,
                      mapping=aes(x=rank, y=-spreadES/16,
                                  xend=rank, yend=spreadES/16),
                      size=0.4) +
         geom_hline(yintercept=posES, colour="red", linetype="dashed", size=1.2) + #size1
         geom_hline(yintercept=negES, colour="red", linetype="dashed", size=1.2) + #size1
         geom_hline(yintercept=0, colour="black") +
         annotate("text", x=8000, y=0.34, label=paste0('padj=',pval), size=8, color="gray26") + #size 6
         annotate("text", x=8000, y=0.30, label=paste0('NES=',NES), size=8, color="gray26") + #size 6
         theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
         ) +
         theme_minimal(base_size = 24)+
         labs(x="rank", y="enrichment score", title=plot_pathway, subtitle='CTR_NEU'))
p2

# %%
ggsave(paste0(figdir,prefix,'_model1_GSEA_ALSH_GMDCs_UP.pdf'), p2, w=7, h=7)

# %%
fig(8,8)
plot_pathway<-gsea_plot_CTR$pathway[gsea_plot_CTR$pathway=='KIRK_M3']

pd <- plotEnrichmentData(
    pathway = hallmarks[[plot_pathway]],
    stats = ranked_model1_NEU_CTRL
)

pval <- gsea_plot_CTR %>% 
  filter(pathway == plot_pathway) %>% 
  pull(padj) %>% round(6)
NES<- gsea_plot_CTR %>% 
  filter(pathway == plot_pathway) %>% 
  pull(NES) %>% round(3)

p2<- with(pd,
     ggplot(data=curve) +
         geom_line(aes(x=rank, y=ES), color="chartreuse3", size=2) + #chartreuse2, size1.5
         geom_ribbon(data=stats,
                     mapping=aes(x=rank, ymin=0,
                                 ymax=stat/maxAbsStat*(spreadES/4)),
                     fill="grey") +
         geom_segment(data=ticks,
                      mapping=aes(x=rank, y=-spreadES/16,
                                  xend=rank, yend=spreadES/16),
                      size=0.4) +
         geom_hline(yintercept=posES, colour="red", linetype="dashed", size=1.2) + #size1
         geom_hline(yintercept=negES, colour="red", linetype="dashed", size=1.2) + #size1
         geom_hline(yintercept=0, colour="black") +
         annotate("text", x=6000, y=0.55, label=paste0('padj=',pval), size=8, color="gray26") + #size 6
         annotate("text", x=6000, y=0.50, label=paste0('NES=',NES), size=8, color="gray26") + #size 6
         theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
         ) +
         theme_minimal(base_size = 24)+
         labs(x="rank", y="enrichment score", title=plot_pathway, subtitle='CTR_NEU'))
p2

# %%
ggsave(paste0(figdir,prefix,'_model1_GSEA_Kirchberger_M3.pdf'), p2, w=12,h=7)

# %%
fig(8,8)
plot_pathway<-gsea_plot_CTR$pathway[gsea_plot_CTR$pathway=='MONTADO_MATURE_BM']

pd <- plotEnrichmentData(
    pathway = hallmarks[[plot_pathway]],
    stats = ranked_model1_NEU_CTRL
)

pval <- gsea_plot_CTR %>% 
  filter(pathway == plot_pathway) %>% 
  pull(padj) %>% round(6)
NES<- gsea_plot_CTR %>% 
  filter(pathway == plot_pathway) %>% 
  pull(NES) %>% round(3)

p2<- with(pd,
     ggplot(data=curve) +
         geom_line(aes(x=rank, y=ES), color="chartreuse3", size=2) + #chartreuse2, size1.5
         geom_ribbon(data=stats,
                     mapping=aes(x=rank, ymin=0,
                                 ymax=stat/maxAbsStat*(spreadES/4)),
                     fill="grey") +
         geom_segment(data=ticks,
                      mapping=aes(x=rank, y=-spreadES/16,
                                  xend=rank, yend=spreadES/16),
                      size=0.4) +
         geom_hline(yintercept=posES, colour="red", linetype="dashed", size=1.2) + #size1
         geom_hline(yintercept=negES, colour="red", linetype="dashed", size=1.2) + #size1
         geom_hline(yintercept=0, colour="black") +
         annotate("text", x=6000, y=0.55, label=paste0('padj=',pval), size=8, color="gray26") + #size 6
         annotate("text", x=6000, y=0.50, label=paste0('NES=',NES), size=8, color="gray26") + #size 6
         theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
         ) +
         theme_minimal(base_size = 24)+
         labs(x="rank", y="enrichment score", title=plot_pathway, subtitle='CTR_NEU'))
p2

# %%
ggsave(paste0(figdir,prefix,'_model1_GSVEA_Montado_MatureBM.pdf'), p2, w=12,h=7)
