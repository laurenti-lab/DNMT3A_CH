# %% [markdown]
# # 03.02_DGE with EdgeR and GSEA
# Author:GM <br>
# useful links: https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf <br>
# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html and https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html <br>
# https://www.sc-best-practices.org/conditions/differential_gene_expression.html
#

# %%
suppressPackageStartupMessages({library(tidyverse)
                                library(gridExtra)
                                library(RColorBrewer)
                                library(edgeR)
                                library(MAST)
                                library(data.table)
                                library(fgsea)
                                library(anndata)})

# %%
cat("R version:", R.version$version.string, "\n")

installed.packages() %>%
  as_tibble() %>%
  select(Package, Version) %>%
  filter(Package %in% c('tidyverse','RColorBrewer', 'edgeR', 'MAST', 'anndata')) %>%
  print(row.names = FALSE)

# %%
basedir <- paste(laurenti, user, project, sep='/')
basename <- '_CHIP27_'
datadir<-paste(basedir, 'output', sep='/')

# %% [markdown]
# ## define functions

# %% [markdown]
# ### 1. Fit model

# %% [markdown]
# #### fit_model1: MUT vs WT untreated: design = ~0 + mut + sort
# group=mut

# %%
fit_model1 <- function(adata_){
    # create an edgeR object with counts and grouping factor
    y <- DGEList(counts=as.matrix(t(adata_$X)), group = adata_$obs$genotype)
    # filter out genes with low counts 
    print("Dimensions before subsetting:")
    print(dim(y))
    print("")
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE] #recalculate library size post filtering
    print("Dimensions after subsetting:")
    print(dim(y))
    print("")
    # normalize
    y <- calcNormFactors(y)
    
    group <- adata_$obs$genotype
    sort <- adata_$obs$FACS
    
    design <- model.matrix(~ 0 + group + sort)
    # estimate dispersion : based on the model derived from the glm model that we choose
    y <- estimateDisp(y, design = design)
    # fit the model
    fit <- glmQLFit(y, design)

    #plot to explore Explore sample relations based on multidimensional scaling and dispersion
    plotMDS(y, col=ifelse(y$samples$group == "MUT", "red", "blue"))
    title(cluster)
    plotBCV(y)
    title(cluster)
    
    return(list("fit"=fit, "design"=design, "y"=y))
}

# %% [markdown]
# #### fit_model2: MUT vs WT treated: design = ~0 + mut
# Same model as above but sort exlcuded because redundant: all treated samples have been FACS enriched

# %%
fit_model2 <- function(adata_){
    # create an edgeR object with counts and grouping factor
    y <- DGEList(counts=as.matrix(t(adata_$X)), group = adata_$obs$genotype)
    # filter out genes with low counts -- performed with a function and not manually
    print("Dimensions before subsetting:")
    print(dim(y))
    print("")
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE] #recalculate library size post filtering
    print("Dimensions after subsetting:")
    print(dim(y))
    print("")
    # normalize
    y <- calcNormFactors(y)
    
    group <- adata_$obs$genotype
    #sort <- adata_$obs$FACS
    
    design <- model.matrix(~ 0 + group)
    # estimate dispersion : based on the model derived from the glm model that we choose
    y <- estimateDisp(y, design = design)
    # fit the model
    fit <- glmQLFit(y, design)

    #plot to explore Explore sample relations based on multidimensional scaling and dispersion
    plotMDS(y, col=ifelse(y$samples$group == "MUT", "red", "blue"))
    title(cluster)
    plotBCV(y)
    title(cluster)
    
    return(list("fit"=fit, "design"=design, "y"=y))
}

# %% [markdown]
# #### fit_model3: MUT:IFN interaction : design = ~0 + mut + sort + treatment + mut:treatment
# to evaluate the effect of interferon treatment independently of mutational status

# %%
fit_model3 <- function(adata_){
    # create an edgeR object with counts and grouping factor
    
    y <- DGEList(counts=as.matrix(t(adata_$X)), group = adata_$obs$genotype)
    
    # filter out genes with low counts -- performed with a function and not manually
    print("Dimensions before subsetting:")
    print(dim(y))
    print("")
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE] #recalculate library size post filtering
    print("Dimensions after subsetting:")
    print(dim(y))
    print("")
    # normalize
    y <- calcNormFactors(y)
    
    group <- adata_$obs$genotype
    sort <- adata_$obs$FACS
    treatment <-adata_$obs$treatment

    
    design <- model.matrix(~ 0 + group + sort + treatment + group:treatment)
    
    # estimate dispersion : based on the model derived from the glm model that we choose
    y <- estimateDisp(y, design = design)
    # fit the model
    fit <- glmQLFit(y, design)

    #plot to explore Explore sample relations based on multidimensional scaling and dispersion
    plotMDS(y, col=ifelse(y$samples$group == "MUT", "red", "blue"))
    title(cluster)
    plotBCV(y)
    title(cluster)
    
    return(list("fit"=fit, "design"=design, "y"=y))
}

# %% [markdown]
# ## 2. calculate contrasts and ranked table (with adjusted fdr)

# %%
# VERY IMPORTANT you need to setcontrast before running this
#eg. setcontrast <-'groupMUT-groupWT'

#condition='untreated'only important to name the files


contrast <- function(outs, outputname, condition){ #where outs is the output of the fit_model function

                fit <- outs$fit
                y <- outs$y
                
                myContrast <- makeContrasts(setcontrast, levels = y$design)
                
                qlf <- glmQLFTest(fit, contrast=myContrast)
                
                group_res <- topTags(qlf, sort.by = "PValue", adjust.method="BH", n=dim(qlf$table)[1])
                
                print(head(group_res))
                    
                    saveRDS(qlf, paste0(datadir,"/objects/edgeR/", prefix, "_", outputname, "_", condition, "_", setcontrast, "_edgeR.rds"))
                
                    write.table(group_res, 
                            file = paste0(datadir,"/tables/edgeR/", prefix, "_", outputname, "_", condition, "_", setcontrast, '_DGE_results.txt'), 
                            row.names = TRUE, quote = FALSE, sep='\t')
                    
                    #for GSVA
                    group_res$table$FDR[group_res$table$FDR == 0.000000e+00 ] <- 1e-323
                    group_res$table$PValue[group_res$table$PValue == 0.000000e+00 ] <- 1e-323
                    
                    er <- group_res$table
                    er$genes <- rownames(er)
                
                    er$fcsign <- sign(er$logFC)
                    er$logP=-log10(er$PValue)
                    er$metric= er$logP/er$fcsign
                
                
                    final<-er[,c("genes", "metric")]
                    
                    write.table(na.exclude( final[order(final$metric, decreasing = TRUE), ] ), 
                            file = paste0(datadir,"/tables/edgeR/", prefix, "_", outputname, "_", condition, "_", setcontrast,'.rnk'), 
                            row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')

                    plotSmear(qlf, de.tags = rownames(group_res)[which(group_res$FDR<0.01)])
                    abline(h = c(-1.5, 1.5), col = "blue")
}

# %% [markdown]
# #### model 3 contrasts

# %%
#adata_$obs$genotype <-relevel(factor(adata_$obs$genotype), ref = "MUT") #essential to get groupMUT-groupWT in the correct order
contrast_model3 <- function(outs, outputname, condition){ #where outs is the output of the fit_model function

                fit <- outs$fit
                y <- outs$y
                colnames(y$design)<- c('groupMUT', 'groupWT', 'sortenriched', 'treatmentIFNy','groupWT_treatmentIFNy')
                
                myContrast <- makeContrasts(setcontrast, levels = y$design)
                
                qlf <- glmQLFTest(fit, contrast=myContrast)
                
                group_res <- topTags(qlf, sort.by = "PValue", adjust.method="BH", n=dim(qlf$table)[1])
                
                print(head(group_res))
                    
                    saveRDS(qlf, paste0(datadir,"/objects/edgeR/", prefix, "_", outputname, "_", condition, "_", setcontrast, "_edgeR.rds"))
                
                    write.table(group_res, 
                            file = paste0(datadir,"/tables/edgeR/", prefix, "_", outputname, "_", condition, "_", setcontrast, '_DGE_results.txt'), 
                            row.names = TRUE, quote = FALSE, sep='\t')
                    
                    #for GSVA
                    group_res$table$FDR[group_res$table$FDR == 0.000000e+00 ] <- 1e-323
                    group_res$table$PValue[group_res$table$PValue == 0.000000e+00 ] <- 1e-323
                    
                    er <- group_res$table
                    er$genes <- rownames(er)
                
                    er$fcsign <- sign(er$logFC)
                    er$logP=-log10(er$PValue)
                    er$metric= er$logP/er$fcsign
                
                
                    final<-er[,c("genes", "metric")]
                    
                    write.table(na.exclude( final[order(final$metric, decreasing = TRUE), ] ), 
                            file = paste0(datadir,"/tables/edgeR/", prefix, "_", outputname, "_", condition, "_", setcontrast,'.rnk'), 
                            row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')

                    plotSmear(qlf, de.tags = rownames(group_res)[which(group_res$FDR<0.01)])
                    abline(h = c(-1.5, 1.5), col = "blue")
}

# %%
#condition <- "interaction"
#adata_$obs$genotype <-relevel(factor(adata_$obs$genotype), ref = "WT") #essential to get MUT:IFNy interaction
#coefficient<-"groupMUT:treatmentIFNy" OR "treatmentIFNy"

contrast_model3_coefficient <- function(outs, outputname, condition){ #where outs is the output of the fit_model function

                fit <- outs$fit
                y <- outs$y
                
                qlf <- glmQLFTest(fit, coef = coefficient)
                
                group_res <- topTags(qlf, sort.by = "PValue", adjust.method="BH", n=dim(qlf$table)[1])
                
                print(head(group_res))

                    coefficient <- if (coefficient == "groupMUT:treatmentIFNy") {
                        gsub(":", "-", coefficient)
                    } else {
                        coefficient
                    }
                                                            
                    saveRDS(qlf, paste0(datadir,"/objects/edgeR/", prefix, "_", outputname, "_all_", coefficient, "_edgeR.rds"))
                
                    write.table(group_res, 
                            file = paste0(datadir,"/tables/edgeR/", prefix, "_", outputname, "_all_", coefficient, "_DGE_results.txt"), 
                            row.names = TRUE, quote = FALSE, sep='\t')
                    
                    #for GSVA
                    group_res$table$FDR[group_res$table$FDR == 0.000000e+00 ] <- 1e-323
                    group_res$table$PValue[group_res$table$PValue == 0.000000e+00 ] <- 1e-323
                    
                    er <- group_res$table
                    er$genes <- rownames(er)
                
                    er$fcsign <- sign(er$logFC)
                    er$logP=-log10(er$PValue)
                    er$metric= er$logP/er$fcsign
                
                
                    final<-er[,c("genes", "metric")]
                    
                    write.table(na.exclude( final[order(final$metric, decreasing = TRUE), ] ), 
                            file = paste0(datadir,"/tables/edgeR/", prefix, "_", outputname, "_all_", coefficient, ".rnk"), 
                            row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')

                    plotSmear(qlf, de.tags = rownames(group_res)[which(group_res$FDR<0.01)])
                    abline(h = c(-1.5, 1.5), col = "blue")
}

# %% [markdown]
# ## 3. visualisation fxs

# %%
options(repr.matrix.max.cols=Inf)

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %%
DEGs_top_plot <- function(data, top_n_genes){
            top_genes <- data %>% 
                         slice_head(n = top_n_genes) %>%
                         pull(gene)
            
            #this mtx is only to get hclust of the genes
            DEGs.mtx <- data %>%
                          filter(gene %in% top_genes) %>%
                          select(gene, cluster, treatment, logFC) %>%
                          pivot_wider(names_from=gene, values_from=logFC, values_fill=0) %>%
                          column_to_rownames(var = "group")
            
            hc_order = hclust(dist(t(DEGs.mtx %>% select(-cluster,-treatment))))
            
            DEGs_plot <- data %>% 
                            filter(gene %in% top_genes) %>%
                            mutate(gene=factor(gene, levels=hc_order$labels[hc_order$order]))
            assign(paste0("DEGs_plot"), DEGs_plot, envir = .GlobalEnv)
            }

# %%
bub_plot <- function(data, title) {
  p <- ggplot(data, aes(cluster, gene, size = -log10(FDR), col = logFC_capped)) +
    geom_point() +
    scale_color_gradient2(
      low = 'deepskyblue4', mid = 'gray90', high = 'firebrick3', 
      name = 'logFC_capped'
    ) +
    xlab('') + 
    ylab('') +
    geom_vline(xintercept = seq(1, 4, 1) + 0.5, alpha = 0.2) +
    ggtitle(title) +
    facet_wrap(~treatment, nrow = 1) +
    theme(
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
      axis.text = element_text(size = 15, color = 'black'),
      text = element_text(size = 15),
      panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
      panel.background = element_blank(),
      axis.line = element_blank(),
      legend.position = 'bottom',
      strip.background = element_blank(),
      strip.text = element_text(size = 15, face = "bold", color = "black"),
      legend.title = element_text(vjust = 1),
      plot.margin = margin(1, 3, 1, 1, "lines")
    )
  clean_title <- gsub(" ", "_", title)
  assign(paste0("p_", clean_title), p, envir = .GlobalEnv)
}


# %% [markdown]
# ## run EdgeR

# %% [markdown]
# ### import anndata

# %%
anndata_pb <- anndata::read_h5ad(paste0(basedir, '/output/objects/', '20241120_CHIP27_pseudobulks_cluster_1120.h5ad'))

# %%
anndata_pb$obs  <- anndata_pb$obs %>%
                  rownames_to_column(var = "row_index") %>%
                  mutate(row_index = str_split(row_index, "-", n = 2) %>% 
                                           sapply(pluck, 1),
                         row_index= str_split(row_index, "_", n=2) %>%
                                          sapply(pluck,2),
                     row_index = paste0('cluster',cluster_1120,'_', row_index)) %>%
                  column_to_rownames(var = "row_index")

anndata_pb$obs  %>% head()

# %%
clusters <- anndata_pb$obs$cluster_1120	%>% unique

# %% [markdown]
# ### Model1: MUT vs WT untreated: design = ~0 + mut + sort
#

# %%
adata_CTR <- anndata_pb[anndata_pb$obs["treatment"] == "CTR", ]

# %% jupyter={"outputs_hidden": true}
setcontrast <- 'groupMUT-groupWT'

for (cluster in clusters) {
    print('untreated dataset only') 
    print(paste0('results for cluster ',cluster, ' cluster_1120'))
    adata_subset <- adata_CTR[adata_CTR$obs["cluster_1120"] == cluster, ]
    outs <- fit_model1(adata_subset)
    contrast(outs, outputname=paste0(as.character(cluster), '_model1'), condition='untreated')
    }

# %% [markdown]
# ### Model2: MUT vs WT treated: design = ~0 + mut

# %%
adata_TX <- anndata_pb[anndata_pb$obs["treatment"] == "IFNy", ]

# %%
clusters_no8<- clusters[clusters != "8"] #cluster 8removed because too few cells

# %% jupyter={"outputs_hidden": true}
setcontrast <- 'groupMUT-groupWT'

for (cluster in clusters_no8) {
    print('treated dataset only') 
    print(paste0('results for cluster ',cluster, ' cluster_1120'))
    adata_subset <- adata_TX[adata_TX$obs["cluster_1120"] == cluster, ]
    outs <- fit_model2(adata_subset)
    contrast(outs, outputname=paste0(as.character(cluster), '_model2'), condition='treated')
    }

# %% [markdown]
# ### Model3: MUT:IFN interaction : design = ~0 + mut + sort + treatment + mut:treatment

# %%
adata<-anndata_pb

# %% jupyter={"outputs_hidden": true}
adata$obs$genotype <-relevel(factor(adata$obs$genotype), ref = "MUT")
#the order of the column is essential, I want the same as above (MUT-WT) to get comparable results

for (cluster in clusters_no8) {
    print('all dataset') 
    print(paste0('results for cluster ',cluster, ' cluster_1120'))
    adata_subset <- adata[adata$obs["cluster_1120"] == cluster, ]
    outs <- fit_model3(adata_subset)
    setcontrast <- 'groupMUT-groupWT'
    contrast_model3(outs, outputname=paste0(as.character(cluster), '_model3'), condition='all')
    }

# %% jupyter={"outputs_hidden": true}
#here I change the order of the columns of my design matrix to obtain MUT:IFNy interaction values
adata$obs$genotype <-relevel(factor(adata$obs$genotype), ref = "WT")
#coefficient<-"groupMUT:treatmentIFNy" OR "treatmentIFNy"

for (cluster in clusters_no8) {
    print('all dataset') 
    print(paste0('results for cluster ',cluster, ' cluster_1120'))
    adata_subset <- adata[adata$obs["cluster_1120"] == cluster, ]
    outs <- fit_model3(adata_subset)
    coefficient<-'groupMUT:treatmentIFNy'
    contrast_model3_coefficient(outs, outputname=paste0(as.character(cluster), '_model3'))
    coefficient<-'treatmentIFNy'
    contrast_model3_coefficient(outs, outputname=paste0(as.character(cluster), '_model3'))
    }

# %% [markdown]
# ## DEGs visualisation

# %% [markdown]
# ### DEGs table 
# Import DEGs table output from code above

# %%
model_date <- '20241127'

# %%
clusters <- c('0', '1_7', '2_5', '3_4', '6', '8')

# %% [markdown]
# ##### import DGE list model1: untreated only; design = ~0 + mut + sort ; contrast= 1'*groupMUT -1'*groupWT 

# %%
treatment <- 'untreated'
model <- 'model1'
contrast<-'groupMUT-groupWT'

# %%
for (cluster in clusters) {
    data <- read.delim(paste0(datadir, '/tables/edgeR/', model_date, '_', cluster, '_',
                              model, '_', treatment, '_', contrast, '_DGE_results.txt'), header = TRUE, stringsAsFactors = FALSE)
    data$treatment <- treatment
    data$cluster <-cluster
    data$model<- model
    data$coeff<- contrast
    data<-rownames_to_column(data, var = 'gene')
    assign(paste0(treatment, '_cl', cluster), data)
  }


# %% [markdown]
# ##### import DGE list model2: treated only; design = ~0 + mut ; contrast= 1'*groupMUT -1'*groupWT 

# %%
clusters_no8<- clusters[clusters != "8"] #cl8 only present in model1 
treatment <- 'treated'
model <- 'model2'
contrast<-'groupMUT-groupWT'

# %%
for (cluster in clusters_no8) {
    data <- read.delim(paste0(datadir, '/tables/edgeR/', model_date, '_', cluster, '_',
                              model, '_', treatment, '_', contrast, '_DGE_results.txt'), header = TRUE, stringsAsFactors = FALSE)
    data$treatment <- treatment
    data$cluster <-cluster
    data$model<- model
    data$coeff<- contrast
    data<-rownames_to_column(data, var = 'gene')
    assign(paste0(treatment, '_cl', cluster), data)
  }


# %% [markdown]
# ##### import DGE list model3: all dataset; design=~0 + mut + sort + treatment + mut:treatment; contrast= 1'*groupMUT -1'*groupWT; interaction=groupMUT:treatmentIFNy; coefficient: treatment IFNy

# %%
clusters_no8<- clusters[clusters != "8"]
treatment <- 'all'
model <- 'model3'
contrast<-c('groupMUT-groupWT','groupMUT-treatmentIFNy','treatmentIFNy')

# %%
for (c in contrast) {
    for (cluster in clusters_no8) {
        data <- read.delim(paste0(datadir, '/tables/edgeR/', model_date, '_', cluster, '_',
                                  model, '_', treatment, '_', c, '_DGE_results.txt'), header = TRUE, stringsAsFactors = FALSE)
        data$treatment <- treatment
        data$cluster <-cluster
        data$model<- model
        data$coeff<- c
        data<-rownames_to_column(data, var = 'gene')
        name <- if (c %in% c('groupMUT-groupWT', 'groupMUT-treatmentIFNy')) {
          gsub('-', '_', c)
        } else {
          c
        }
        assign(paste0(treatment,'_', name, '_cl', cluster), data)
      }
    }


# %% [markdown]
# ##### merge data

# %%
prefixes <- c('untreated','treated','all_groupMUT_groupWT','all_groupMUT_treatmentIFNy', 'all_treatmentIFNy')

# %%
all_data <- list()

for (p in prefixes) {
  for (cluster in clusters_no8) {
    data_name <- paste0(p, "_cl", cluster)
      data <- get(data_name)
      data$group <- paste0(data$cluster, '_', data$treatment) 
      data <- data %>%
              mutate(group = case_when(
                coeff == 'treatmentIFNy' ~ paste0(group, 'IFNy'),
                coeff == 'groupMUT-treatmentIFNy' ~ paste0(group, 'MUTIFNy'),
                TRUE ~ group
              ))
      #data <- rownames_to_column(data, var = 'gene')
      all_data[[data_name]] <- data
  }
}

DEGs <- bind_rows(all_data)


# %% [markdown]
# ### Establish significance thresholds

# %%
FDR_thres = 0.05
logFC_thres = 1

# %%
DEGs_sig <- DEGs %>%
            filter(FDR < FDR_thres, abs(logFC) > logFC_thres) %>%
            group_by(group, cluster) %>%
            arrange(desc(abs(logFC)), .by_group=TRUE) %>%
            select(-treatment) %>%
            mutate(treatment= str_split(group, "_(?=[^_]+$)", n = 2)[[1]][2]) %>% 
            mutate(treatment=factor(treatment, levels= c('untreated', 'treated', 'all', 'allMUTIFNy', 'allIFNy')))%>%
            mutate(cluster=factor(cluster, levels=c('2_5','1_7','0','6','3_4'))) %>%
            mutate(logFC_capped = ifelse(logFC >= 5, 5, ifelse(logFC <= -5, -5, logFC)))

# %% [markdown]
# ### DEGs count per cluster

# %%
DEGs_count <- DEGs_sig %>%
  mutate(direction = ifelse(logFC > 0, "UP", "DOWN")) %>%
  group_by(treatment, coeff, cluster, direction, model) %>%
  summarise(count = n(), .groups = 'drop')

# %%
DEGs_wide <- DEGs_count %>%
  pivot_wider(names_from = direction, values_from = count)

# %%
DEGs_wide <- DEGs_wide %>% 
  mutate(cluster = factor(cluster, levels = c('2_5','1_7','0', '6', '3_4')))

# %%
options(repr.plot.width=6, repr.plot.height=2)

p3 <- ggplot(DEGs_wide %>% filter(treatment %in% c('treated', 'untreated')),
             aes(x = factor(cluster), y = treatment, fill = UP, label=UP)) +
  geom_tile() +
  geom_text(color = "grey20", size = 4, fontface = "bold", show.legend = FALSE) + 
  scale_fill_gradient(low = "white", high = "firebrick3", guide="none") +
  labs(title = "UP genes",
       x = "",
       y = "",
       fill = "Count") +
  expand_limits(x = 0.5:length(unique(DEGs_wide$cluster)))+
  theme_classic(base_size=16)

p4 <- ggplot(DEGs_wide %>% filter(treatment %in% c('treated', 'untreated')),
             aes(x = factor(cluster), y = treatment, fill = DOWN, label = DOWN)) +
  geom_tile() +
  geom_text(color = "grey20", size = 4, fontface = "bold", show.legend = FALSE) + 
  scale_fill_gradient(low = "white", high = "dodgerblue4", guide="none") +
  labs(title = "DOWN genes",
       x = "Cluster",
       y = "",
       fill = "Count") +
  theme_classic(base_size=16)

fig(5,5)
p7<- grid.arrange(p3, p4, ncol = 1, top = paste0('FDR < ', FDR_thres, '  abs(logFC)>', logFC_thres, '\nmodel1,2 groupMUT-groupWT'))

# %% [markdown]
# ### top DEGs

# %%
DEGs_sig <- DEGs_sig %>% 
  mutate(cluster = factor(cluster, levels = c('2_5','1_7','0', '6', '3_4')))

# %%
fig(6,30)
DEGs_top_plot(data=DEGs_sig %>% filter(treatment %in% c('untreated', 'treated')), top_n_genes=20)
bub_plot(DEGs_plot, 'DNMT3A MUT vs WT')

# %%
#top genes UP

fig(6,30)

DEGs_top_plot(data=DEGs_sig %>% 
              filter(treatment %in% c('untreated', 'treated'), logFC>0), top_n_genes=20)
bub_plot(DEGs_plot, 'DNMT3A MUT vs WT up')
p_DNMT3A_MUT_vs_WT_up

# %%
#top genes DOWN

fig(6,30)
DEGs_top_plot(data=DEGs_sig %>% 
              filter(treatment %in% c('untreated', 'treated'), logFC<0), top_n_genes=20)
bub_plot(DEGs_plot, 'DNMT3A MUT vs WT up')
p_DNMT3A_MUT_vs_WT_up

# %% [markdown]
# ### selected DEGs

# %%
#transformation genes
genes_sel2 <- c('HOXA9','HOXA10','HOPX','FGF13','HGF')  %>% rev()

fig(8,4.5)
logFC_capping= 3 
log10_FDR_capping = (10)
p1<- ggplot(DEGs_sig %>% 
         filter(gene %in% genes_sel2) %>% filter(treatment %in% c('untreated', 'treated')) %>% 
         mutate(gene = factor(gene, levels = genes_sel2)) %>%
         mutate(logFC_capped = ifelse(logFC >= logFC_capping, logFC_capping, ifelse(logFC <= -logFC_capping, -logFC_capping, logFC))) %>%
         mutate(minus_log10_FDR=-log10(FDR)) %>%
         mutate(capped_minus_log10_FDR = ifelse(minus_log10_FDR >= log10_FDR_capping, log10_FDR_capping, minus_log10_FDR))
       , aes(cluster, gene, size = capped_minus_log10_FDR, col = logFC_capped)) +
    geom_point() +
    scale_color_gradient2(
      low = 'deepskyblue4', mid = 'gray90', high = 'firebrick3', 
      name = 'logFC_capped'
    ) +
    xlab('') + 
    ylab('') +
    geom_vline(xintercept = seq(1, 4, 1) + 0.5, alpha = 0.2) +
    ggtitle('Selected DEGs') +
    facet_wrap(~treatment, nrow = 1) +
    theme(
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
      axis.text = element_text(size = 15, color = 'black'),
      text = element_text(size = 15),
      panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
      panel.background = element_blank(),
      axis.line = element_blank(),
      legend.position = 'bottom',
      legend.box= 'vertical',
      strip.background = element_blank(),
      strip.text = element_text(size = 15, face = "bold", color = "black"),
      legend.title = element_text(vjust = 1),
      plot.margin = margin(1, 10, 1, 1, "lines")
    )
p1

# %%
#inflammation genes
genes_sel2 <- c("DEFA3", "RETN","IL2RA",
                    "SERPINB" , "CXCL10", "SDC4", "CCRL2", "TNF", "IL1B","IL1R1","CCL20","COLEC12","BATF3","PTGER2",
                    "IL1A","IL7R","TLR7","PROCR") %>% rev()

fig(8,8)
logFC_capping= 3 
log10_FDR_capping = (10)
p3<- ggplot(DEGs_sig %>% 
         filter(gene %in% genes_sel2) %>% filter(treatment %in% c('untreated', 'treated')) %>% 
         mutate(gene = factor(gene, levels = genes_sel2)) %>%
         mutate(logFC_capped = ifelse(logFC >= logFC_capping, logFC_capping, ifelse(logFC <= -logFC_capping, -logFC_capping, logFC))) %>%
         mutate(minus_log10_FDR=-log10(FDR)) %>%
         mutate(capped_minus_log10_FDR = ifelse(minus_log10_FDR >= log10_FDR_capping, log10_FDR_capping, minus_log10_FDR))
       , aes(cluster, gene, size = capped_minus_log10_FDR, col = logFC_capped)) +
    geom_point() +
    scale_color_gradient2(
      low = 'deepskyblue4', mid = 'gray90', high = 'firebrick3', 
      name = 'logFC_capped'
    ) +
    xlab('') + 
    ylab('') +
    geom_vline(xintercept = seq(1, 4, 1) + 0.5, alpha = 0.2) +
    ggtitle('Selected DEGs inflammation') +
    facet_wrap(~treatment, nrow = 1) +
    theme(
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
      axis.text = element_text(size = 15, color = 'black'),
      text = element_text(size = 15),
      panel.border = element_rect(color = 'black', linewidth = 1, fill = NA),
      panel.background = element_blank(),
      axis.line = element_blank(),
      legend.position = 'bottom',
      legend.box= 'vertical',
      strip.background = element_blank(),
      strip.text = element_text(size = 15, face = "bold", color = "black"),
      legend.title = element_text(vjust = 1),
      plot.margin = margin(1, 10, 1, 1, "lines")
    )
p3

# %% [markdown]
# ## GSEA

# %% [markdown]
# ### ranked vectors

# %%
groups <- unique(DEGs$group)

# %%
for (group in groups) {
    data <- DEGs[DEGs$group == group, ]
    data$fcsign <- sign(data$logFC)
    data$logP <- -log10(data$PValue) 
    data$metric <- data$logP/data$fcsign 
    data <- data[order(data$metric, decreasing= TRUE), ]
    gene_vector <- data$metric
    names(gene_vector) <- data$gene
    #assign(paste0('table_ranked_', group), data) #generated the table just to check that it was all ok but then the vector is enough
    assign(paste0('ranked_', group), gene_vector)
  }

# %% [markdown]
# ### gene lists

# %% [markdown]
# #### import Hallmarks

# %%
hallmark_genesets_path <- "~/giovanna/references/h.all.v2023.2.Hs.symbols.gmt"
reactome_genesets_path <- "~/giovanna/references/c2.cp.reactome.v2023.2.Hs.symbols.gmt"

# %%
hallmarks <- fgsea::gmtPathways(hallmark_genesets_path)

# %% [markdown]
# #### import genes from Hackert et al
# get to their "core inflammatory programme" as performed in the paper https://www.nature.com/articles/s41467-023-43573-9

# %%
hackert_allgenes <- read_delim('~/giovanna/references/CHIP27/Hackert2023_neut_inflamm.txt', delim="\t")

# %%
hackert_allgenes <- hackert_allgenes %>% arrange(fisher_adjusted)

# %% [markdown]
# the table was already ranked

# %%
hackert_top500<- hackert_allgenes %>% slice_head(n = 500)

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
hallmarks$HACKERT_NEUTS_INFLAMMATION_UP <- hackert_selected_up$symbol

# %% [markdown]
# #### import genes from Van Galen et al
# https://www.sciencedirect.com/science/article/pii/S0092867419300947?via%3Dihub
# <br> from supplementary tables

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
# https://www.science.org/doi/10.1126/sciimmunol.aay6017 <br>
# already human, filtered for abs(logFC)>1 and adjusted P value < 0.05 

# %%
Alsh_GMDCs <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/references/Alshetaiwi2020/20250204_aay6017_table_s3_GMDSCs_versus_neutrophils_filtered_human.csv')

# %%
Alsh_GMDCs %>% filter(avg_logFC>1, p_val_adj<0.05) %>% pull(external_gene_name) %>% unique() %>% length()

# %%
hallmarks$ALSH_GMDCs_UP <- Alsh_GMDCs %>% filter(avg_logFC>1, p_val_adj<0.05) %>% unique() %>% pull(external_gene_name)

# %%
hallmarks$ALSH_GMDCs_DOWN <- Alsh_GMDCs %>% filter(avg_logFC<(-1), p_val_adj<0.05) %>% unique() %>% pull(external_gene_name)

# %% [markdown]
# #### import maturation genes Montaldo
# https://www.nature.com/articles/s41590-022-01311-1

# %%
Montado <- read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/ref/Montado2022_suppltab29.csv')

# %%
hallmarks$MONTADO_MATURE_BM <- Montado %>% filter(Neutrophil.maturation.stage=='Mature BM neutrophils')  %>% pull(Genes)

# %% [markdown]
# #### import maturation genes Kirchberger-Shoeb
# https://www.nature.com/articles/s41467-024-45802-1

# %%
kirk<-read.csv('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/ref/20250518_Kirchberger2024_pansignature_human_GM.csv')

# %%
hallmarks$KIRK_M3<- kirk %>% filter(gene_module == 'M3') %>% pull(genes)

# %% [markdown]
# ## visualise

# %%
gsea_hall<-list()

for (group in groups) {
      ranked_list <- get(paste0('ranked_', group))
      cluster <- str_split(group, "_(?=[^_]+$)", n = 2)[[1]][1]
      treatment <- str_split(group, "_(?=[^_]+$)", n = 2)[[1]][2]
      
      fgseaRes.hallmarks <- (fgsea(pathways=hallmarks, stats=ranked_list, minSize=15, maxSize=400, nperm=10000) %>% as.data.frame)[1:5]
      fgseaRes.hallmarks.to.plot <- fgseaRes.hallmarks[fgseaRes.hallmarks$padj<0.05,]
      fgseaRes.hallmarks.to.plot$pathway <- gsub("HALLMARK_","", fgseaRes.hallmarks.to.plot$pathway)
      fgseaRes.hallmarks.to.plot$cluster <- cluster
      fgseaRes.hallmarks.to.plot$treatment <- treatment
    
      assign(paste0("gsea_hall_", group), fgseaRes.hallmarks.to.plot)
      gsea_hall[[group]] <- fgseaRes.hallmarks.to.plot
      }
    

gsea_hallmarks <- bind_rows(gsea_hall, .id="group")


# %% [markdown]
# ### bubble plot

# %%
gsea_hallmarks <- gsea_hallmarks %>%
                  mutate(treatment= as_factor(treatment)) %>%
                  mutate(cluster=factor(cluster, levels= c('2_5','1_7','0','6','3_4'))) %>%
                  arrange(cluster, treatment) %>%
                  mutate(group = factor(group, levels = unique(group)))

# %%
fig(20,15)
p0<- ggplot(gsea_hallmarks, aes(cluster, pathway, size=-log10(padj), col=NES)) + #col=logFC_capped
    geom_point() + 
    scale_color_gradient2(low='deepskyblue4', mid = 'gray90', high='firebrick3', name='NES') +
    xlab('') + ylab('') +
    geom_vline(xintercept = seq(1, 4, 1) + 0.5, alpha = 0.2)+
    #geom_vline(xintercept = seq(1, 4, 1) + 0.5, color = "grey20", alpha = 1, linewidth=0.5)+
    ggtitle(paste0('DNMT3A MUT vs WTx3, MUT:IFN, IFN') ) + 
    facet_wrap(~treatment, nrow=1)+
    theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=0.5),
         axis.text=element_text(size=15, color='black'),
         text=element_text(size=15),
            panel.border = element_rect(color='black', linewidth=1, fill = NA),
            panel.background = element_blank(),
            axis.line = element_blank(),
          legend.position = 'bottom',
          strip.background = element_blank(),
          strip.text = element_text(size = 15, face = "bold", color = "black"),
          legend.title = element_text(vjust=1),
        plot.margin=margin(1, 3, 1, 1, "lines")
)
p0

# %%
sel_path<- c('E2F_TARGETS','MYC_TARGETS_V1','MTORC1_SIGNALING', 'UNFOLDED_PROTEIN_RESPONSE', 'G2M_CHECKPOINT',
                       'APOPTOSIS', 'P53_PATHWAY', 
                       'TNFA_SIGNALING_VIA_NFKB', 'INFLAMMATORY_RESPONSE', 'INTERFERON_ALPHA_RESPONSE', 
                       'INTERFERON_GAMMA_RESPONSE', 'IL6_JAK_STAT3_SIGNALING', 
                       'HACKERT_NEUTS_INFLAMMATION_UP', 'HYPOXIA', 
                       'OXIDATIVE_PHOSPHORYLATION', 'CHOLESTEROL_HOMEOSTASIS', 'FATTY_ACID_METABOLISM') %>% rev()

# %%
fig(10,8)
p2<- ggplot(gsea_hallmarks %>% filter(treatment %in% c('untreated', 'treated')) %>%
            filter(pathway %in% sel_path) %>% 
            mutate(pathway = factor(pathway, levels = sel_path)),
            aes(cluster, pathway, size=-log10(padj), col=NES)) +
    geom_point() + 
    scale_color_gradient2(low='deepskyblue4', mid = 'gray90', high='firebrick3', name='NES') +
    xlab('') + ylab('') +
    geom_vline(xintercept = seq(1, 4, 1) + 0.5, alpha = 0.2)+
    #geom_vline(xintercept = seq(1, 4, 1) + 0.5, color = "grey20", alpha = 1, linewidth=0.5)+
    ggtitle(paste0('DNMT3A MUT vs WT') ) + 
    facet_wrap(~treatment, nrow=1)+
    theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=0.5),
         axis.text=element_text(size=15, color='black'),
         text=element_text(size=15),
            panel.border = element_rect(color='black', linewidth=1, fill = NA),
            panel.background = element_blank(),
            axis.line = element_blank(),
          legend.position = 'bottom',
          strip.background = element_blank(),
          strip.text = element_text(size = 15, face = "bold", color = "black"),
          legend.title = element_text(vjust=1),
        plot.margin=margin(1, 3, 1, 1, "lines")
)
p2

# %%
#interferon treatment transcription check

sel_path_infl<- c('E2F_TARGETS','MYC_TARGETS_V1','G2M_CHECKPOINT',
                       'TNFA_SIGNALING_VIA_NFKB', 'INFLAMMATORY_RESPONSE', 'INTERFERON_ALPHA_RESPONSE', 
                       'INTERFERON_GAMMA_RESPONSE', 'IL6_JAK_STAT3_SIGNALING', 'HACKERT_NEUTS_INFLAMMATION_UP','HYPOXIA') %>% rev()
9.5*length(sel_path_infl)/length(sel_path) #for fig size

fig(8.5,6)
p3<- ggplot(gsea_hallmarks %>% filter(treatment %in% c('allIFNy')) %>%
            filter(pathway %in% sel_path_infl) %>% 
            mutate(pathway = factor(pathway, levels = sel_path)),
            aes(cluster, pathway, size=-log10(padj), col=NES)) + #col=logFC_capped
    geom_point() + 
    scale_color_gradient2(low='deepskyblue4', mid = 'gray90', high='firebrick3', name='NES') +
    xlab('') + ylab('') +
    geom_vline(xintercept = seq(1, 4, 1) + 0.5, alpha = 0.2)+
    #geom_vline(xintercept = seq(1, 4, 1) + 0.5, color = "grey20", alpha = 1, linewidth=0.5)+
    ggtitle(paste0('DNMT3A MUT vs WT') ) + 
    facet_wrap(~treatment, nrow=1)+
    theme(axis.text.x=element_text(angle=-45, hjust=0, vjust=0.5),
         axis.text=element_text(size=15, color='black'),
         text=element_text(size=15),
            panel.border = element_rect(color='black', linewidth=1, fill = NA),
            panel.background = element_blank(),
            axis.line = element_blank(),
          legend.position = 'bottom',
          strip.background = element_blank(),
          strip.text = element_text(size = 15, face = "bold", color = "black"),
          legend.title = element_text(vjust=1),
        plot.margin=margin(1, 8, 1, 1, "lines")
)
p3

# %% [markdown]
# ### GSEA plot on neut cluster

# %%
fgsea_try <- fgsea(pathways=hallmarks, stats=ranked_6_treated, minSize=15, maxSize=400, nperm=10000) #I am only interested in cl6

# %%
plot_pathway<-fgsea_try$pathway[51]


pd <- plotEnrichmentData(
    pathway = hallmarks[[plot_pathway]],
    stats = ranked_6_treated
)

pval <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(padj) %>% round(6)

NES <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(NES) %>% round(3)

p1<- with(pd,
     ggplot(data=curve) +
         geom_line(aes(x=rank, y=ES), color="chartreuse3", size=2) +
         geom_ribbon(data=stats,
                     mapping=aes(x=rank, ymin=0,
                                 ymax=stat/maxAbsStat*(spreadES/4)),
                     fill="grey") +
         geom_segment(data=ticks,
                      mapping=aes(x=rank, y=-spreadES/16,
                                  xend=rank, yend=spreadES/16),
                      size=0.4) +
         geom_hline(yintercept=posES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=negES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=0, colour="black") +
         annotate("text", x=6000, y=0.55, label=paste0('padj=',pval), size=8, color="gray26") +
         annotate("text", x=6000, y=0.45, label=paste0('NES=',NES), size=8, color="gray26") + 
         theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
         ) +
         theme_minimal(base_size = 24)+
         labs(x="rank", y="enrichment score", title=plot_pathway))
p1

# %%
plot_pathway<-fgsea_try$pathway[52]


pd <- plotEnrichmentData(
    pathway = hallmarks[[plot_pathway]],
    stats = ranked_6_treated
)

pval <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(padj) %>% round(6)

NES <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(NES) %>% round(3)

p2<- with(pd,
     ggplot(data=curve) +
         geom_line(aes(x=rank, y=ES), color="chartreuse3", size=2) + 
         geom_ribbon(data=stats,
                     mapping=aes(x=rank, ymin=0,
                                 ymax=stat/maxAbsStat*(spreadES/4)),
                     fill="grey") +
         geom_segment(data=ticks,
                      mapping=aes(x=rank, y=-spreadES/16,
                                  xend=rank, yend=spreadES/16),
                      size=0.4) +
         geom_hline(yintercept=posES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=negES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=0, colour="black") +
         annotate("text", x=3000, y=-0.6, label=paste0('padj=',pval), size=8, color="gray26") + 
         annotate("text", x=3000, y=-0.75, label=paste0('NES=',NES), size=8, color="gray26") +
         theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
         ) +
         theme_minimal(base_size = 24)+
         labs(x="rank", y="enrichment score", title=plot_pathway))+
         coord_cartesian(ylim = c(NA, 0.0))
p2

# %%
fig(5,5)


plot_pathway<-fgsea_try$pathway[50]


pd <- plotEnrichmentData(
    pathway = hallmarks[[plot_pathway]],
    stats = ranked_6_treated
)

pval <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(padj) %>% round(6)

NES <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(NES) %>% round(3)

p3<- with(pd,
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
         geom_hline(yintercept=posES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=negES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=0, colour="black") +
         annotate("text", x=6000, y=0.55, label=paste0('padj=',pval), size=8, color="gray26") + 
         annotate("text", x=6000, y=0.46, label=paste0('NES=',NES), size=8, color="gray26") + 
         theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
         ) +
         theme_minimal(base_size = 24)+
         labs(x="rank", y="enrichment score", title=plot_pathway))
p3

# %%
plot_pathway<-fgsea_try$pathway[54]


pd <- plotEnrichmentData(
    pathway = hallmarks[[plot_pathway]],
    stats = ranked_6_treated
)

pval <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(padj) %>% round(6)

NES <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(NES) %>% round(3)

p4<- with(pd,
     ggplot(data=curve) +
         geom_line(aes(x=rank, y=ES), color="chartreuse3", size=2) + 
         geom_ribbon(data=stats,
                     mapping=aes(x=rank, ymin=0,
                                 ymax=stat/maxAbsStat*(spreadES/4)),
                     fill="grey") +
         geom_segment(data=ticks,
                      mapping=aes(x=rank, y=-spreadES/16,
                                  xend=rank, yend=spreadES/16),
                      size=0.4) +
         geom_hline(yintercept=posES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=negES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=0, colour="black") +
         annotate("text", x=6880, y=0.38, label=paste0('padj=',pval), size=8, color="gray26") +
         annotate("text", x=6880, y=0.25, label=paste0('NES=',NES), size=8, color="gray26") +
         theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
         ) +
         theme_minimal(base_size = 24)+
         labs(x="rank", y="enrichment score", title=plot_pathway))
p4

# %%
plot_pathway<-fgsea_try$pathway[55]


pd <- plotEnrichmentData(
    pathway = hallmarks[[plot_pathway]],
    stats = ranked_6_treated
)

pval <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(padj) %>% round(6)

NES <- fgsea_try %>% 
  filter(pathway == plot_pathway) %>% 
  pull(NES) %>% round(3)

p5<- with(pd,
     ggplot(data=curve) +
         geom_line(aes(x=rank, y=ES), color="chartreuse3", size=2) + 
         geom_ribbon(data=stats,
                     mapping=aes(x=rank, ymin=0,
                                 ymax=stat/maxAbsStat*(spreadES/4)),
                     fill="grey") +
         geom_segment(data=ticks,
                      mapping=aes(x=rank, y=-spreadES/16,
                                  xend=rank, yend=spreadES/16),
                      size=0.4) +
         geom_hline(yintercept=posES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=negES, colour="red", linetype="dashed", size=1.2) + 
         geom_hline(yintercept=0, colour="black") +
         annotate("text", x=6880, y=0.48, label=paste0('padj=',pval), size=8, color="gray26") + 
         annotate("text", x=6880, y=0.38, label=paste0('NES=',NES), size=8, color="gray26") + 
         theme(
             panel.background = element_blank(),
             panel.grid.major=element_line(color="grey92")
         ) +
         theme_minimal(base_size = 24)+
         labs(x="rank", y="enrichment score", title=plot_pathway))
p5
