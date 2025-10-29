# %% [markdown]
# # 02.02_Azimuth label transfer
#
# 1. Azimuth PBMC reference map
# 2. dataset from Kwok 2023 paper Nat Immunology
#
# Author:GM <br>
# Date started:2024/10/04 <br>
# useful links: 
# > https://satijalab.github.io/azimuth/reference/RunAzimuth.html
# > https://azimuth.hubmapconsortium.org/
# > https://zenodo.org/records/4546839
# > https://github.com/satijalab/azimuth/wiki/Azimuth-Reference-Format
# > https://satijalab.org/seurat/articles/seurat5_integration
# > https://satijalab.org/seurat/articles/seurat5_essential_commands

# %%
suppressPackageStartupMessages({library(Seurat)
                                library(Azimuth)
                                library(SeuratData)
                                library(patchwork)
                                library(ggplot2)
                                library(dplyr)})

# %%
cat("R version:", R.version$version.string, "\n")

installed.packages() %>%
  as_tibble() %>%
  select(Package, Version) %>%
  filter(Package %in% c('dplyr','Seurat', 'Azimuth', 'SeuratData', 'patchwork','ggplot2')) %>%
  print(row.names = FALSE)

# %%
basedir <- paste(laurenti, user, project, sep='/')
basename <- '_CHIP27_'
print(basedir)

# %%
outdir<-paste(basedir, 'output', 'tables', '02.02_labeltx', sep='/')

# %% [markdown]
# ### import data
# importing anndata with raw counts

# %%
ann <- anndata::read_h5ad('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/output/objects/20241008_CHIP27_countsforazimuth.h5ad')

# %%
target <- CreateSeuratObject(counts = t( as.matrix(ann$X) ),      
                           meta.data = ann$obs)

# %%
rm(ann)

# %% [markdown]
# ### reference: Kwok 2023 Nat immunology
# https://www.nature.com/articles/s41590-023-01490-5

# %%
reference.path = paste0( basedir, '/references/kwokWB/' )

# %%
Sys.time()

annotated <- RunAzimuth(query=target,
                        reference=reference.path
                       )

Sys.time()

# %%
options(repr.plot.width = 16, repr.plot.height = 8)

Idents(annotated) <- 'predicted.fine_annot'

p1 <- DimPlot(annotated, label = TRUE) + NoLegend() + ggtitle('fine')

Idents(annotated) <- 'predicted.broad_annot'

p2 <- DimPlot(annotated, label = TRUE) + NoLegend() + ggtitle('coarse')

p1|p2

# %% [markdown]
# #### export annotations

# %%
write.table(file = paste0(outdir, '/', prefix, '_CHIP27_vs_kwokWB_azimuth_metadata.txt'), 
            annotated@meta.data, sep = '\t')

# %% [markdown]
# ### reference: Azimuth pbmc
# https://azimuth.hubmapconsortium.org/references/

# %%
reference.path = paste0( basedir, '/references/human_pbmc/' )

# %%
Sys.time()

annotated <- RunAzimuth(query=target,
                        reference=reference.path
                       )

Sys.time()

# %%
options(repr.plot.width = 18, repr.plot.height = 8)

Idents(annotated) <- 'predicted.celltype.l1'

p1 <- DimPlot(annotated, label = TRUE) + NoLegend() + ggtitle('l1')

Idents(annotated) <- 'predicted.celltype.l2'

p2 <- DimPlot(annotated, label = TRUE) + NoLegend() + ggtitle('l2')

Idents(annotated) <- 'predicted.celltype.l3'

p3 <- DimPlot(annotated, label = TRUE) + NoLegend() + ggtitle('l3')


p1|p2|p3

# %%
write.table(file = paste0(outdir, '/', prefix, '_CHIP27_vs_pbmc_azimuth_metadata.txt'), 
            annotated@meta.data, sep = '\t')

# %% [markdown]
# ### reference: Azimuth bmmc

# %%
reference.path = paste0( basedir, '/references/human_bmmc/' )

# %%
Sys.time()

annotated <- RunAzimuth(query=target,
                        reference=reference.path
                       )

Sys.time()

# %%
options(repr.plot.width = 16, repr.plot.height = 8)

Idents(annotated) <- 'predicted.celltype.l1'

p1 <- DimPlot(annotated, label = TRUE) + NoLegend() + ggtitle('l1')

Idents(annotated) <- 'predicted.celltype.l2'

p2 <- DimPlot(annotated, label = TRUE) + NoLegend() + ggtitle('l2')


p1|p2

# %%
write.table(file = paste0(outdir, '/', prefix, '_CHIP27_vs_bmmc_azimuth_metadata.txt'), 
            annotated@meta.data, sep = '\t')
