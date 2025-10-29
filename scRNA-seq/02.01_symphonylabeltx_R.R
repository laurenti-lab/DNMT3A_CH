# %% [markdown]
# # B_02.04_CHIP27 Symphony label transfer
# Author:GM <br>
# From: 
# > https://htmlpreview.github.io/?https://github.com/andygxzeng/BoneMarrowMap_Extras/blob/main/BoneMarrowMap_Tutorial.nb.html

# %%
suppressPackageStartupMessages({library(Seurat)
                                library(tidyverse)
                                library(symphony)
                                library(patchwork)
                                library(ggpubr)
                                library(BoneMarrowMap)})

# %%
cat("R version:", R.version$version.string, "\n")

installed.packages() %>%
  as_tibble() %>%
  select(Package, Version) %>%
  filter(Package %in% c('Seurat', 'tidyverse', 'symphony', 'patchwork','ggpubr', 'BoneMarrowMap')) %>%
  print(row.names = FALSE)

# %%
basedir <- paste(laurenti, user, project, sep='/')

outdir<-paste(basedir, 'output', 'tables', '02.04', sep='/')

# %% [markdown]
# ### load reference

# %%
ref_path="/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/references/Zeng_BoneMarrowMap/"
ref <- readRDS(paste0(ref_path,"BoneMarrow_RefMap_SymphonyRef.rds"))
ref$save_uwot_path <- paste0(ref_path,"BoneMarrow_RefMap_uwot_model.uwot")

# %% [markdown]
# #### visualise reference

# %%
ReferenceSeuratObj <- create_ReferenceObject(ref)

# %%
options(repr.plot.width = 18, repr.plot.height = 8)
DimPlot(ReferenceSeuratObj, reduction = 'umap', group.by = 'CellType_Annotation_formatted', raster=FALSE, label=TRUE, label.size = 4)

# %% [markdown]
# ### import query data

# %%
ann <- anndata::read_h5ad('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/output/objects/20241008_CHIP27_countsforazimuth.h5ad')

# %%
query <- CreateSeuratObject(counts = t( as.matrix(ann$X) ),      
                           meta.data = ann$obs)

# %%
Sys.time()

batchvar <- 'library'

# Map query dataset using Symphony (Kang et al 2021)
query <- map_Query(
    exp_query = query@assays$RNA@counts, 
    metadata_query = query@meta.data,
    ref_obj = ref,
    vars = batchvar
)

Sys.time()


# %%
# Run QC based on mapping error score, flag cells with mapping error >= 3 MADs above median
query <- query %>% calculate_MappingError(., reference = ref, MAD.threshold = 3)
QC_plots <- plot_MappingErrorQC(query)
patchwork::wrap_plots(QC_plots, ncol = 4, widths = c(0.8, 0.3, 0.8, 0.3))

# %%
Sys.time()

# Predict Hematopoietic Cell Types by KNN classification
query <- predict_CellTypes(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_CellType', # celltype assignments before filtering on mapping QC I don't have this but it's ok
  final_label = 'predicted_CellType'  # celltype assignments with map QC failing cells assigned as NA
) 

DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap', group.by = c('predicted_CellType'), raster=FALSE, label=TRUE, label.size = 4)

Sys.time()

# %%
options(repr.plot.width = 14, repr.plot.height = 10)
DimPlot(subset(query, mapping_error_QC == 'Pass'), reduction = 'umap', group.by = c('predicted_CellType'), raster=FALSE, label=TRUE, label.size = 4)

# %%
Sys.time()

# Predict Pseudotime values by KNN
query <- predict_Pseudotime(
  query_obj = query, 
  ref_obj = ref, 
  initial_label = 'initial_Pseudotime', 
  final_label = 'predicted_Pseudotime'
)

Sys.time()

# %%
# Visualize Hematopoietic Pseudotime in query data
options(repr.plot.width = 10, repr.plot.height = 10)

FeaturePlot(subset(query, mapping_error_QC == 'Pass'), features = c('predicted_Pseudotime'))

# %%
# Save CellType Annotations and Projected UMAP coordinates
save_ProjectionResults(
  query_obj = query, 
  celltype_label = 'predicted_CellType', 
  celltype_KNNprob_label = 'predicted_CellType_prob', 
  pseudotime_label = 'predicted_Pseudotime', 
  file_name = paste0(outdir, '/', prefix, '_CHIP27_vs_andyBM.csv'))

# %%
# Set batch/condition to be visualized individually
batch_key <- 'predicted_CellType'

# returns a list of plots for each donor from a pre-specified batch variable
projection_plots <- plot_Projection_byDonor(
  query_obj = query, 
  batch_key = batch_key, 
  ref_obj = ref, 
  save_folder = '../output/projections/'
)

# show plots together with patchwork
patchwork::wrap_plots(projection_plots)

# %%
options(repr.plot.width = 18, repr.plot.height = 18)

patchwork::wrap_plots(projection_plots)

# %%
query_composition <- get_Composition(
  query_obj = query, 
  donor_key = 'library', 
  celltype_label = 'predicted_CellType', 
  mapQC_col = 'mapping_error_QC', 
  knn_prob_cutoff = NULL, 
  return_type = 'long')


# %%
write.table(file = paste0(outdir, '/', prefix, 'CHIP27_vs_andyBM_ref_annots_size.tsv'), sep='\t', query_composition)
