# %% [markdown]
# # 04.01_CellChat
# useful links: https://rdrr.io/github/sqjin/CellChat/f/tutorial/CellChat-vignette.Rmd <br>
# Author:GM <br>

# %%
suppressMessages({
    library(CellChat)
    library(patchwork)
})

options(stringsAsFactors = FALSE)

# %%
cat("R version:", R.version$version.string, "\n")

installed.packages() %>%
  as_tibble() %>%
  select(Package, Version) %>%
  filter(Package %in% c('anndata', 'CellChat')) %>%
  print(row.names = FALSE)

# %%
basedir <- paste(laurenti, user, project, sep='/')
basename <- '_CHIP27_'
print(basedir)

# %%
datadir<-paste(basedir, 'output', sep='/')
figdir<- paste(datadir, 'figures', '10.01/', sep='/')
refdir<- paste(basedir, 'ref', sep='/')

# %% [markdown]
# ### fxs

# %%
options(repr.matrix.max.cols=Inf)

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %% [markdown]
# ### import and prep data

# %%
ann.path <-'/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP27/output/objects/'

# %%
ann <- anndata::read_h5ad(paste0(ann.path, '20250814_CHIP27_after_lognorm_BBKNN_Kwok_Zeng_cellcycle_clustering1115_1120_pseudotime.h5ad'))

# %% [markdown]
# ##### create grouping info

# %%
ann$obs<- ann$obs %>%
  mutate(cluster_1120_geno = paste0(cluster_1120, "_", DNMT3A_genotype))

# %% [markdown]
# ##### isolate only clusters of interest
# clusters with enough cells

# %%
ann_filt <- ann[ann$obs$cluster_1120 %in% c('2_5','1_7','0','3_4','6'), ]

# %% [markdown]
# ##### segregate untreated and treated datasets

# %%
ann_CTR <- ann_filt[ann_filt$obs$treatment == "CTR", ]

# %%
ann_IFN <- ann_filt[ann_filt$obs$treatment == "IFNy", ]

# %% [markdown]
# ## process untreated data

# %% [markdown]
# ### create cellchat object

# %%
# access raw counts
counts <- t( as.matrix(ann_CTR$layers$get("counts")) )

# normalize the count data as they do in the vignette (probably lognorm would work but let's stick to it)
library.size <- Matrix::colSums(counts)
data.input <- as(log1p(Matrix::t(Matrix::t(counts)/library.size) * 10000), "dgCMatrix")

# %%
cellchat_CTR <- createCellChat(object = data.input, 
                           meta = ann_CTR$obs, 
                           group.by = "cluster_1120_geno")

# %%
saveRDS(cellchat_CTR, file = paste0("../output/objects/", prefix, '_cellchat_CTR_object.rds'))

# %% [markdown]
# ### process - Secreted Signalling

# %%
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# set the used database in the object
cellchat_CTR@DB <- CellChatDB.use

# %%
# subset the expression data of signaling genes for saving computation cost
cellchat_CTR <- subsetData(cellchat_CTR) # This step is necessary even if using the whole database
future::plan("multisession", workers = 12) # do parallel

# %%
ptm = Sys.time()
cellchat_CTR <- identifyOverExpressedGenes(cellchat_CTR) #, do.fast = FALSE)
cellchat_CTR <- identifyOverExpressedInteractions(cellchat_CTR)

# %%
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# %%
options(future.globals.maxSize = 1300 * 1024^2)

# %%
ptm = Sys.time()

cellchat_CTR <- computeCommunProb(cellchat_CTR, type = "triMean", population.size = TRUE)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# %%
df.net <- subsetCommunication(cellchat_CTR)

# %%
write.table(file= paste0('../output/tables/', prefix, 'cellchat_CTR_inferred_cell_cell_comms.txt'), df.net)

# %%
cellchat_CTR <- computeCommunProbPathway(cellchat_CTR)

# %%
ptm = Sys.time()

cellchat_CTR <- aggregateNet(cellchat_CTR)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# %%
saveRDS(cellchat_CTR, file = paste0("../output/objects/", prefix, '_cellchat_CTR_object_processed.rds'))

# %% [markdown]
# ### plot

# %%
groupSize <- as.numeric(table(cellchat_CTR@idents))

# %%
levels(cellchat_CTR@idents)

# %%
pathways.show<-cellchat_CTR@netP$pathways

# %%
i <- 1
for (pathway in pathways.show) {
  fname <- paste0(figdir, prefix, "_CTR_p", i, ".pdf")
  grDevices::pdf(fname, width = 8, height = 8)
  netVisual_chord_cell(
    cellchat_CTR,
    signaling = pathway,
    group = group.cellType,
    title.name = paste0('CTR_',pathway))
  rp <- recordPlot()
  grDevices::dev.off()
  assign(paste0("_CTR_p", i), rp, envir = .GlobalEnv)
  i <- i + 1
}


# %% [markdown]
# ## process IFN data

# %% [markdown]
# ### create cellchat object

# %%
# access raw counts
counts <- t( as.matrix(ann_IFN$layers$get("counts")) )

# normalize the count data as they do in the vignette (probably lognorm would work but let's stick to it)
library.size <- Matrix::colSums(counts)
data.input <- as(log1p(Matrix::t(Matrix::t(counts)/library.size) * 10000), "dgCMatrix")

# %%
cellchat_IFN <- createCellChat(object = data.input, 
                           meta = ann_IFN$obs, 
                           group.by = "cluster_1120_geno")

# %%
saveRDS(cellchat_IFN, file = paste0("../output/objects/", prefix, '_cellchat_IFN_object.rds'))

# %% [markdown]
# ### process - Secreted signaling

# %%
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling

# set the used database in the object
cellchat_IFN@DB <- CellChatDB.use

# %%
# subset the expression data of signaling genes for saving computation cost
cellchat_IFN <- subsetData(cellchat_IFN) # This step is necessary even if using the whole database
future::plan("multisession", workers = 12) # do parallel

# %%
ptm = Sys.time()
cellchat_IFN <- identifyOverExpressedGenes(cellchat_IFN) #, do.fast = FALSE)
cellchat_IFN <- identifyOverExpressedInteractions(cellchat_IFN)

# %%
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# %%
options(future.globals.maxSize = 700 * 1024^2)

# %%
ptm = Sys.time()

cellchat_IFN <- computeCommunProb(cellchat_IFN, type = "triMean", population.size = TRUE)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# %%
df.net <- subsetCommunication(cellchat_IFN)

# %%
write.table(file= paste0('../output/tables/', prefix, 'cellchat_IFN_inferred_cell_cell_comms.txt'), df.net)

# %%
cellchat_IFN <- computeCommunProbPathway(cellchat_IFN)

# %%
ptm = Sys.time()

cellchat_IFN <- aggregateNet(cellchat_IFN)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# %%
saveRDS(cellchat_IFN, file = paste0("../output/objects/", prefix, '_cellchat_IFN_object_processed.rds'))

# %%
#explore
unique( df.net[df.net$source == '6_MUT',]$pathway_name )

# %% [markdown]
# #### plotting

# %%
groupSize <- as.numeric(table(cellchat_IFN@idents))


# %%
levels(cellchat_IFN@idents)

# %%
pathways.show<-cellchat_IFN@netP$pathways

# %%
i <- 1
for (pathway in pathways.show) {
  fname <- paste0(figdir, prefix, "_IFN_p", i, ".pdf")
  grDevices::pdf(fname, width = 8, height = 8)
  netVisual_chord_cell(
    cellchat_IFN,
    signaling = pathway,
    group = group.cellType,
    title.name = paste0('IFN_',pathway))
  rp <- recordPlot()
  grDevices::dev.off()
  assign(paste0("_IFN_p", i), rp, envir = .GlobalEnv)
  i <- i + 1
}

