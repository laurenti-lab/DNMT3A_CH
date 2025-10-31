# %% [markdown]
# # 01. Extract DNMT3A-mut and DNMT3A-wt PDIDIs of the colonies
#

# %%
Sys.getenv('SINGULARITY_CONTAINER')
Sys.getenv('PYTHONUSERBASE')

# %%
tidyverse <- c('dplyr', 'readr', 'ggplot2', 'readr', 'tidyr', 'purrr', 'tibble', 'stringr', 'forcats')
lapply(tidyverse, library, character.only = TRUE)

library('ape')
library('phangorn')
library('ggtree')
library('treeio')

# %% [markdown]
# R version: R version 4.3.1 (2023-06-16) 
# A tibble: 5 × 2
#   Package  Version
#   <chr>    <chr>  
# 1 ape      5.7-1  
# 2 dplyr    1.1.4  
# 3 ggtree   3.10.0 
# 4 phangorn 2.11.1 
# 5 treeio   1.26.0 

# %%
system("echo $SLURM_TASKS_PER_NODE", intern = TRUE)

# %%
slurm_nprocs <- Sys.getenv("$SLURM_NPROCS")
slurm_nprocs

# %%
slurm_job_partition <- Sys.getenv("SLURM_JOB_PARTITION")
slurm_job_partition

# %% [markdown]
# ---

# %% [markdown]
# ## Set dir and load object

# %%
outputdir <- ('/.../CHIP28/output')

# %%
today <- format(Sys.Date(), "%Y%m%d")

# %%
load('/.../annotated_mut_set_LEUK4_1_01_standard_rho01')

# %%
tree_raw <- read.newick("/.../tree_LEUK4_1_01_standard_rho01.tree")

# %% [markdown]
# ## find DNMT3A colony IDs

# %% [markdown]
# ### find DNMT3A mutant node number

# %%
#create a smaller object d (data) to work on
d <- filtered_muts$COMB_mats.tree.build$mat

# %% [markdown]

# %%
DNMT3A_R882H <- d %>%
                filter(Chrom == 2,
                       Pos %in% 25457241:25457243)
head(DNMT3A_R882H)

# %%
DNMT3A_R882H$node

# %% [markdown]
# #### visualise tree

# %%
ggtree(tree_raw) +
   geom_highlight(data = DNMT3A_R882H,
                 aes(node = node),
                 fill = "lightgoldenrod1") +
   geom_tiplab(hjust = 0, size = 1.5)

# %% [markdown]
# #### get node IDs of the tips

# %% jupyter={"outputs_hidden": true}
DNMT3A_tips <- offspring(tree_raw, DNMT3A_R882H$node, tiponly=TRUE)

# %% [markdown]
# #### get label IDs of the tips

# %%
DNMT3A_tips_table <- tree_raw %>%
                     as_tibble %>%
                     filter(node %in% DNMT3A_tips) 

# %%
DNMT3A_tips_list <- DNMT3A_tips_table$label

# %% [markdown]
# ### tips counts

# %%
#total tips of the tree:
Ntip(tree_raw)

# %% [markdown]
# **DNMT3A mutant tips:**

# %%
length(DNMT3A_tips_list)

# %% [markdown]
# **non mutant tips:**

# %%
Ntip(tree_raw)-length(DNMT3A_tips_list)

# %% [markdown]
# ## find other myeloid drivers node IDs

# %% [markdown]
# ### TET2

# %%
TET2_all <- d %>%
            filter(Gene == 'TET2')

# %%
ggtree(tree_raw) +
   geom_highlight(data = TET2_all,
                 aes(node = node),
                 fill = "blue") +
   geom_tiplab(hjust = 0, size = 1.5)

# %% [markdown]
# Exclude the synonimous variants and the ones in the intronic regions:

# %%
TET2_M1333K <- d %>%
            filter(Gene == 'TET2') %>%
            filter(Protein =='p.M1333K')

# %%
ggtree(tree_raw) +
   geom_highlight(data = TET2_M1333K,
                 aes(node = node),
                 fill = "skyblue") +
   geom_tiplab(hjust = 0, size = 1.5)

# %%
ggtree(tree_raw) +
   geom_highlight(data = DNMT3A_R882H,
                  aes(node = node),
                  fill = "lightgoldenrod1") +
   geom_highlight(data = TET2_M1333K,
                 aes(node = node),
                 fill = "skyblue") +
   geom_tiplab(hjust = 0, size = 1.5)

# %%
TET2_tips <- offspring(tree_raw, TET2_M1333K$node, tiponly=TRUE)

# %% [markdown]
# ### JAK2

# %%
filtered_muts$COMB_mats.tree.build$mat %>%
            filter(Gene == 'JAK2')

# %%
JAK2_V617F <- d %>%
              filter(Gene == 'JAK2') %>%
              filter(Protein == 'p.V617F')

# %%
ggtree(tree_raw) +
   geom_highlight(data = JAK2_V617F,
                 aes(node = node),
                 fill = "olivedrab3") +
   geom_tiplab(hjust = 0, size = 1.5)

# %%
ggtree(tree_raw) +
   geom_highlight(data = DNMT3A_R882H,
                  aes(node = node),
                  fill = "lightgoldenrod1") +
   geom_highlight(data = TET2_M1333K,
                 aes(node = node),
                 fill = "skyblue") +
   geom_highlight(data = JAK2_V617F,
                 aes(node = node),
                 fill = "olivedrab3") +
   geom_tiplab(hjust = 0, size = 1.5)

# %%
JAK2_tips <- offspring(tree_raw, JAK2_V617F$node, tiponly=TRUE)

# %% [markdown]
# ### ASXL1
#

# %%
d %>%       filter(Gene == 'ASXL1')

# %%
ASXL1_fs <- d %>%
            filter(Gene == 'ASXL1') %>%
            filter(Protein =='p.S1079fs*8')

# %%
ggtree(tree_raw) +
   geom_highlight(data = ASXL1_fs,
                 aes(node = node),
                 fill = "orange") +
   geom_tiplab(hjust = 0, size = 1.5)

# %%
ASXL1_tips <- offspring(tree_raw, ASXL1_fs$node, tipsonly=TRUE)
length(ASXL1_tips)

# %%
ASXL1_label <- nodelab(tree_raw, ASXL1_fs$node)

# %% [markdown]
# ## visualise and export tree with final annotations:
#

# %%
fig1 <- ggtree(tree_raw) +
   geom_highlight(data = DNMT3A_R882H,
                  aes(node = node),
                  fill = "lightgoldenrod1") +
   geom_highlight(data = TET2_M1333K,
                 aes(node = node),
                 fill = "skyblue") +
   geom_highlight(data = JAK2_V617F,
                 aes(node = node),
                 fill = "olivedrab3") +
   geom_highlight(data = ASXL1_fs,
                 aes(node = node),
                 fill = "orange") +
   geom_tiplab(hjust = 0, size = 2)
fig1

# %%
DNMT3A_only_tips <- DNMT3A_tips[!(DNMT3A_tips %in% c(TET2_tips, JAK2_tips))]

# %%
DNMT3A_only_tips_table <- tree_raw %>%
                     as_tibble %>%
                     filter(node %in% DNMT3A_only_tips )  %>%
                     mutate(mutation = 'DNMT3A_only')
dim(DNMT3A_only_tips_table)
head(DNMT3A_only_tips_table)

# %% [markdown]
# #### DNMT3A TET2

# %%
intersect(DNMT3A_tips, TET2_tips) %>% length()
length(TET2_tips)

# %%
DNMT3A_TET2_tips_table <- tree_raw %>%
                     as_tibble %>%
                     filter(node %in% TET2_tips ) %>%
                     mutate(mutation = 'DNMT3A_TET2')


# %% [markdown]
# #### DNMT3A JAK2

# %%
intersect(DNMT3A_tips, JAK2_tips) %>% length()
length(JAK2_tips)

# %%
DNMT3A_JAK2_tips_table <- tree_raw %>%
                     as_tibble %>%
                     filter(node %in% JAK2_tips ) %>%
                     mutate(mutation = 'DNMT3A_JAK2')


# %% [markdown]
# #### WT AXL1

# %%
WT_ASXL1_tips_table <- tree_raw %>%
                     as_tibble %>%
                     filter(label == ASXL1_label ) %>%
                     mutate(mutation = 'WT_ASXL1')


# %% [markdown]
# #### WT
# colonies not offspring of the DNMT3A R882H node 133

# %%
WT_tips <-  tree_raw %>%
            as_tibble %>%
            filter(!(node %in% DNMT3A_tips)) %>%
            filter(str_detect(label, 'PD'))


# %%
WT_tips_table <-     tree_raw %>%
                     as_tibble %>%
                     filter(node %in% WT_tips$node) %>%
                     filter(label != ASXL1_label) %>%
                     mutate(mutation = 'WT')


# %% [markdown]
# ## merge and export in a csv:

# %%
ARCH1_colonyinfo <- bind_rows(DNMT3A_only_tips_table, DNMT3A_TET2_tips_table, DNMT3A_JAK2_tips_table, WT_tips_table, WT_ASXL1_tips_table)
head(ARCH1_colonyinfo)
tail(ARCH1_colonyinfo)
dim(ARCH1_colonyinfo)


# %%
write.csv(ARCH1_colonyinfo, paste0(outputdir,'/tables/', today, '_ARCH1_colonies_mut_info.csv'), row.names = FALSE)
