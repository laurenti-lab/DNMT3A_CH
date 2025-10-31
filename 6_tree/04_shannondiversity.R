# %% [markdown]
# # 05b. Shannon diversity index calculation
# GM <br>
# useful links: https://www.nature.com/articles/s41586-022-04786-y 

# %%
tidyverse <- c('dplyr', 'readr', 'ggplot2', 'readr', 'tidyr', 'purrr', 'tibble', 'stringr', 'forcats')

invisible(lapply(tidyverse, function(pkg) {suppressPackageStartupMessages(library(pkg, character.only = TRUE))}))

suppressPackageStartupMessages({library('ape')
                                library('phangorn')
                                library('ggtree')
                                library('treeio')
                                library('phytools')
                                library('vegan') })

# %%
R version: R version 4.3.2 (2023-10-31) 
# A tibble: 7 × 2
  Package  Version
  <chr>    <chr>  
1 ape      5.7-1  
2 dplyr    1.1.4  
3 ggtree   3.10.0 
4 phangorn 2.11.1 
5 phytools 2.3-0  
6 treeio   1.26.0 
7 vegan    2.6-6.1

# %% [markdown]
# ---

# %% [markdown]
# ## Define functions

# %%
getTips<-function(tree,node){
getDescendants(tree,node)[which(getDescendants(tree,node)<(length(tree$tip.label)+1))]
                   }

# %%
fig<- function(width,height){options(repr.plot.width=width, repr.plot.height=height)}


# %% [markdown]
# ---

# %% [markdown]
# ## Set dir and load objects
# CHIP1(tree_raw); AX001 and KX003 from Mitchell et al (link above)

# %%
outputdir <- ('/home/gm686/giovanna/CHIP28/output')

# %%
figdir <- paste0(outputdir, '/figures/05/')

# %%
load('/home/gm686/giovanna/CHIP28/data/tree_data/annotated_mut_set_LEUK4_1_01_standard_rho01')

# %%
tree_raw <- read.newick("/home/gm686/giovanna/CHIP28/data/tree_data/tree_LEUK4_1_01_standard_rho01.tree")

# %%
AX001<- read.tree("/home/gm686/giovanna/CHIP28/data/tree_data/tree_AX001_4_01_standard_rho01.tree")

# %%
KX003<- read.tree("/home/gm686/giovanna/CHIP28/data/tree_data/tree_KX003_5_01_standard_rho01_final.tree")

# %% [markdown]
# ---

# %% [markdown]
# ### env info

# %%
today <- format(Sys.Date(), "%Y%m%d")

# %%
Sys.getenv('SINGULARITY_CONTAINER')
Sys.getenv('PYTHONUSERBASE')
Sys.getenv('SLURM_JOB_PARTITION')
Sys.getenv('SLURM_NPROCS')

# %% [markdown]
# ### find DNMT3A and WT mutant node numbers
#

# %%
head(filtered_muts$COMB_mats.tree.build$mat)

# %%
#create a smaller object d (data) to work on (not strictly necessary)
d <- filtered_muts$COMB_mats.tree.build$mat


# %%
DNMT3A_R882H <- d %>%
                filter(Chrom == 2,
                       Pos %in% 25457241:25457243)
head(DNMT3A_R882H)

# %%
DNMT3A_R882H$node

# %% [markdown]
# this matches the R markdown from EM

# %% [markdown]
# #### viz

# %%
options(repr.plot.width=7, repr.plot.height=9)
p0 + geom_cladelab(node=133, label="DNMT3A R882", angle=270, #textcolor = 'blue',
                   hjust='center', offset.text= 15 , barsize= 0.5, fontsize= 5)

# %% [markdown]
# #### tip numbers (MUT, WT, total)

# %%
DNMT3A_tips <- getTips(tree_raw, 133)

# %% [markdown]
# **total number of tips**

# %%
Ntip(tree_raw)

# %% [markdown]
# **number of mutant tips**

# %%
length(DNMT3A_tips)

# %% [markdown]
# **number of WT tips**

# %%
Ntip(tree_raw)-length(DNMT3A_tips)

# %%
# the root node is usually numbered as the num_tips + 1 (the indexing first goes through all the tips and then back to the node)
num_tips <- length(tree_raw$tip.label)


root_node <- num_tips + 1

print(root_node)


# %%
#I confirm that 128 is indeed the root node
options(repr.plot.width=7, repr.plot.height=9)
ggtree(tree_raw) + geom_cladelab(node=128, label="root", angle=270, #textcolor = 'blue',
                   hjust='center', offset.text= 15 , barsize= 0.5, fontsize= 5)

# %%
all_tips <- offspring(tree_raw, 128 , type='tips')
all_tips %>% length()
Ntip(tree_raw)

# %%
setdiff(all_tips, DNMT3A_tips) %>% length

# %%
WT_tips <- setdiff(all_tips, DNMT3A_tips)

# %% [markdown]
# #### tip IDs

# %% [markdown]
# ##### WT

# %%
tree_raw_WT_tips_table <- tree_raw %>%
                          as_tibble %>%
                          filter(node %in% WT_tips) 

# %% [markdown]
# ##### MUT

# %%
tree_raw_DNMT3A_tips_table <- tree_raw %>%
                          as_tibble %>%
                          filter(node %in% DNMT3A_tips) 

# %%
tree_raw_DNMT3A_tips_table %>% head

# %% [markdown]
# #### final visualisation

# %%
fig(10,13)
ggtree(tree_raw) + 
    geom_tiplab(color='firebrick', size=3)+
    geom_highlight(node = DNMT3A_tips,fill = "lightgoldenrod1") +
    geom_highlight(node = WT_tips, fill= "lightblue")+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %% [markdown]
# ### Shannon DI on the WT tree

# %% [markdown]
# #### isolate WT tree

# %%
WT_tree <- drop.tip(tree_raw, DNMT3A_tips)

# %%
WT_tree

# %%
fig(10,5)
ggtree(WT_tree) + 
    geom_tiplab(color='firebrick', size=3)+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %%
setequal(tree_raw_WT_tips_table$label, WT_tree$tip.label) #sanitycheck

# %% [markdown]
# #### set clade threshold on WT tree + calculate clade size

# %%
time_point = 60
nodeheights <- nodeHeights(WT_tree)
clades_post_cutoff=WT_tree$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]
clade_sizes=sapply(clades_post_cutoff,function(node) {length(getTips(WT_tree, node))}) 

# %%
clade_sizes

# %%
fig(12,7)
p1<- ggtree(WT_tree) + 
  geom_text(aes(label = ifelse(isTip, label, '')), hjust = -0.3, color = 'firebrick') +
  geom_highlight(node = clades_post_cutoff, fill = "lightgreen") +
  geom_text2(aes(subset = node %in% clades_post_cutoff, label = node), 
             vjust = -1, color = 'blue')+
    theme(plot.margin=margin(1, 5, 1, 1, "lines"))
p1

# %%
ggsave(paste0(figdir, today, "_p1_WT.pdf"), p1, dpi=300, width= 10, height=20)

# %% [markdown]
# #### calculate Shannon DI on WT tree

# %%
clades_post_cutoff %>% length()

# %%
table(clade_sizes)

# %%
Shannon_index_WT <- diversity(clade_sizes, index = "shannon", MARGIN = 1, base = exp(1))
Shannon_index_WT

# %%
num_cols=length(clade_sizes)
try_matrix <- matrix(NA, ncol = num_cols , nrow = 1)

try_matrix[, 1:num_cols] <- clade_sizes

try_matrix


# %%
Shannon_index <- diversity(try_matrix, index = "shannon", MARGIN = 1, base = exp(1))
Shannon_index

# %% [markdown]
# ---

# %% [markdown]
# ### Shannon on MUT tree

# %% [markdown]
# #### isolate MUT tree

# %%
DNMT3A_tree <- drop.tip(tree_raw, WT_tips)

# %%
DNMT3A_tree

# %%
fig(10,10)
ggtree(DNMT3A_tree) + 
    geom_tiplab(color='firebrick', size=3)+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %%
setequal(tree_raw_DNMT3A_tips_table$label, DNMT3A_tree$tip.label)#sanitycheck

# %% [markdown]
# #### generate random trees subsets of the DNMT3A tree that have 33 tips (equal number to WT tree)

# %%
DNMT3A_tips %>% length() - 33

# %%
set.seed(123)
tip_combinations <- replicate(100, sample(DNMT3A_tips, 61), simplify = FALSE)
tip_combinations[1:5]

# %%
MUT_tree_subsets <- list()
for(i in seq_along(tip_combinations)) {
      tip_set <- tip_combinations[[i]]
      MUT_tree_subset <- drop.tip(DNMT3A_tree, tip_set)
      tree_name <- paste0("MUT_tree_subset_", i)
      assign(tree_name, MUT_tree_subset)
      MUT_tree_subsets[[tree_name]] <- MUT_tree_subset
}

# %% [markdown]
# ##### few checks

# %%
MUT_tree_subsets %>% head(n=3)

# %%
fig(10,6)
ggtree(MUT_tree_subset_1) + 
    geom_tiplab(color='firebrick', size=3)+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %%
fig(10,6)
ggtree(MUT_tree_subset_100) + 
    geom_tiplab(color='firebrick', size=3)+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %% [markdown]
# #### calculate Shannon DI for these randomly generated trees

# %%
time_point=60

# %%
MUT_clades_post_cutoff <- list()
MUT_clade_sizes <- list()
MUT_Shannon_index <- list()

for(i in seq_along(MUT_tree_subsets)) {
      tree <- MUT_tree_subsets[[i]]
      nodeheights <- nodeHeights(tree)
    
      clades_post_cutoff=tree$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]
      clades_post_cutoff_name=paste0("clades_post_cutoff_tree_subset_", i)
      MUT_clades_post_cutoff[[clades_post_cutoff_name]] <- clades_post_cutoff

      clade_sizes=sapply(clades_post_cutoff,function(node) {length(getTips(tree, node))}) 
      clade_sizes_name=paste0("clade_sizes_tree_subset_", i)
      MUT_clade_sizes[[clade_sizes_name]] <- clade_sizes

      Shannon_index <- diversity(clade_sizes, index = "shannon", MARGIN = 1, base = exp(1))
      Shannon_index_name=paste0("Shannon_index_tree_subset_", i)
      MUT_Shannon_index[[Shannon_index_name]] <- Shannon_index
}

# %% jupyter={"outputs_hidden": true}
MUT_clades_post_cutoff %>% head (n=3)
MUT_clade_sizes %>% head(n=3)
MUT_Shannon_index %>% head (n=10)

# %%
Shannon_index_WT

# %% [markdown]
# ### Plot WT vs MUT

# %%
shannon_table <- tibble(
  index = names(MUT_Shannon_index),
  value = unlist(MUT_Shannon_index)
)

# %%
shannon_table <- shannon_table %>% 
                 mutate(DNMT3A='MUT')

# %%
shannon_table %>% head

# %%
shannon_table <- shannon_table %>% 
                 add_row(index='WT_tree_whole', value=Shannon_index_WT, DNMT3A='WT')

# %%
shannon_table %>% tail

# %%
shannon_table <- shannon_table %>% 
                 mutate(DNMT3A=as_factor(DNMT3A))

# %%
fig(7,5)
shannon_table$DNMT3A <- fct_relevel(shannon_table$DNMT3A, "WT", "MUT")

p2 <- ggplot(shannon_table, aes(x = DNMT3A, y = value, fill = DNMT3A)) +
  geom_boxplot(alpha = 0.2) +
  #geom_dotplot(binaxis = "y", stackdir = "center", binwidth=0.01, dotsize=5) +
  theme_classic(base_size = 25) +
  labs(x = "", y = "Shannon diversity index") +
  guides(fill = "none") +
  scale_fill_manual(values = c("blue", "orange")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold")) +
  coord_flip() +
  ylim(0, 4) 


p2 <- p2 + geom_dotplot(data = subset(shannon_table, DNMT3A == "WT"),
                        aes(x = DNMT3A, y = value),
                        binaxis = "y", stackdir = "center", binwidth = 0.01, dotsize = 10)

p2


# %% [markdown]
# ### Shannon DI on AX001

# %%
p0<- ggtree(AX001) +
       geom_tiplab(hjust = 0, size = 1.5)
p0

# %%
Ntip(AX001)

# %%
# the root node is usually numbered as the num_tips + 1 (the indexing first goes through all the tips and then back to the node)
num_tips <- length(AX001$tip.label)
print(num_tips)

root_node <- num_tips + 1

print(root_node)


# %%
#I confirm that 329 is indeed the root node
fig(10,12)

ggtree(AX001) + geom_cladelab(node=root_node, label="root", angle=270, #textcolor = 'blue',
                   hjust='center', offset.text= 15 , barsize= 0.5, fontsize= 5)

# %%
AX001_tips <- getTips(AX001,root_node)

# %% [markdown]
# #### generate random trees subsets of KX003 that have 33 tips

# %%
n_to_exclude <- AX001_tips %>% length() - 33
n_to_exclude

# %%
set.seed(123)
tip_combinations <- replicate(100, sample(AX001_tips, n_to_exclude), simplify = FALSE)
tip_combinations[5:7]

# %%
AX001_tree_subsets <- list()
for(i in seq_along(tip_combinations)) {
      tip_set <- tip_combinations[[i]]
      AX001_tree_subset <- drop.tip(AX001, tip_set)
      tree_name <- paste0("AX001_tree_subset_", i)
      assign(tree_name, AX001_tree_subset)
      AX001_tree_subsets[[tree_name]] <- AX001_tree_subset
}

# %% [markdown]
# ##### few checks

# %%
AX001_tree_subsets %>% head(n=3)

# %%
fig(10,6)
ggtree(AX001_tree_subset_1) + 
    geom_tiplab(color='firebrick', size=3)+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %%
fig(10,6)
ggtree(AX001_tree_subset_12) + 
    geom_tiplab(color='firebrick', size=3)+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %% [markdown]
# #### calculate Shannon DI for these randomly generated trees

# %%
nodeheights <- nodeHeights(AX001_tree_subset_25)
clades_post_cutoff=AX001_tree_subset_25$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]

ggtree(AX001_tree_subset_25) + 
    geom_text(aes(label = ifelse(isTip, node, '')), hjust = -0.3, color = 'firebrick')+
    geom_highlight(node = clades_post_cutoff,fill = "coral1", alpha=0.2) +
    geom_text2(aes(subset = node %in% clades_post_cutoff, label = node), 
             vjust = -1, color = 'blue')

# %% [markdown]
# ##### I calculate Shannon DI

# %%
AX001_clades_post_cutoff <- list()
AX001_clade_sizes <- list()
AX001_Shannon_index <- list()

for(i in seq_along(AX001_tree_subsets)) {
      tree <- AX001_tree_subsets[[i]]
      nodeheights <- nodeHeights(tree)
    
      clades_post_cutoff=tree$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]
      clades_post_cutoff_name=paste0("clades_post_cutoff_AX001_tree_subset_", i)
      AX001_clades_post_cutoff[[clades_post_cutoff_name]] <- clades_post_cutoff

      clade_sizes=sapply(clades_post_cutoff,function(node) {length(getTips(tree, node))}) 
      clade_sizes_name=paste0("clade_sizes_AX001_tree_subset_", i)
      AX001_clade_sizes[[clade_sizes_name]] <- clade_sizes

      Shannon_index <- diversity(clade_sizes, index = "shannon", MARGIN = 1, base = exp(1))
      Shannon_index_name=paste0("Shannon_index_AX001_tree_subset_", i)
      AX001_Shannon_index[[Shannon_index_name]] <- Shannon_index
}

# %% [markdown]
# #### add values to table

# %%
shannon_table <- shannon_table %>%
                      bind_rows(tibble(index = names(AX001_Shannon_index),
                                           value = unlist(AX001_Shannon_index),
                                           DNMT3A = 'AX001'))
shannon_table <- shannon_table %>% 
                 mutate(DNMT3A=as_factor(DNMT3A))

# %% [markdown]
# #### Plot

# %%
fig(7,5)
shannon_table$DNMT3A <- fct_relevel(shannon_table$DNMT3A, "WT", "MUT")

p3 <- ggplot(shannon_table, aes(x = DNMT3A, y = value, fill = DNMT3A)) +
  geom_boxplot(alpha = 0.2) +
  #geom_dotplot(binaxis = "y", stackdir = "center", binwidth=0.01, dotsize=5) +
  theme_classic(base_size = 25) +
  labs(x = "", y = "Shannon diversity index") +
  guides(fill = "none") +
  scale_fill_manual(values = c("blue", "orange", "green")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold")) +
  coord_flip() +
  ylim(0, 4) 


p3 <- p3 + geom_dotplot(data = subset(shannon_table, DNMT3A == "WT"),
                        aes(x = DNMT3A, y = value),
                        binaxis = "y", stackdir = "center", binwidth = 0.01, dotsize = 10)

p3


# %% [markdown]
# ### Shannon DI on KX003

# %%
p0<- ggtree(KX003) +
       geom_tiplab(hjust = 0, size = 1.5)
p0

# %%
Ntip(KX003)

# %%
# the root node is usually numbered as the num_tips + 1 (the indexing first goes through all the tips and then back to the node)
num_tips <- length(KX003$tip.label)
print(num_tips)

root_node <- num_tips + 1

print(root_node)


# %%
#I confirm that 329 is indeed the root node
fig(10,12)

ggtree(KX003) + geom_cladelab(node=329, label="root", angle=270, #textcolor = 'blue',
                   hjust='center', offset.text= 15 , barsize= 0.5, fontsize= 5)

# %%
KX003_tips <- getTips(KX003,329)

# %%
KX003_tips %>% length()

# %% [markdown]
# #### generate random trees subsets of KX003 that have 33 tips

# %%
KX003_tips %>% length() - 33

# %% jupyter={"outputs_hidden": true}
set.seed(123)
tip_combinations <- replicate(100, sample(KX003_tips, 295), simplify = FALSE)
tip_combinations[1:2]

# %%
KX003_tree_subsets <- list()
for(i in seq_along(tip_combinations)) {
      tip_set <- tip_combinations[[i]]
      KX003_tree_subset <- drop.tip(KX003, tip_set)
      tree_name <- paste0("KX003_tree_subset_", i)
      assign(tree_name, KX003_tree_subset)
      KX003_tree_subsets[[tree_name]] <- KX003_tree_subset
}

# %% [markdown]
# ##### few checks:

# %%
KX003_tree_subsets %>% head(n=3)

# %%
fig(10,6)
ggtree(KX003_tree_subset_1) + 
    geom_tiplab(color='firebrick', size=3)+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %%
fig(10,6)
ggtree(KX003_tree_subset_12) + 
    geom_tiplab(color='firebrick', size=3)+
    theme(plot.margin=margin(1, 2, 1, 1, "lines"))

# %% [markdown]
# #### calculate Shannon DI for these randomly generated trees

# %%
nodeheights <- nodeHeights(KX003_tree_subset_25)
clades_post_cutoff=KX003_tree_subset_25$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]

ggtree(KX003_tree_subset_25) + 
    geom_text(aes(label = ifelse(isTip, node, '')), hjust = -0.3, color = 'firebrick')+
    geom_highlight(node = clades_post_cutoff,fill = "coral1", alpha=0.2) +
    geom_text2(aes(subset = node %in% clades_post_cutoff, label = node), 
             vjust = -1, color = 'blue')

# %% [markdown]
# ##### I calculate Shannon DI

# %%
KX003_clades_post_cutoff <- list()
KX003_clade_sizes <- list()
KX003_Shannon_index <- list()

for(i in seq_along(KX003_tree_subsets)) {
      tree <- KX003_tree_subsets[[i]]
      nodeheights <- nodeHeights(tree)
    
      clades_post_cutoff=tree$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]
      clades_post_cutoff_name=paste0("clades_post_cutoff_KX003_tree_subset_", i)
      KX003_clades_post_cutoff[[clades_post_cutoff_name]] <- clades_post_cutoff

      clade_sizes=sapply(clades_post_cutoff,function(node) {length(getTips(tree, node))}) 
      clade_sizes_name=paste0("clade_sizes_KX003_tree_subset_", i)
      KX003_clade_sizes[[clade_sizes_name]] <- clade_sizes

      Shannon_index <- diversity(clade_sizes, index = "shannon", MARGIN = 1, base = exp(1))
      Shannon_index_name=paste0("Shannon_index_KX003_tree_subset_", i)
      KX003_Shannon_index[[Shannon_index_name]] <- Shannon_index
}

# %% [markdown]
# #### add values to table

# %%
shannon_table <- shannon_table %>%
                      bind_rows(tibble(index = names(KX003_Shannon_index),
                                           value = unlist(KX003_Shannon_index),
                                           DNMT3A = 'KX003'))
shannon_table <- shannon_table %>% 
                 mutate(DNMT3A=as_factor(DNMT3A))

# %% [markdown]
# #### Plot

# %%
fig(5,7)
ggplot(shannon_table, aes(x = DNMT3A, y = value, fill = DNMT3A)) +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(binaxis = "y", stackdir = "center", binwidth=0.01, dotsize=2) +
  theme_classic(base_size = 25) +
  labs(x = "", y = "Shannon diversity index") +
  guides(fill = "none") +
  scale_fill_manual(values = c("red", "blue", "green", "orange")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold"))+
  ylim(2, 4) 




# %%
fig(7,5)
shannon_table$DNMT3A <- fct_relevel(shannon_table$DNMT3A, "WT", "MUT")

p3 <- ggplot(shannon_table, aes(x = DNMT3A, y = value, fill = DNMT3A)) +
  geom_boxplot(alpha = 0.2) +
  #geom_dotplot(binaxis = "y", stackdir = "center", binwidth=0.01, dotsize=5) +
  theme_classic(base_size = 25) +
  labs(x = "", y = "Shannon diversity index") +
  guides(fill = "none") +
  scale_fill_manual(values = c("blue", "orange", "green", "red")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold")) +
  coord_flip() +
  ylim(0, 4) 


p3 <- p3 + geom_dotplot(data = subset(shannon_table, DNMT3A == "WT"),
                        aes(x = DNMT3A, y = value),
                        binaxis = "y", stackdir = "center", binwidth = 0.01, dotsize = 10)

p3


# %%
fig(6.67,7)
shannon_table$DNMT3A <- fct_relevel(shannon_table$DNMT3A, "WT", "MUT")

p4 <- ggplot(shannon_table, aes(x = value, y = DNMT3A, fill = DNMT3A)) +
  geom_boxplot(alpha = 0.2) +
  #geom_dotplot(binaxis = "y", stackdir = "center", binwidth=0.01, dotsize=5) +
  theme_classic(base_size = 25) +
  labs(x = "", y = "Shannon diversity index") +
  guides(fill = "none") +
  scale_fill_manual(values = c("blue", "red", "green", "orange")) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip() +
  xlim(0, 4) 


p4 <- p4 + geom_dotplot(data = subset(shannon_table, DNMT3A == "WT"),
                        aes(x = value, y = DNMT3A),
                        binaxis = "y", stackdir = "center", binwidth = 0.01, dotsize = 10)

p4


# %%
ggsave(paste0(figdir, today, "_shannon.pdf"), p4, dpi=300, width= 6, height=7)
