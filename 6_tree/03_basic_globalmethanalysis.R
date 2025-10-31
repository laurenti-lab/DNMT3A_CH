# %% [markdown]
# # 04. basic_global meth levels
#

# %%
tidyverse <- c('dplyr', 'readr', 'ggplot2', 'tidyr', 'purrr', 'tibble', 'stringr', 'forcats')

invisible(suppressPackageStartupMessages(lapply(tidyverse, library, character.only = TRUE)))

suppressPackageStartupMessages({library(data.table)
                                library('cpgAccessData')
                                library(gridExtra)
                                library(RColorBrewer)})

# %%
R version: R version 4.3.2 (2023-10-31) 
 A tibble: 4 × 2
  Package       Version
  <chr>         <chr>  
1 cpgAccessData 2.1.0  
2 data.table    1.14.10
3 dplyr         1.1.4  
4 RColorBrewer  1.1-3

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
# ### define options

# %%
CpgConfig$DATA="/home/gm686/giovanna/CHIP28/data"

# %%
CpgConfig$METADATA="/home/gm686/giovanna/CHIP28/data/PD48544"

# %%
chrs= c(1:22)

# %%
outdir <- paste0('/.../figures/04/')

# %%
fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

options(repr.matrix.max.cols=Inf)

# %% [markdown]
# ### load data

# %% [markdown]
# #### load object

# %%
obj=cpg_get_cpgdata('PD48544',
                         chrom = chrs)

# %%
names(obj)

# %%
head(obj$details)

# %%
dim(obj$meta)

# %% [markdown]
# #### create methylation count matrices

# %%
meth_mtx <- obj$MF + obj$MR
unmeth_mtx <- obj$UF + obj$UR
totreads_mtx <- meth_mtx + unmeth_mtx

# %%
dim(meth_mtx) 
dim(unmeth_mtx)
dim(totreads_mtx)

# %% [markdown]
# #### create methylation pct matrices

# %%
meth_pct_mtx <- round(meth_mtx/totreads_mtx *100)
head(meth_pct_mtx)

# %% [markdown]
# ### edit metadata matrix

# %% [markdown]
# #### include ASXL1 as driver

# %%
head(obj$meta)

# %%
table(obj$meta$driver)

# %%
obj$meta %>% filter(sample=='PD48544b_lo0072') #this has to be changed to AXL1 accordingly to tree annotation (notebook 01)

# %%
obj$meta <- obj$meta %>%
  mutate(driver = ifelse(sample == 'PD48544b_lo0072', 'ASXL1', driver))

# %% [markdown]
# #### identify Ery colony and high read colony to filter them out
# (accordingly to QC checks, notebook 02)

# %%
obj$meta %>% filter(Ery > 30) %>% select(sample, PD_id, driver)

# %%
obj$meta %>% filter(sample=='PD48544b_lo0119') %>% select(sample, PD_id, driver)

# %% [markdown]
# #### generate colony driver vectors + filtering Ery + filtering abnormal read colony

# %%
DNMT3A_col <- obj$meta %>% filter(driver == 'DNMT3A') %>% 
                           filter(!(sample=='PD48544b_lo0119')) %>%
                           pull(sample)
length(DNMT3A_col)

# %%
JAK2_col <- obj$meta %>% filter(driver == 'JAK2:DNMT3A') %>% pull(sample)
length(JAK2_col)

# %%
TET2_col <- obj$meta %>% filter(driver == 'TET2:DNMT3A') %>% pull(sample)
length(TET2_col)

# %%
ASXL1_col <- obj$meta %>% filter(driver == 'ASXL1') %>% pull(sample)
length(ASXL1_col)

# %%
WT_col<- obj$meta %>% filter(!(Ery > 30)) %>% 
                      filter(driver == 'WT') %>%
                      pull(sample)
length(WT_col)

# %% [markdown]
# #### generate cell type vectors

# %%
HSC_col <- obj$meta %>% filter(Cell_seeded == 'HSC') %>% pull(sample)
length(HSC_col)

# %%
GMP_col <- obj$meta %>% filter(Cell_seeded == 'GMP') %>% pull(sample)
length(GMP_col)

# %% [markdown]
# ## global mean methylation level

# %% [markdown]
# ### per promoter vs not

# %%
meth_pct_long <- meth_pct_mtx %>%
                 as_tibble() %>%
                 pivot_longer(cols=everything(), names_to='colony', values_to='meth_pct')
head(meth_pct_long)

# %%
meth_pct_stats <- meth_pct_mtx %>%
                 as_tibble() %>%
                 pivot_longer(cols=everything(), names_to='colony', values_to='meth_pct') %>%
                 group_by(colony) %>%
                 summarize(mean_meth_pct =mean(meth_pct, na.rm = TRUE))
head(meth_pct_stats)

# %%
prom_pct_mtx <- meth_pct_mtx[obj$details$promoters, ]
dim(prom_pct_mtx)

# %%
prom_pct_stats <- prom_pct_mtx %>%
                 as_tibble() %>%
                 pivot_longer(cols=everything(), names_to='colony', values_to='meth_pct') %>%
                 group_by(colony) %>%
                 summarize(mean_meth_pct =mean(meth_pct, na.rm = TRUE)) %>%
                 mutate(promoters='YES')
head(prom_pct_stats)
dim(prom_pct_stats)

# %%
nonprom_pct_mtx <- meth_pct_mtx[!(obj$details$promoters), ]
dim(nonprom_pct_mtx)

# %%
nonprom_pct_stats <- nonprom_pct_mtx %>%
                     as_tibble() %>%
                     pivot_longer(cols=everything(), names_to='colony', values_to='meth_pct') %>%
                     group_by(colony) %>%
                     summarize(mean_meth_pct =mean(meth_pct, na.rm = TRUE)) %>%
                     mutate(promoters='NO')
head(nonprom_pct_stats)
dim(nonprom_pct_stats)

# %%
prom_type_stats <- bind_rows(prom_pct_stats, nonprom_pct_stats)
prom_type_stats$promoters<-as_factor(prom_type_stats$promoters)
head(prom_type_stats)
dim(prom_type_stats)

# %%
fig(5,7)

p1 <- ggplot(prom_type_stats, aes(x=promoters, mean_meth_pct, fill= promoters )) + 
    #geom_violin(alpha = 0.4) + 
    geom_boxplot(fill='white', width = 0.7) + 
    geom_jitter(height = 0, width = 0.3, alpha = 0.6 , color='black') + 
    ylab('Mean methylation (%)') + 
    xlab('promoter') + 
    ylim(0,100)+
    theme_minimal(base_size=18) +
    theme(legend.position='none')+
    ggtitle("Mean methylation at promoters") 
    

p1

# %% [markdown]
# ### per promoter per CpG context

# %%
table(obj$details$cpg_context, obj$details$promoters)

# %%
contexts <- unique(obj$details$cpg_context)

for (type in contexts) {
    prom_type <- obj$details$promoters & (obj$details$cpg_context == type)
    prom_type_pct_mtx <- meth_pct_mtx[prom_type, ]
    prom_type_pct_stats <- prom_type_pct_mtx %>%
                 as_tibble() %>%
                 pivot_longer(cols=everything(), names_to='colony', values_to='meth_pct') %>%
                 group_by(colony) %>%
                 summarize(mean_meth_pct =mean(meth_pct, na.rm = TRUE)) %>%
                 mutate(promoter='YES') %>%
                 mutate(cpg_context=type)
    
    assign(paste("prom_", type, sep = ""), prom_type)
    assign(paste("prom_", type, "_pct_mtx", sep = ""), prom_type_pct_mtx)
    assign(paste("prom_", type, "_pct_stats", sep = ""), prom_type_pct_stats)
    
}

# %%
prom_cpg_inter_pct_mtx %>% dim()
prom_cpg_islands_pct_mtx %>% dim()
prom_cpg_shelves_pct_mtx %>% dim()
prom_cpg_shores_pct_mtx %>% dim()

# %%
prom_cpg_inter_pct_stats %>% dim()
prom_cpg_islands_pct_stats %>% dim()
prom_cpg_shelves_pct_stats %>% dim()
prom_cpg_shores_pct_stats %>% dim()

# %%
prom_cpg_context_stats <- bind_rows(prom_cpg_inter_pct_stats,
                                   prom_cpg_islands_pct_stats,
                                   prom_cpg_shelves_pct_stats,
                                   prom_cpg_shores_pct_stats)

prom_cpg_context_stats$cpg_context <- as_factor(prom_cpg_context_stats$cpg_context)
head(prom_cpg_context_stats)
dim(prom_cpg_context_stats)

# %%
unique(prom_cpg_context_stats$promoter)

# %%
levels(prom_cpg_context_stats$cpg_context)

# %%
fig(5,7)
p2 <- ggplot(prom_cpg_context_stats, aes(x=cpg_context, y= mean_meth_pct, fill= cpg_context )) + 
    #geom_violin(alpha = 0.4) + 
    geom_boxplot(fill='white', width = 0.7) + 
    geom_jitter(height = 0, width = 0.3, alpha = 0.6 , color='black') + 
    ylab('Mean methylation (%)') + 
    xlab('promoter') + 
    ylim(0,100)+
    theme_minimal(base_size=18) +
    theme(legend.position='none', axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle("Mean methylation at promoters") 
    

p2

# %%
table(obj$details$cpg_context, obj$details$promoters)

# %%
contexts <- unique(obj$details$cpg_context)

for (type in contexts) {
    nonprom_type <- !(obj$details$promoters) & (obj$details$cpg_context == type)
    nonprom_type_pct_mtx <- meth_pct_mtx[nonprom_type, ]
    nonprom_type_pct_stats <- nonprom_type_pct_mtx %>%
                 as_tibble() %>%
                 pivot_longer(cols=everything(), names_to='colony', values_to='meth_pct') %>%
                 group_by(colony) %>%
                 summarize(mean_meth_pct =mean(meth_pct, na.rm = TRUE)) %>%
                 mutate(promoter='NO') %>%
                 mutate(cpg_context=type)
    
    assign(paste("nonprom_", type, sep = ""), nonprom_type)
    assign(paste("nonprom_", type, "_pct_mtx", sep = ""), nonprom_type_pct_mtx)
    assign(paste("nonprom_", type, "_pct_stats", sep = ""), nonprom_type_pct_stats)
    
}

# %%
nonprom_cpg_shores_pct_stats %>% head()

# %%
cpg_context_stats <- bind_rows(prom_cpg_inter_pct_stats,
                                   prom_cpg_islands_pct_stats,
                                   prom_cpg_shelves_pct_stats,
                                   prom_cpg_shores_pct_stats,
                                   nonprom_cpg_inter_pct_stats,
                                   nonprom_cpg_islands_pct_stats,
                                   nonprom_cpg_shelves_pct_stats,
                                   nonprom_cpg_shores_pct_stats)

cpg_context_stats$cpg_context <- forcats::fct_relevel(cpg_context_stats$cpg_context, 'cpg_islands','cpg_shores', 'cpg_shelves', 'cpg_inter')
cpg_context_stats$promoter <- as_factor(cpg_context_stats$promoter)
head(cpg_context_stats)
tail(cpg_context_stats)
dim(cpg_context_stats)

# %%
fig(10,8)
p3 <- ggplot(cpg_context_stats, aes(x=cpg_context, y= mean_meth_pct, fill= cpg_context)) + 
    #geom_violin(alpha = 0.4) + 
    geom_boxplot(fill='white', width = 0.7) + 
    geom_jitter(height = 0, width = 0.3, alpha = 0.6 , color='black') + 
    ylab('Mean methylation (%)') + 
    xlab('promoter') + 
    ylim(0,100)+
    theme_minimal(base_size=18) +
    theme(legend.position='none', axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle("Mean methylation level")  +
    facet_grid(. ~ promoter, scales = "free_x", space = "free_x") +
    theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold"))
    

p3

# %%
fig(10,8)
p4 <- ggplot(cpg_context_stats, aes(x=promoter, y= mean_meth_pct)) + 
    #geom_violin(alpha = 0.4) + 
    geom_boxplot(fill='white', width = 0.7) + 
    geom_jitter(height = 0, width = 0.3, alpha = 0.6 , color='black') + 
    ylab('Mean methylation (%)') + 
    xlab('promoter') + 
    ylim(0,100)+
    theme_minimal(base_size=18) +
    theme(legend.position='none', axis.text.x = element_text(angle = 45, hjust = 1))+
    ggtitle("Mean methylation level")  +
    facet_grid(. ~ cpg_context, scales = "free_x", space = "free_x") +
    theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold"))
    

p4

# %% [markdown]
# ### per driver

# %%
meth_pct_stats <- meth_pct_stats %>%
                    mutate(driver=if_else(colony %in% WT_col, 'WT', NA)) %>%
                    mutate(driver=if_else(colony %in% DNMT3A_col, 'DNMT3A', driver)) %>%
                    mutate(driver=if_else(colony %in% ASXL1_col, 'ASXL1', driver)) %>%
                    mutate(driver=if_else(colony %in% TET2_col, 'DNMT3A_TET2', driver)) %>%
                    mutate(driver=if_else(colony %in% JAK2_col, 'DNMT3A_JAK2', driver))
head(meth_pct_stats)


# %%
fig(5,7)
p5 = ggplot(meth_pct_stats %>% filter(!is.na(driver)), 
            aes(x= driver, y=mean_meth_pct, fill= driver)) + 
        geom_violin(alpha = 0.4) + 
        geom_boxplot(fill='white', width = 0.2) + 
        geom_jitter(height = 0, width = 0.3, alpha = 0.6 , color='black') + 
        ylab('Mean methylation level (%)') + 
        xlab('') + 
        theme_classic(base_size=18) +
        theme(legend.position='none', axis.text.x = element_text(angle = 45, hjust = 1))+
        scale_fill_brewer(palette="Set1") +
        ylim(75,85)

p5

# %%
DNMT3A_meth_rate <- meth_pct_stats %>%
                    filter(driver=='DNMT3A') %>%
                    pull(mean_meth_pct)
length(DNMT3A_meth_rate)

WT_meth_rate <- meth_pct_stats %>%
                    filter(driver=='WT') %>%
                    pull(mean_meth_pct)
length(WT_meth_rate)

wtest5 <- wilcox.test(DNMT3A_meth_rate, WT_meth_rate, paired=FALSE, correct=FALSE)
print(wtest5)
pv5 <- round(wtest5$p.value, 15)
pv5


# %%
wtest6 <- wilcox.test(DNMT3A_meth_rate, JAK2_meth_rate, paired=FALSE, correct=FALSE)
print(wtest6)
pv6 <- round(wtest6$p.value, 7)
pv6

# %% [markdown]
# #### summary plot with stats

# %%
p5<- p5 + 
        annotate("segment", x = 1, xend = 2, y = 84, yend = 84, color = "black") +
        annotate("text", x = 1.5, y = 84 + 0.4, label = paste("p =", pv5), color = "black", size=6) +
        annotate("segment", x = 2, xend = 3, y = 83, yend = 83, color = "black") +
        annotate("text", x = 2.5, y = 83 + 0.4, label = paste("p =", pv6), color = "black", size=6) 
p5
