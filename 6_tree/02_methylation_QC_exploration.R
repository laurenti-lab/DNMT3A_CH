# %% [markdown]
# # 02. Methylation QC exploration
#
# GM On top of QC already performed by Nick Williams<br>


# %%
tidyverse <- c('dplyr', 'readr', 'ggplot2', 'tidyr', 'purrr', 'tibble', 'stringr', 'forcats')

invisible(suppressPackageStartupMessages(lapply(tidyverse, library, character.only = TRUE)))

suppressPackageStartupMessages({library(data.table)
                                library('cpgAccessData')
                                library(gridExtra)
                                library(RColorBrewer)})

# %% [markdown]
# cat("R version:", R.version$version.string, "\n")
#
# installed.packages() %>%
#   as_tibble() %>%
#   select(Package, Version) %>%
#   filter(Package %in% c('cpgAccessData', 'dplyr', 'data.table', 'RColorBrewer')) %>%
#   print(row.names = FALSE)

# %%
Sys.getenv('SINGULARITY_CONTAINER')
Sys.getenv('PYTHONUSERBASE')
Sys.getenv('SLURM_JOB_PARTITION')
Sys.getenv('SLURM_NPROCS')

# %% [markdown]
# ### Set dirs and fx

# %%
CpgConfig$DATA="/home/gm686/giovanna/CHIP28/data"

# %%
CpgConfig$METADATA="/home/gm686/giovanna/CHIP28/data/PD48544"

# %%
chrs= c(1:22)

# %%
out_dir <- ('.../figures/03/')

# %%
fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

options(repr.matrix.max.cols=Inf)

# %% [markdown]
# ### load data
#

# %% [markdown]
# #### load object

# %%
obj=cpg_get_cpgdata('PD48544',
                         chrom = chrs)


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

# %%
meth_pct_mtx <- round(meth_mtx/totreads_mtx *100)
head(meth_pct_mtx)

# %% [markdown]
# ## Coverage

# %% [markdown]
# ### Coverage distribution per colony

# %%
head(totreads_mtx)

# %%
totreads_long <- as_tibble(totreads_mtx) %>%
  rownames_to_column(var = "CpG_Site") %>%
  pivot_longer(-CpG_Site, names_to = "colony", values_to = "coverage")
tail(totreads_long)

# %%
coverage_counts <- totreads_long %>%
                   group_by(colony, coverage) %>%
                   summarise(count = n(), .groups = 'drop')
head(coverage_counts)

# %%
totreads_long %>% filter(colony=='PD48544b_lo0004') %>% filter(coverage==0) %>% dim()

# %%
total_counts_per_sample <- coverage_counts %>%
  group_by(colony) %>%
  summarise(total_CpGs = sum(count))

head(total_counts_per_sample)

# %%
coverage_pct <- coverage_counts %>%
                left_join(total_counts_per_sample, by = "colony") %>%
                mutate(pct = (count / total_CpGs) * 100)
head(coverage_pct)

# %%
fig(9,7)

p1 = ggplot(coverage_pct, aes(coverage, pct, color=colony)) + 
    geom_line() + 
    labs(title = "Coverage distribution per colony") +
    xlim(0, 100) +
    ylab('Percentage of CpGs') + 
    xlab('Reads') +
    theme_minimal(base_size=20) + theme(legend.position = "none")
p1

# %% [markdown]
# ### Coverage statistics per colony

# %%
head(totreads_mtx)

# %%
head(totreads_long)

# %%
totreads_stat_summary <- totreads_long %>%
                         group_by(colony) %>%
                         summarize(mean_coverage= mean(coverage), 
                                   median_coverage= median(coverage),
                                   tot_reads= sum(coverage))
head(totreads_stat_summary)
dim(totreads_stat_summary)

# %%
min(totreads_stat_summary$median_coverage)
min(totreads_stat_summary$tot_reads)

# %%
fig(15,5)
p3= ggplot(totreads_stat_summary, aes('', y= median_coverage)) +
    geom_violin(fill='green', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=colony), height = 0, width = 0.4, alpha = 2)+
    theme_minimal(base_size=20)+
    theme(legend.position = "none",
        axis.title.x = element_blank())

p4= ggplot(totreads_stat_summary, aes('', y= mean_coverage)) +
    geom_violin(fill='red', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=colony), height = 0, width = 0.4, alpha = 2)+
    theme_minimal(base_size=20)+
    theme(legend.position = "none",
        axis.title.x = element_blank())

p5= ggplot(totreads_stat_summary, aes('', y= tot_reads)) +
    geom_violin(fill='blue', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=colony), height = 0, width = 0.4, alpha = 2)+
    theme_minimal(base_size=20)+
    theme(legend.position = "none",
        axis.title.x = element_blank())

grid.arrange(p3, p4, p5, nrow=1)

p6 <- arrangeGrob(p3, p4, p5, nrow=1)

#ggsave(paste0(out_dir, today,'globalmethylation_allcolonies.png'), plot = p5, height=5, width=15, dpi=300)


# %% [markdown]
# #### based on mutational status
# - edit metadata matrix
# - add ASXL1 mutation information to driver column (based on tree annotation)

# %%
table(obj$meta$driver)

# %%
obj$meta %>% filter(sample=='PD48544b_lo0072') #this has to be changed to AXL1 accordingly to tree annotation (notebook 01)

# %%
obj$meta <- obj$meta %>%
  mutate(driver = ifelse(sample == 'PD48544b_lo0072', 'ASXL1', driver))

# %%
table(obj$meta$driver)

# %% [markdown]
# **2.** extract info name for each mutational status <br>

# %%
DNMT3A_col <- obj$meta %>% filter(driver == 'DNMT3A') %>% pull(sample)
length(DNMT3A_col)
JAK2_col <- obj$meta %>% filter(driver == 'JAK2:DNMT3A') %>% pull(sample)
length(JAK2_col)
TET2_col <- obj$meta %>% filter(driver == 'TET2:DNMT3A') %>% pull(sample)
length(TET2_col)
ASXL1_col <- obj$meta %>% filter(driver == 'ASXL1') %>% pull(sample)
length(ASXL1_col)
WT_col<- obj$meta %>% filter(driver == 'WT') %>% pull(sample)
length(WT_col)

# %%
obj$meta %>% filter(is.na(driver)) 

# %% [markdown]
# **3.** add mut information to summary stat matrix

# %%
totreads_stat_summary_edit <- totreads_stat_summary %>%
                             mutate(driver = case_when(colony %in% DNMT3A_col ~ 'DNMT3A',
                                                      colony %in% JAK2_col ~ 'JAK2:DNMT3A',
                                                      colony %in% TET2_col ~ 'TET2:DNMT3A',
                                                      colony %in% WT_col ~ 'WT',
                                                      colony %in% ASXL1_col ~ 'ASXL1',
                                                      TRUE ~ NA))
head(totreads_stat_summary_edit)

# %%
fig(15,8)
p3= ggplot(totreads_stat_summary_edit, aes('', y= median_coverage)) +
    geom_violin(fill='grey', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=driver), height = 0, width = 0.4, alpha = 2, size=3)+
    scale_color_brewer(palette = "Set1")+
    theme_minimal(base_size=20)+
    theme(legend.position = "bottom",
        axis.title.x = element_blank())

p4= ggplot(totreads_stat_summary_edit, aes('', y= mean_coverage)) +
    geom_violin(fill='grey', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=driver), height = 0, width = 0.4, alpha = 2, size=3)+
    scale_color_brewer(palette = "Set1")+
    theme_minimal(base_size=20)+
    theme(legend.position = "bottom",
        axis.title.x = element_blank())

p5= ggplot(totreads_stat_summary_edit, aes('', y= tot_reads)) +
    geom_violin(fill='grey', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=driver), height = 0, width = 0.4, alpha = 2, size=3)+
    scale_color_brewer(palette = "Set1")+
    theme_minimal(base_size=20)+
    theme(legend.position = "bottom",
        axis.title.x = element_blank())

grid.arrange(p3, p4, p5, nrow=1)

p7 <- arrangeGrob(p3, p4, p5, nrow=1)



# %%
fig(15,5)

p3= ggplot(totreads_stat_summary_edit %>% filter(driver %in% c('WT', 'DNMT3A')), 
           aes(x=driver, y= median_coverage)) +
    geom_violin(fill='grey', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=driver), height = 0, width = 0.4, alpha = 2, size=3)+
    scale_color_brewer(palette = "Set1")+
    theme_minimal(base_size=20)+
    theme(legend.position = "bottom",
        axis.title.x = element_blank())

p4= ggplot(totreads_stat_summary_edit %>% filter(driver %in% c('WT', 'DNMT3A')), 
           aes(x=driver, y= mean_coverage)) +
    geom_violin(fill='grey', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=driver), height = 0, width = 0.4, alpha = 2, size=3)+
    scale_color_brewer(palette = "Set1")+
    theme_minimal(base_size=20)+
    theme(legend.position = "bottom",
        axis.title.x = element_blank())

p5= ggplot(totreads_stat_summary_edit %>% filter(driver %in% c('WT', 'DNMT3A')), 
           aes(x=driver, y= tot_reads)) +
    geom_violin(fill='grey', alpha=0.2)+
    geom_boxplot(width=0.5, outlier.shape = NA)+
    geom_jitter(aes(color=driver), height = 0, width = 0.4, alpha = 2, size=3)+
    scale_color_brewer(palette = "Set1")+
    theme_minimal(base_size=20)+
    theme(legend.position = "bottom",
        axis.title.x = element_blank())

grid.arrange(p3, p4, p5, nrow=1)

p8 <- arrangeGrob(p3, p4, p5, nrow=1)

# %%
#median coverage
DNMT3A_median_coverage <- totreads_stat_summary_edit %>%
                          filter(driver=='DNMT3A') %>%
                          pull(median_coverage)
length(DNMT3A_median_coverage)

WT_median_coverage <- totreads_stat_summary_edit %>%
                          filter(driver=='WT') %>%
                          pull(median_coverage)
length(WT_median_coverage)
wtest1 <- wilcox.test(DNMT3A_median_coverage, WT_median_coverage, paired=FALSE, correct=FALSE)
print(wtest1)
pv1 <- round(wtest1$p.value, 2)
pv1

# %%
#meancoverage
DNMT3A_mean_coverage <- totreads_stat_summary_edit %>%
                          filter(driver=='DNMT3A') %>%
                          pull(mean_coverage)
length(DNMT3A_mean_coverage)

WT_mean_coverage <- totreads_stat_summary_edit %>%
                          filter(driver=='WT') %>%
                          pull(mean_coverage)
length(WT_mean_coverage)
wtest2 <- wilcox.test(DNMT3A_mean_coverage, WT_mean_coverage, paired=FALSE, correct=FALSE)
print(wtest2)
pv2 <- round(wtest2$p.value, 2)
pv2

# %%
#totreads
DNMT3A_tot_reads <- totreads_stat_summary_edit %>%
                          filter(driver=='DNMT3A') %>%
                          pull(tot_reads)
length(DNMT3A_tot_reads)

WT_tot_reads <- totreads_stat_summary_edit %>%
                          filter(driver=='WT') %>%
                          pull(tot_reads)
length(WT_tot_reads)
wtest3 <- wilcox.test(DNMT3A_tot_reads, WT_tot_reads, paired=FALSE, correct=FALSE)
print(wtest3)
pv3 <- round(wtest2$p.value, 2)
pv3

# %%
fig(15,5)

p3= p3 +
    annotate("segment", x = 1, xend = 2, y = 37, yend = 37, color = "black") +
    annotate("text", x = 1.5, y = 37 + 2, label = paste("p =", pv1), color = "black", size=6)

p4= p4 +
    annotate("segment", x = 1, xend = 2, y = 37, yend = 37, color = "black") +
    annotate("text", x = 1.5, y = 37 + 2, label = paste("p =", pv2), color = "black", size=6)

p5= p5 +
    annotate("segment", x = 1, xend = 2, y = 1.05e9, yend = 1.05e9, color = "black") +
    annotate("text", x = 1.5, y = 1.1e9 , label = paste("p =", pv3), color = "black", size=6)

grid.arrange(p3, p4, p5, nrow=1)

p8 <- arrangeGrob(p3, p4, p5, nrow=1)

# %% [markdown]
# ### Coverage statistics per CpGs

# %%
head(totreads_long)

# %%
read_threshold = 6

# %%
kept_CpG <- totreads_long %>%
            mutate(CpG_Site = as.numeric(CpG_Site)) %>%
            mutate(kept = coverage > read_threshold) %>%
            group_by(CpG_Site) %>%
            summarise(col_pct = round(mean(kept) * 100))
head(kept_CpG)

# %%
kept_CpG_pct <- kept_CpG %>% 
            as_tibble() %>%
            count(col_pct) %>% 
            mutate(CpG_pct = n / sum(n) * 100) %>%
            mutate(CpG_pct_plus1=CpG_pct+1)
head(kept_CpG_pct)

# %%
fig(8,6)
p9 <- ggplot(kept_CpG_pct, aes(x = col_pct, y = CpG_pct)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      xlab("Percentage of colonies") +
      ylab(paste0("Percentage of CpG sites with cov >", read_threshold)) +
      ggtitle("Distribution of CpG detection across colonies") +
      geom_vline(xintercept = 10, color = "red", linewidth = 1)+
      theme_minimal(base_size = 18)
p9

# %%
fig(8,6)
p10 <- ggplot(kept_CpG_pct, aes(x = col_pct, y = CpG_pct_plus1)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      xlab("Percentage of colonies") +
      ylab(paste0("Percentage of CpG sites+1 with cov >", read_threshold)) +
      ggtitle("Distribution of CpG detection across colonies") +
      geom_vline(xintercept = 10, color = "red", linewidth = 1)+
      theme_minimal(base_size = 18) +
      scale_y_log10()
p10

# %% [markdown]
# ### Global meth profile

# %%
head(meth_pct_mtx)

# %%
meth_pct_long <- meth_pct_mtx %>%
                 as_tibble() %>%
                 pivot_longer(cols=everything(), names_to='colony', values_to='meth_pct')
head(meth_pct_long)

# %%
sum(is.na(meth_pct_long$meth_pct))

# %%

p11 <- ggplot(meth_pct_long, aes(x = meth_pct)) +
      geom_histogram(binwidth = 5, fill = "yellow", color = "gray") +
      xlab("Methylation Percentage") +
      ylab("Number of CpGs") +
      ggtitle("Distribution of Methylation Percentages") +
      scale_y_continuous(expand=c(0,0))+
      theme_classic(base_size=18)
p11

# %% [markdown]
# ### individual colony meth profile

# %%

fig(9,7)

p12 <- ggplot(meth_pct_long, aes(x = meth_pct)) +
      geom_freqpoly(aes(color= colony),binwidth = 10) +
      xlab("Methylation Percentage") +
      ylab("Number of CpGs") +
      ggtitle("Distribution of Methylation Percentages") +
      theme_minimal(base_size=18) + theme(legend.position = "none")
p12

# %%
p13 <- ggplot(meth_pct_long %>% filter(colony=='PD48544b_lo0004'), aes(x = meth_pct)) +
      geom_histogram(binwidth = 5, fill = "lightblue", color = "gray") +
      xlab("Methylation Percentage") +
      ylab("Number of CpGs") +
      ggtitle("Distribution of Met Pct - PD-lo0004") +
      scale_y_continuous(expand=c(0,0))+
      theme_classic(base_size=18)
p13

# %% [markdown]
# #### colony examples (Mono/Gran)

# %%
ggplot(meth_pct_long %>% filter(colony=='PD48544b_lo0011'), aes(x = meth_pct)) +
  geom_histogram(binwidth = 5, fill = "lightblue", color = "gray") +
  xlab("Methylation Percentage") +
  ylab("Number of CpGs") +
  ggtitle("Distribution of Met Pct – lo0011") +
  scale_y_continuous(expand=c(0,0))+
  theme_classic(base_size=18)

# %% [markdown]
# #### colony examples (Ery)

# %%
obj$meta %>% as_tibble() %>% filter(Ery>=10) %>% head()

# %%
p14 <- ggplot(meth_pct_long %>% filter(colony=='PD48544b_lo0128'), aes(x = meth_pct)) +
      geom_histogram(binwidth = 5, fill = "red", color = "gray", alpha=0.2) +
      xlab("Methylation Percentage") +
      ylab("Number of CpGs") +
      ggtitle("Distribution of Methylation Percentages- Ery col PD-lo0128") +
      scale_y_continuous(expand=c(0,0))+
      theme_classic(base_size=18)
p14

# %%
fig(9,7)

p15 <- ggplot(meth_pct_long %>% filter(colony=='PD48544b_lo0128'), aes(x = meth_pct)) +
      geom_freqpoly(aes(color= colony),binwidth = 10) +
      xlab("Methylation Percentage") +
      ylab("Number of CpGs") +
      ggtitle("Distribution of Methylation Percentages - Ery col PD-lo0128") +
      theme_minimal(base_size=18) + theme(legend.position = "none")
p15

# %%
p16 <- ggplot(meth_pct_long %>% filter(colony=='PD48544b_lo0026'), aes(x = meth_pct)) +
  geom_histogram(binwidth = 5, fill = "orange", color = "gray", alpha=0.2) +
  xlab("Methylation Percentage") +
  ylab("Number of CpGs") +
  ggtitle("Distribution of Met Pct - 10 Ery lo0026") +
  scale_y_continuous(expand=c(0,0))+
  theme_classic(base_size=18)
p16

# %%
fig(9,7)

ggplot(meth_pct_long %>% filter(colony=='PD48544b_lo0026'), aes(x = meth_pct)) +
  geom_freqpoly(aes(color= colony),binwidth = 10) +
  xlab("Methylation Percentage") +
  ylab("Number of CpGs") +
  ggtitle("Distribution of Methylation Percentages") +
  theme_minimal(base_size=18) + theme(legend.position = "none")
