# %% [markdown]
# # 2_CHIP10_phenotypeanalysis_committment
# Title:  Phenotype analysis, all experiments combined <br>
# Author: Giovanna Mantica<br>

# %%
suppressPackageStartupMessages({library(tidyverse)
                               library(RColorBrewer)
                               library(gridExtra)})

# %%
R version: R version 4.2.3 (2023-03-15) 
A tibble: 3 × 2
  Package      Version
  <chr>        <chr>  
1 gridExtra    2.3    
2 RColorBrewer 1.1-3  
3 tidyverse    2.0.0 

# %%
project <- "/.../CHIP_GM/CHIP/CHIP10"
analysis<-"phenotypeanalysis_GM"
step<-"3_phenotypeanalysis_GM"

basedir <- paste(project, analysis, step, sep='/')


# %%
options(repr.matrix.max.cols=Inf)

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %%
committment <- function(data, lineage, ylim, py){
                d <- data
                containing <- d %>% group_by(exp) %>%
                              summarise(R882 = first(R882),
                                        WT_count = sum(Colony.type %in% lineage & DNMT3A == "WT"),
                                        MUT_count = sum(Colony.type %in% lineage & DNMT3A == "MUT"),
                                        WT_Tot = sum(DNMT3A == "WT"),
                                        MUT_Tot = sum(DNMT3A == "MUT")) %>%
                              mutate(WT = 100* WT_count / WT_Tot,
                                     MUT = 100* MUT_count / MUT_Tot)
                containing_pivot <- containing %>% select(exp, R882, WT, MUT) %>% pivot_longer(cols = c(MUT, WT), names_to = "DNMT3A", values_to = "Count")
                containing_pivot$DNMT3A <- forcats::fct_relevel(containing_pivot$DNMT3A, "WT", "MUT")
                table_name= paste0(lineage, '_containing')
                assign(table_name, containing_pivot, .GlobalEnv) 
               
                d<-containing_pivot %>% filter(Count!=0)
                stats<- wilcox.test(containing$WT, containing$MUT, paired = TRUE, exact= TRUE)
                pval<-stats$p.value
    
                p <- ggplot(d, aes(x = DNMT3A, y = Count)) +
                            geom_line(aes(group = exp), linewidth = 1) +
                            geom_point(aes(shape = R882, fill= DNMT3A), stroke= 0.9, size = 5)+
                            scale_shape_manual(values = c(21, 24))+
                            scale_fill_manual(values = c(brewer.pal(3, "Set1")[2], brewer.pal(3, "Set1")[1])) +
                            theme_classic(base_size = 26) +
                            labs(x = "DNMT3A", y = "percentage of colonies", fill = NULL, shape = "R882", title = paste0(lineage, " containing colonies")) +
                            guides(fill = "none")+
                            ylim(0, ylim) + 
                            annotate("text", x = 2, y = py, label = paste("p =", round(pval, 3)), hjust = 1, size = 9)
                plot_name= paste0('p_',lineage, '_containing')
                assign(plot_name, p, .GlobalEnv)
                return(p)
                }

# %% [markdown]
# ## import data

# %%
data <- read.csv(paste0(basedir, "/output/", "20241112_allexp_statsandcolonies.csv"))

# %%
data$DNMT3A <- factor(data$DNMT3A, levels = c("MUT", "WT"))
data$exp <- factor(data$exp, levels = unique(data$exp))
data$TrueFalse <- factor(data$TrueFalse, levels= unique(data$TrueFalse))

# %% [markdown]
# ## HSC/MPP committment
# count of colonies that give rise to that lineage <br>
# Statistic performed: paired Wilcoxon test two tailed <br>
# important: HSCs/MPPs that didn't give rise to a colony have already been filtered out at the genotyping stage

# %%
Ery <- c("Ery", "EryNeut", "EryMono", "EryMonoNeut")
Mono <- c("Mono", "EryMono", "MonoNeut", "EryMonoNeut")
Neut <- c("Neut", "EryNeut", "MonoNeut", "EryMonoNeut")

# %%
committment(data,Ery, 70, 70)

# %%
committment(data,Ery, 100, 100)

# %%
committment(data,Mono, 100, 100)

# %%
committment(data, Neut, 100, 100)

# %% [markdown]
# ### export plots

# %%
lineages <- c('Ery', 'Mono', 'Neut')
for (l in lineages){
    plot <- paste0('p_',l,'_containing')
    filename <- paste0(basedir, "/output/figures/commitment/", prefix, "_committment_", l, ".pdf")
    ggsave(filename, plot = get(plot),  dpi = 300)
    }
