# %% [markdown]
# # 3_CHIP10_phenotypeanalysis_HSC_differentiation output_and_index_sort
# Title:  Phenotype analysis, all experiments combined <br>
# Author: Giovanna Mantica<br>

# %%
suppressPackageStartupMessages({library(tidyverse)
                               library(RColorBrewer)
                               library(gridExtra)})

# %%
R version: R version 4.5.0 (2025-04-11) 
# A tibble: 3 × 2
  Package      Version
  <chr>        <chr>  
1 gridExtra    2.3    
2 RColorBrewer 1.1-3  
3 tidyverse    2.0.0 

# %%
project <- "/Users/giovannamantica/Library/CloudStorage/OneDrive-UniversityofCambridge/Laurentilab/CHIP_GM/CHIP/CHIP10"
analysis<-"phenotypeanalysis_GM"
step<-"3_phenotypeanalysis_GM"

basedir <- paste(project, analysis, step, sep='/')

outdir <-paste(basedir,'output','figures','differentiation/', sep='/')

# %%
options(repr.matrix.max.cols=Inf)

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %%
violin <- function(data, parameter){
                p <-ggplot(data, aes(x = DNMT3A, y = !!sym(parameter), fill = DNMT3A))+
                      geom_violin(alpha = 0.4) + 
                      # geom_dotplot(binaxis = "y", stackdir = "center",
                      #              dotsize = 0.9,
                      #              binwidth = 1,
                      #              position = position_jitter(0.03), alpha = 0.5) +
                      stat_summary(fun = median, geom = "crossbar", width = 0.8, fatten = 0.8, color = "black") +
                      theme_classic(base_size = 16) +
                      labs(x = "", y = parameter) +
                      guides(fill = "none") +
                      scale_fill_manual(values = c("blue", "red")) +
                      facet_grid(. ~ exp, scales = "free_x", space = "free_x") +
                      theme(strip.background = element_blank(),
                            strip.text = element_text(size = 16, face = "bold"))
    
                plot_name= paste0('p_',parameter)
                assign(plot_name, p, .GlobalEnv) 
    
                return(p)
                           }

# %%
violin_log <- function(data, parameter){
                p <-ggplot(data, aes(x = DNMT3A, y = !!sym(parameter), fill = DNMT3A))+
                      geom_violin(alpha = 0.4) + 
                      # geom_dotplot(binaxis = "y", stackdir = "center",
                      #              dotsize = 1.6,
                      #              binwidth = 0.015,
                      #              position = position_jitter(0.03), alpha = 0.5) +
                      stat_summary(fun = median, geom = "crossbar", width = 0.8, fatten = 0.8, color = "black") +
                      theme_classic(base_size = 16) +
                      labs(x = "", y = parameter) +
                      guides(fill = "none") +
                      scale_fill_manual(values = c("blue", "red")) +
                      facet_grid(. ~ exp, scales = "free_x", space = "free_x") +
                      scale_y_log10() +
                      theme(strip.background = element_blank(),
                            strip.text = element_text(size = 16, face = "bold"))
    
                plot_name= paste0('p_log_',parameter)
                assign(plot_name, p, .GlobalEnv) 
    
                return(p)
                           }

# %% [markdown]
# ## import data

# %%
data <- read.csv(paste0(basedir, "/output/", "20241112_allexp_statsandcolonies.csv"))

# %%
data$DNMT3A <- factor(data$DNMT3A, levels = c("WT", "MUT"))
data$exp <- factor(data$exp, levels = unique(data$exp))
data$TrueFalse <- factor(data$TrueFalse, levels= unique(data$TrueFalse))

# %%
Ery <- c("Ery", "EryNeut", "EryMono", "EryMonoNeut")
Mono <- c("Mono", "EryMono", "MonoNeut", "EryMonoNeut")
Neut <- c("Neut", "EryNeut", "MonoNeut", "EryMonoNeut")

# %% [markdown]
# ### visualise HSC output

# %%
fig(7,6)
violin_log(data, 'SingleCells.count')

# %%
table(data$DNMT3A, data$exp)

# %%
plot_name= paste0(outdir, prefix, '_log_SingleCells.count.pdf')
ggsave(plot_name, p_log_SingleCells.count, dpi = 300, width = 7, height = 4)

# %%
violin(data %>% filter(Colony.type %in% Mono), 'Monocytes.freq.of.CD56neg')

# %%
table(data[data$Colony.type == "Mono", c("DNMT3A", "exp")])

# %%
violin(data %>% filter(Colony.type %in% Neut), 'NeutLin.freq.of.CD56neg')

# %%
table(data[data$Colony.type == "Neut", c("DNMT3A", "exp")])

# %% [markdown]
# ### visualise MFI of index sorting parameters

# %%
fig(7,6)
violin_log(data, 'CD49f')

# %%
fig(7,6)
violin_log(data, 'CD34')

# %%
fig(7,6)
violin_log(data, 'CD90')

# %% [markdown]
# ### Other flow cytometry parameters - supplementary figures

# %%
fig(7,6)
violin_log(data %>% filter(Colony.type %in% Neut), 'CD66bNeut.SSC.A')

# %%
fig(7,6)
violin_log(data %>% filter(Colony.type %in% Neut), 'CD66bNeut.FSC.A')
