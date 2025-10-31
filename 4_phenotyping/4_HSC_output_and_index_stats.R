# %% [markdown]
# # 4_CHIP10_phenotypeanalysis_stats
# Title:  Phenotype analysis, all experiments combined<br>
# Author: Giovanna Mantica, with support of the University of Cambridge Statistics Clinic (Centre for Mathematical Sciences) <br>

# %%
suppressPackageStartupMessages({library(tidyverse)
                               library(RColorBrewer)
                               library(gridExtra)
                                library(nlme)
                                library(lme4)
                                library(lmerTest)
                                library(openxlsx)})

# %% [markdown]
# R version: R version 4.5.0 (2025-04-11) 
#  A tibble: 6 × 2
#   Package      Version
#   <chr>        <chr>  
# 1 gridExtra    2.3    
# 2 lme4         1.1-37 
# 3 lmerTest     3.1-3  
# 4 nlme         3.1-168
# 5 RColorBrewer 1.1-3  
# 6 tidyverse    2.0.0  

# %%
project <- "/.../CHIP_GM/CHIP/CHIP10"
analysis<-"phenotypeanalysis_GM"
step<-"3_phenotypeanalysis_GM"

basedir <- paste(project, analysis, step, sep='/')

outdir <-paste(basedir,'output','figures','stats/', sep='/')

# %%
options(repr.matrix.max.cols=Inf)

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %%
stat_table <- data.frame()

forest_table <- function(model, summary, effects, scale=15){
                CIs <- confint(model, level=0.95)
                df <- data.frame(details= as.character(summary$call)[2],
                         var=rownames(summary$coefficients),
                         mid = summary$coefficients[, "Estimate"],
                         lo = CIs[rownames(summary$coefficients), "2.5 %"],
                         hi = CIs[rownames(summary$coefficients), "97.5 %"],
                         pval= summary$coefficients[rownames(summary$coefficients),"Pr(>|t|)"],
                         st.er= summary$coefficients[rownames(summary$coefficients),"Std. Error"])
                df$sig <- df$pval<0.05
                stat_table <<- rbind(stat_table, df)
    
                p <- ggplot(df %>% filter(var %in% effects), aes(x = var, y = mid, ymin = lo, ymax = hi)) +
                     geom_errorbar(aes(color=sig), width=0.2, linewidth=1.5) +                     
                     geom_point(aes(fill=sig, color=sig), size=3, pch=21, stroke=0.8) +
                     geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth=0.5) + 
                     scale_color_manual(values=c("black", "red2")) +
                     scale_fill_manual(values=c("black", "red2")) +
                     coord_flip() +                         
                     labs(
                       title = df$details, x = "", y = "") +
                     theme_classic(base_size=20)+
                     theme(legend.position = "none",
                            plot.title = element_text(size = 10))+
                     ylim(-scale,scale)
    
                
                parameter=strsplit(df$details, "~")[[1]][1]
                plot_name= paste0('p_',parameter)
                assign(plot_name, p, .GlobalEnv)             

                    
                return(p)
                
                         }


# %% [markdown]
# ## import data

# %%
data <- read.csv(paste0(basedir, "/output/", "20241101_allexp_statsandcolonies.csv"))

# %%
data$DNMT3A <- factor(data$DNMT3A, levels = c("WT", "MUT"))
data$exp <- factor(data$exp, levels = unique(data$exp))
data$TrueFalse <- factor(data$TrueFalse, levels= unique(data$TrueFalse))

# %%
data %>% head()

# %%
data$DNMT3A %>% unique()

# %%
Ery <- c("Ery", "EryNeut", "EryMono", "EryMonoNeut")
Mono <- c("Mono", "EryMono", "MonoNeut", "EryMonoNeut")
Neut <- c("Neut", "EryNeut", "MonoNeut", "EryMonoNeut")

# %%
data<- data %>% mutate(R882H = case_when(DNMT3A == 'MUT' & R882 == 'H' ~1, TRUE ~0)) %>%
                relocate(R882H, .before=Plate_well)

# %%
table(data$Sex, data$Age)

# %% [markdown]
# # lmer

# %% [markdown]
# ## Size

# %% [markdown]
# #### DNMT3A + exp

# %%
glm <- lmer(log(SingleCells.count)  ~ DNMT3A + (1| exp), data = data)
res <- summary(glm)
res


# %% [markdown]
# ##### model check

# %%
fig(5,5)
qqnorm(residuals(glm), main = "Q-Q Plot of Residuals")
qqline(residuals(glm), col = "red")

# %%
hist(residuals(glm))

# %% [markdown]
# ##### figure

# %%
fig(6,4)
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=1)

# %% [markdown]
# ## Monocytes %

# %% [markdown]
# #### DNMT3A + exp

# %%
glm <- lmer(Monocytes.freq.of.CD56neg  ~ DNMT3A + (1| exp), data = (data %>% filter(Colony.type %in% Mono)))
res <- summary(glm)
res


# %%
confint(glm, level=0.95)

# %% [markdown]
# ##### model check

# %%
fig(5,5)
qqnorm(residuals(glm), main = "Q-Q Plot of Residuals")
qqline(residuals(glm), col = "red")

# %%
hist(residuals(glm))

# %% [markdown]
# ##### figure

# %%
fig(6,4)
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=10)

# %% [markdown]
# ## Neut %

# %% [markdown]
# #### DNMT3A + exp

# %%
glm <- lmer(NeutLin.freq.of.CD56neg  ~ DNMT3A + (1| exp), data = (data %>% filter(Colony.type %in% Neut)))
res <- summary(glm)
res


# %%
confint(glm, level=0.95)

# %% [markdown]
# ##### model check

# %%
fig(5,5)
qqnorm(residuals(glm), main = "Q-Q Plot of Residuals")
qqline(residuals(glm), col = "red")

# %%
hist(residuals(glm))

# %% [markdown]
# ##### figure

# %%
fig(6,4)
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=10)

# %% [markdown]
# ## Index sorting

# %% [markdown]
# ### CD49f

# %% [markdown]
# #### DNMT3A + exp

# %%
table(data$DNMT3A, data$exp)

# %%
glm <- lmer(log(CD49f)  ~ DNMT3A + (1| exp), data = data)
res <- summary(glm)
res

# %%
fig(6,4)
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=0.5)

# %% [markdown]
# ### CD34

# %% [markdown]
# #### DNMT3A + exp

# %%
glm <- lmer(log(CD34)  ~ DNMT3A + (1| exp), data = data)
res <- summary(glm)
res

# %% [markdown]
# ##### model check

# %%
fig(5,5)
qqnorm(residuals(glm), main = "Q-Q Plot of Residuals")
qqline(residuals(glm), col = "red")

# %%
hist(residuals(glm))

# %% [markdown]
# ##### figure

# %%
fig(6,4)
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=0.25)

# %% [markdown]
# ### CD90

# %% [markdown]
# #### DNMT3A + exp

# %%
glm <- lmer(log(CD90)  ~ DNMT3A + (1| exp), data = data %>% filter(CD90>0))
res <- summary(glm)
res

# %% [markdown]
# ##### model check

# %%
fig(5,5)
qqnorm(residuals(glm), main = "Q-Q Plot of Residuals")
qqline(residuals(glm), col = "red")

# %%
hist(residuals(glm))

# %% [markdown]
# ##### figure

# %%
fig(6,4)
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=0.25)

# %% [markdown]
# ## Other flow cytometry parameters for supplementary

# %% [markdown]
# ### FSC-A CD66b neutrophils

# %% [markdown]
# #### DNMT3A + exp

# %%
glm <- lmer(log(CD66bNeut.FSC.A)  ~ DNMT3A + (1| exp), data = data)
res <- summary(glm)
res


# %% [markdown]
# ##### model check

# %%
fig(5,5)
qqnorm(residuals(glm), main = "Q-Q Plot of Residuals")
qqline(residuals(glm), col = "red")

# %%
hist(residuals(glm))

# %% [markdown]
# ##### figure

# %%
fig(6,4)
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=0.4)

# %% [markdown]
# ### SSC-A CD66b neutrophils

# %% [markdown]
# #### DNMT3A + exp

# %%
glm <- lmer(log(CD66bNeut.SSC.A)  ~ DNMT3A + (1| exp), data = data)
res <- summary(glm)
res


# %% [markdown]
# ##### model check

# %%
fig(5,5)
qqnorm(residuals(glm), main = "Q-Q Plot of Residuals")
qqline(residuals(glm), col = "red")

# %%
hist(residuals(glm))

# %% [markdown]
# ## p value Multiple correction

# %% [markdown]
# ### Phenotyping: Monocytes and Neuts and FSC-A and SSC-A

# %%
raw_p<- c(6.15e-05, 0.000204, 0.452, 0.004)

# %%
p.adjust(raw_p, method = "BH")

# %% [markdown]
# ### Index sorting

# %%
raw_p<-c(0.006, 0.341, 0.024)

# %%
p.adjust(raw_p, method = "BH")
