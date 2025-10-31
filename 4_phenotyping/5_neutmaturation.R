# %% [markdown]
# # 5_Neut maturation_plots
# Title:  Phenotype analysis, all experiments combined; statistics <br>
# Author: Giovanna Mantica, with support of the University of Cambridge Statistics Clinic (Centre for Mathematical Sciences)<br>
#

# %%
suppressPackageStartupMessages({library(tidyverse)
                               library(RColorBrewer)
                               library(gridExtra)
                                library(nlme)
                                library(lme4)
                                library(lmerTest)
                                library(openxlsx)})

# %%
R version: R version 4.2.3 (2023-03-15) 
 A tibble: 6 × 2
  Package      Version
  <chr>        <chr>  
1 gridExtra    2.3    
2 lme4         1.1-34 
3 lmerTest     3.1-3  
4 nlme         3.1-163
5 RColorBrewer 1.1-3  
6 tidyverse    2.0.0

# %%
project <- "/.../CHIP_GM/CHIP/CHIP10"
analysis<-"phenotypeanalysis_GM"
step<-"3_phenotypeanalysis_GM"

basedir <- paste(project, analysis, step, sep='/')

outdir <-paste(basedir,'output','figures','neutmaturation/', sep='/')

# %%
stat_table <- data.frame()

forest_table <- function(model, summary, effects, scale){
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
options(repr.matrix.max.cols=Inf)

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %% [markdown]
# ## Calculate % of Neut subpopulations

# %%
data2 <- data %>%
        mutate(Early.freq.of.NeutLin = (EarlyNeut.count / NeutLin.count) * 100,
        Intermediate.freq.of.NeutLin = (IntermediateNeut.count / NeutLin.count) * 100,
        Mature.freq.of.NeutLin= (MatureNeut.count / NeutLin.count) * 100) #%>%


# %%
file_name <- paste0(basedir,"/output/",prefix, "_allexp_statsandcolonies_neutdata.csv")

write.csv(data2, file = file_name, row.names = FALSE)

# %% [markdown]
# ## stacked bar plot with three tiers of Neut maturation

# %%
data_neutstat <- data2 %>%
    filter(Colony.type %in% Neut)%>%
    group_by(exp, DNMT3A) %>%
    summarise(R882=first(R882),
        Early.mean = mean(Early.freq.of.NeutLin),
        Intermediate.mean = mean(Intermediate.freq.of.NeutLin),
        Mature.mean = mean(Mature.freq.of.NeutLin),
        Early.se = sd(Early.freq.of.NeutLin) / sqrt(n()),
        Intermediate.se = sd(Intermediate.freq.of.NeutLin) / sqrt(n()),
        Mature.se = sd(Mature.freq.of.NeutLin) / sqrt(n())
  )


# %%
neutpop <- data_neutstat %>%
  pivot_longer(
    cols = -c(DNMT3A, exp, R882),
    names_to = c("maturation", "stat"),
    names_pattern = "(.*)\\.(.*)"
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  )
neutpop$maturation <- factor(neutpop$maturation, levels = c("Early", "Intermediate", "Mature"))


# %%
neutpop$maturation <- forcats::fct_relevel(neutpop$maturation, "Mature", "Intermediate", "Early")

neutpop <- neutpop %>%
  arrange(DNMT3A, desc(maturation)) %>%
  group_by(DNMT3A, exp) %>%
  mutate(mean_cumulative = cumsum(mean),
         ymax = mean_cumulative + se) %>%
  arrange(DNMT3A, maturation)
neutpop

filename <- paste0(basedir, "/output/", prefix, "_Neutpctforprism.csv")

write_csv(neutpop, filename)

# %%
fig(6,5)
neutpop$exp <- forcats::fct_relevel(neutpop$exp,'AML2', 'AML3', 'AML4', 'CHIP4')
p20 <- ggplot(neutpop, aes(x = DNMT3A, y = mean, fill = maturation)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_linerange(aes(ymin = mean_cumulative, ymax = ymax), colour = brewer.pal(11, "RdGy")[10], size = 1) +
  geom_segment(aes(x = as.numeric(DNMT3A) - 0.15, xend = as.numeric(DNMT3A) + 0.15, y = ymax, yend = ymax), colour =  brewer.pal(11, "RdGy")[10], linewidth = 1) +
  labs(y = "Mean Neutrophil composition") +
  scale_fill_brewer(palette = "Paired") +
  theme_classic(base_size=16)+
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100, 125), expand = c(0, 0.5))+
  facet_grid(. ~ exp, scales = "free_x", space = "free_x") +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16, face = "bold"))
p20  

# %% [markdown]
# ## lmer

# %% [markdown]
# #### Early

# %%
glm <- lmer(Early.freq.of.NeutLin  ~ DNMT3A + (1| exp), data = (data2 %>% filter(Colony.type %in% Neut)))
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
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=20)

# %% [markdown]
# #### Intermediate

# %%
glm <- lmer(Intermediate.freq.of.NeutLin  ~ DNMT3A + (1| exp), data = (data2 %>% filter(Colony.type %in% Neut)))
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
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=20)

# %% [markdown]
# #### Mature

# %%
glm <- lmer(Mature.freq.of.NeutLin  ~ DNMT3A + (1| exp), data = (data2 %>% filter(Colony.type %in% Neut)))
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
forest_table(glm, res, effects=c('DNMT3AMUT'), scale=20)

# %% [markdown]
# ## Multiple correction
# Flow cytometry output tests

# %%
raw_p<- c(6.15e-05, 0.000204,2.15e-08, 1.44e-06, 1.75e-05)

# %%
p.adjust(raw_p, method = "BH")
