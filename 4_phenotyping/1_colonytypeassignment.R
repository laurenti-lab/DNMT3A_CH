# %% [markdown]
# # 1_CHIP10_phenotypeanalysis_colonytypeassignment
# Title:  Phenotype analysis, all experiments combined <br>
# Author: Giovanna Mantica<br>
# Description: This script analyzes flow cytometry data (exported from FlowJo 10.8.1) <br>
#

# %%
suppressPackageStartupMessages({library(tidyverse)
                               library(RColorBrewer)
                               library(gridExtra)})

# %% [markdown]
# R version: R version 4.2.3 (2023-03-15) 
# A tibble: 3 × 2
#   Package      Version
#   <chr>        <chr>  
# 1 gridExtra    2.3    
# 2 RColorBrewer 1.1-3  
# 3 tidyverse    2.0.0 

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
prefix <- format(Sys.Date(), "%Y%m%d")

# %% [markdown]
# ### import data

# %%
data <- read.csv(paste0(basedir, "/input/", "20241112_metadata_merged_cleaned.csv"))

# %%
table(data$DNMT3A,data$exp)

# %%
data$DNMT3A <- factor(data$DNMT3A, levels = c("MUT", "WT"))
data$exp <- factor(data$exp, levels = unique(data$exp))
data$TrueFalse <- factor(data$TrueFalse, levels= unique(data$TrueFalse))

# %%
unique(data$TrueFalse) #the false colony filtering has been performed at the step of genotype assignment (colonies <30 cells were defined false and excluded)

# %%
data <- data %>% 
        mutate(NeutLin.count = EarlyNeut.count + IntermediateNeut.count + MatureNeut.count) %>%
        mutate(NeutLin.freq.of.CD56neg = 100* NeutLin.count/CD56neg.count) %>%
        relocate(NeutLin.count, NeutLin.freq.of.CD56neg, .after = CD66bNeut.SSC.A)

# %% [markdown]
# ## assign colony type

# %% [markdown]
# The thresholds were established by **visual inspection of all colonies** (in FlowJo 10.8.1), in ordert to allow the simultaneous inclusion of small colonies and the exclusion of noise.
# <br>
# Here we are just visualising again that said threshold make sense also when looking at the data in its entirety.

# %%
fig(10,4)
population='Monocytes'
thresh=400

ggplot(data, aes(x = .data[[paste0(population, '.count')]])) +
  geom_histogram(binwidth = 15) +
  labs(title = paste0(population, ' count'),
       x = paste0(population, ' count')) +
  geom_vline(xintercept = thresh, color = "blue", size = 1) +
  theme_minimal() +
  #xlim(0, 500) +
  ylim(0, 20)

# %%
fig(18,6)

ggplot(data, aes(x = Monocytes.count, 
                        y = Monocytes.freq.of.NonEry, 
                        color = Monocytes.freq.of.NonEry > 2.5)) +  # Use a logical condition directly
  geom_point() +
  labs(x = "Monocytes.count",
       y = "Monocytes.freq.of.NonEry") +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"), 
                     labels = c("FALSE" = "Below Threshold", "TRUE" = "Above Threshold")) +  
  theme_minimal(base_size=16) +
  geom_hline(yintercept = 2.5, color = "blue", size = 1) +  
  geom_vline(xintercept = 150, color = "blue", size = 1) +  
  geom_vline(xintercept = 400, color = "blue", size = 1) +
  #xlim(0, 2500) +
  ylim(0, 100)

#COMMENT: we are doing this to be sure to include also smaller monocyte colonies (but exclude noise that falls in the gate)

# %%
#double check fine tuning of monocyte thresholds
overview1 <- data %>%
  group_by(exp) %>%
  summarise(
    Count_E = sum(Ery.count > 50),
    Count_M = sum(Monocytes.count > 400 | (Monocytes.count <= 400 & Monocytes.count > 150 & Monocytes.freq.of.NonEry > 2.5)),
    Count_N = sum(EarlyNeut.count > 50 | IntermediateNeut.count > 50 | MatureNeut.count > 50)
  )

print(overview1$Count_M)

overview2 <- data %>%
  group_by(exp) %>%
  summarise(
    Count_E = sum(Ery.count > 50),
    Count_M = sum(Monocytes.count > 400 | (Monocytes.count <= 400 & Monocytes.count > 150 & Monocytes.freq.of.NonEry > 5)),
    Count_N = sum(EarlyNeut.count > 50 | IntermediateNeut.count > 50 | MatureNeut.count > 50)
  )

print(overview2$Count_M)

overview3 <- data %>%
  group_by(exp) %>%
  summarise(
    Count_E = sum(Ery.count > 50),
    Count_M = sum(Monocytes.count > 500 | (Monocytes.count <= 500 & Monocytes.count > 150 & Monocytes.freq.of.NonEry > 2.5)),
    Count_N = sum(EarlyNeut.count > 50 | IntermediateNeut.count > 50 | MatureNeut.count > 50)
  )

print(overview3$Count_M)


overview5 <- data %>%
  group_by(exp) %>%
  summarise(
    Count_E = sum(Ery.count > 50),
    Count_M = sum(Monocytes.count > 500 | (Monocytes.count <= 500 & Monocytes.count > 50 & Monocytes.freq.of.NonEry > 2.5)),
    Count_N = sum(EarlyNeut.count > 50 | IntermediateNeut.count > 50 | MatureNeut.count > 50)
  )

print(overview5$Count_M)

#conclusion: fine tuning the thresholds only leads to minimal changes: these are probably minimal classification errors (either on the too stringent or on the too permissive side) that won't change the outcome

# %%
#double checking the lower threshold for small monocyte colonies
#I explore which ones are these colonies and then check visually on FlowJo that the assignment is true

filtered_M1 <- data %>%
  filter(Monocytes.count > 400 | (Monocytes.count <= 400 & Monocytes.count > 150 & Monocytes.freq.of.NonEry > 2.5))

filtered_M2 <- data %>%
  filter(Monocytes.count > 400 | (Monocytes.count <= 400 & Monocytes.count > 150 & Monocytes.freq.of.NonEry > 5))

excluded <- filtered_M1 %>% 
            anti_join(filtered_M2, by = "fcs.name") %>%
            select(exp, fcs.name, Plate_well, Monocytes.count, Monocytes.freq.of.NonEry)
excluded

# %%
#set thresholds
E <- data$Ery.count > 50
M <- data$Monocytes.count > 400 | (data$Monocytes.count <= 400 & data$Monocytes.count > 150 & data$Monocytes.freq.of.NonEry > 2.5)
N <- ((data$EarlyNeut.count > 50 | data$IntermediateNeut.count > 50 | data$MatureNeut.count > 50) &
       (!is.na(data$EarlyNeut.count) | !is.na(data$IntermediateNeut.count) | !is.na(data$MatureNeut.count)))

data2 <- data %>% mutate(Colony.type = case_when( E & !M & !N    ~ "Ery",
                                                 !E & M & !N ~ "Mono",
                                                 !E & !M & N    ~ "Neut",
                                                  E & M & !N    ~ "EryMono",
                                                  E & !M & N   ~ "EryNeut",
                                                 !E & M & N   ~ "MonoNeut",
                                                  E & M & N   ~ "EryMonoNeut",
                                                  TRUE  ~ "Other"))

data2$Colony.type <- factor(data2$Colony.type)
options(repr.matrix.max.cols=Inf)
head(data2 %>% select(exp, Plate_well,Ery.count,Monocytes.count, EarlyNeut.count, IntermediateNeut.count, MatureNeut.count, Colony.type) , n=12)
tail(data2 %>% select(exp, Plate_well,Ery.count,Monocytes.count, EarlyNeut.count, IntermediateNeut.count, MatureNeut.count, Colony.type), n=12)

levels(data2$Colony.type)

# %%
table(data2$exp, data2$Colony.type)

# %%
table(data2$exp, data2$Colony.type, data2$DNMT3A)

# %%
file_name <- paste0(basedir,"/output/",prefix, "_allexp_statsandcolonies.csv")

write.csv(data2, file = file_name, row.names = FALSE)

