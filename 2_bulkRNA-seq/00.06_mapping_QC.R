# %% [markdown]
# # 00.06_mapping QC
#

# %%
suppressPackageStartupMessages({library(tidyverse)
                                library(RColorBrewer)
                                library(gridExtra)
                                })

# %% [markdown]
# R version: R version 4.3.2 (2023-10-31) 
# A tibble: 4 × 2
#   Package      Version
#   <chr>        <chr>  
# 1 biomaRt      2.58.2 
# 2 gridExtra    2.3    
# 3 RColorBrewer 1.1-3  
# 4 tidyverse    2.0.0 

# %%
basedir <- paste(laurenti, user, project, sep='/')
basename <- '_CHIP30_'
bam_dir <- '/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/processed/star_bam'
outdir<-paste(basedir, 'output', sep='/') 
figdir<- paste(datadir, 'figures', '06.01/', sep='/')

# %%
options(repr.matrix.max.cols=20) #Inf

fig <- function(width, height) {
    options(repr.plot.width = width, repr.plot.height = height)
}

# %% [markdown]
# ### MappingQC

# %%
setwd(bam_dir)
list.files<-list.files(path = '.', pattern="Log.final.out")
list.files<-unique(list.files)
file.names<-sub("Log\\.final\\.out$", "", list.files)

# %%
list.files %>% length()

# %%
num.mapped<-data.frame()
for (i in (1:(length(file.names)))){
  dat <- read.csv(list.files[i], header = FALSE, skip = 1) 
  dat1 <- dat[7,] %>% separate(V1, into=c('mapped','uq_map'), sep='\\|') %>% select(uq_map)
  dat2<- dat[8,] %>% separate(V1, into=c('mapped','uq_map_perc'), sep='\\|') %>% select(uq_map_perc)
  dat2$uq_map_perc <- as.numeric(substr(dat2$uq_map_perc, 1, nchar(dat2$uq_map_perc) - 1))
  dat<-cbind(dat1,dat2)
  dat<-data.frame(cut=file.names[i], uq_mapped_reads=dat$uq_map, uq_mapped_perc=dat$uq_map_perc)
  num.mapped<-rbind(num.mapped, dat)
}
colnames(num.mapped)<-c("sample_id", "uq_mapped_reads", "uq_mapped_perc")

# %%
hist(num.mapped$uq_mapped_perc, main = "Histogram of uniquely_mapped_percentages", xlab = "Mapped Percentage", col = "darkgray", border="white")


# %%
#Plot distribution of uniquely mapped reads
num.mapped$uq_mapped_reads<-as.numeric(num.mapped$uq_mapped_reads)
ggplot(num.mapped, aes(x=uq_mapped_reads)) +
          geom_histogram(binwidth=1000000) + 
          scale_y_continuous(breaks = seq(0, 65, 10)) + 
          scale_x_continuous(breaks = seq(0, 60000000, 1000000)) + 
          theme_classic(base_size=16)+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  + 
          labs(x = "Number of uniquely mapped reads", y="Samples")

# %%
#Plot fractions from mapping output

map.reads<-data.frame()
for (i in (1:(length(file.names)))){
  dat <- read.csv(list.files[i], header = FALSE, skip = 1) 
  dat <- dat[c(8,23,25, 30, 32),] %>% separate(V1, into=c('legend','values'), sep='\\|') %>% select(-V2)
  dat$legend<-factor(c( "Uniq. mapped", "Multi-mapping", "Too many loci", "Unmapped - short", "Unmapped - other"), levels = c("Uniq. mapped", "Multi-mapping", "Too many loci", "Unmapped - short", "Unmapped - other"))
  dat$values<-str_sub(dat$values, end = -2)
  dat$mapping<-c("unique", "not unique",  "not unique","not unique","not unique")
  dat$cut<-c(rep(file.names[i], 5))
  map.reads<-rbind(map.reads, dat)
}

map.reads$values<-as.numeric(map.reads$values)

map.reads$legend <- factor(map.reads$legend, 
                           levels = c("Unmapped - other", "Unmapped - short", "Too many loci", "Multi-mapping", "Uniq. mapped"))

# %%
ggplot(map.reads, aes(x = cut, y = values, fill = legend)) + 
  geom_bar(stat = "identity") + 
  theme_minimal() + 
  scale_fill_manual(values = c("pink3", "palevioletred2", 'gray', "lightpink2", "indianred", "mediumseagreen")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 13), 
        axis.title = element_text(size = 0.5))

# %%
#this useful if you have many samples per patient/mouse, not applicable here
# box<-map.reads %>% filter(mapping=="unique")
# ggplot(box, aes(x = patient, y = values)) +
#   geom_boxplot() +
#   geom_jitter(width = 0.2, size=3) +  # To display individual data points
#   theme_minimal() +
#   labs(x = "patient", y = "Values") +  # Label axes and legend
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
#         axis.title = element_text(size = 12), 
#         legend.position = "none")

# %% [markdown]
# ### Feature counts QC

# %%
setwd('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/output/tables')

# %%
list.files<-list.files(path = '.', pattern="summary")

fc.output<-NULL
for (i in (1:(length(list.files)))){
  dat <- read.csv(list.files[i], header = TRUE, sep = '\t') 
  if (is.null(fc.output)) {
    fc.output <- dat
  } else {
    fc.output <- cbind(fc.output, dat)
  }
}

# %%
fc.output <- read.csv('202503_CHIP30_countmatrix_v2.txt.summary', header = TRUE, sep = '\t') 

# %%
colnames(fc.output) <- gsub("Aligned\\.sortedByCoord\\.out\\.bam|\\.\\.\\/star_bam\\.|\\.\\.star_bam\\.", "", colnames(fc.output))
status<-fc.output$Status
rownames(fc.output) <- fc.output$Status
fc.output <- fc.output[, -1]  

# %%
fc.output

# %%
fc.output = as.data.frame(sapply(fc.output, as.numeric))
total<-colSums(fc.output)
rownames(fc.output)<-status
fc.output<-fc.output[as.logical(rowSums(fc.output != 0)), ]
status<-rownames(fc.output)

assigned.total<-fc.output[1,]

# %%
fc.output <- (mapply('/', fc.output, total))*100
rownames(fc.output)<-status

# %%
fc.output

# %%
fc.output.t<-t(fc.output)
fc.output.t<-as.data.frame(fc.output.t)
fc.output.t$cut<-rownames(fc.output.t)

fc.output.l<-fc.output.t %>% pivot_longer(!cut, names_to = "status", values_to = "count")
fc.output.l<- fc.output.l %>% arrange(cut)

fc.output.l$status<-factor(fc.output.l$status, levels = c("Unassigned_Singleton", "Unassigned_NoFeatures", "Unassigned_MultiMapping", "Unassigned_Ambiguity", "Assigned"))



# %%
fig(9,8)
ggplot(fc.output.l, aes(x = cut, y = count, fill = status)) + 
  geom_bar(stat = "identity") +
  theme_minimal(base_size=16)+ 
  scale_fill_manual(values = c("red", "palevioletred2", "lightpink2", "indianred", "mediumseagreen")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 13), axis.title=element_text(size=0.1))+
  ggtitle('v2 script, exon')

# %% [markdown]
# ### QC parameters

# %%
setwd('/rds/project/rds-0p1wRpSNFIw/users/gm686/CHIP30/output/tables')

# %%
v2 <- read.csv('202503_CHIP30_countmatrix_v2.txt',header = FALSE, sep = '\t',comment.char = "#")
names(v2)<- v2[1,]
v2<-v2[-1,]

# %%
names(v2)<-gsub("Aligned.sortedByCoord.out.bam", "",names(v2))


# %%
v2_cl <- v2 %>% select(- c(Chr,Start,End,Strand))

# %%
write.csv(v2_cl, file='202503_CHIP30_countmatrix_v2_clean.txt', row.names=FALSE)

# %%
counts.df<-v2_cl
rownames(counts.df)<-counts.df$Geneid
counts.df<-counts.df[,-1]
counts.df<-counts.df[,-1]

# %%
counts.df %>% head(2)

# %%
MT_genes.vec <- rownames(counts.df)[grepl("^MT-", rownames(counts.df))]
RB_genes.vec <- rownames(counts.df)[grepl("^RPS|^RPL", rownames(counts.df))]

# %%
#Get nCount without mito genes WHY
dat <- counts.df[!rownames(counts.df) %in% MT_genes.vec, ]
dat[] <- lapply(dat, as.numeric)
nCount.dat<-colSums(dat)

# %%

# %%
#Get nCount, nFeatures, nFeatures (>5 reads), %MT genes in a table

counts.df<- mutate_all(counts.df, function(x) as.numeric(as.character(x)))

#nCount
nCount<-colSums(counts.df)

#nFeatures
nFeatures<-c()
for (i in 1:(ncol(counts.df))){
  dat<-counts.df[,i]
  nFeatures<-c(nFeatures,length(dat[dat>0]))
}

#nFeatures expressed >5                       
nFeatures.filt<-c()
for (i in 1:(ncol(counts.df))){
  dat<-counts.df[,i]
  nFeatures.filt<-c(nFeatures.filt,length(dat[dat>5])) 
}
nFeatures.difference<-nFeatures-nFeatures.filt

#MT
mt.counts<- counts.df%>% filter(rownames(counts.df) %in% MT_genes.vec)
mt.Total<-colSums(mt.counts)
mt.perc<-(mt.Total/nCount)*100

#RB
rb.counts<- counts.df%>% filter(rownames(counts.df) %in% RB_genes.vec)
rb.Total<-colSums(rb.counts)
rb.perc<-(rb.Total/nCount)*100      
                       
                       
feat.count.df<-data.frame(sample_ID=colnames(counts.df), lib_size=nCount, feature_counts=nFeatures, feature_counts_overfive=nFeatures.filt, feature_counts_difference= nFeatures.difference, mito_fraction=mt.perc, ribo_fraction=rb.perc)

feat.count.df<-feat.count.df %>% mutate(pass_QC_standard = ifelse(feature_counts_overfive > 5000 & 
                             lib_size > 500000 & 
                             mito_fraction < 50, "Y", "N"), pass_QC_strict = ifelse(feature_counts_overfive > 7000 & 
                             lib_size > 500000 & 
                             mito_fraction < 15, "Y", "N"))

# %%
feat.count.df

# %%
fig(3.5,5)
ggplot(feat.count.df, aes(x = "", y = mito_fraction)) + 
  geom_violin(fill = "lightpink2", alpha = 0.5) +  
  geom_jitter(width = 0.1, size = 3, color='gray33' ,alpha = 1) +  
  labs(y = "% mito", x = "") +
  ylim(0,5)+
  theme_minimal(base_size=16)+
  ggtitle("v2")


# %%
fig(3.5,5)
ggplot(feat.count.df, aes(x = "", y = ribo_fraction)) + 
  geom_violin(fill = "lightblue2", alpha = 0.5) +  
  geom_jitter(width = 0.1, size = 3, color='gray33' ,alpha = 1) +  
  labs(y = "% ribo", x = "") +
  ylim(0,10)+
  theme_minimal(base_size=16)+
  ggtitle("v2")


# %%
fig(6,5)
ggplot(feat.count.df, aes(x = lib_size)) + 
  geom_histogram(binwidth = 1e6, fill = "mediumseagreen", color = "white", alpha = 1) +
  labs(x = "lib_size", y = "Number of Samples") +
  theme_minimal(base_size=16)+
  ggtitle("v2")

# %%
fig(6,5)
ggplot(feat.count.df, aes(x = feature_counts)) + 
  geom_histogram(binwidth = 500, fill = "turquoise4", color = "white", alpha = 1) +
  labs(x = "Feature Counts", y = "Number of Samples") +
  theme_minimal(base_size=16)+
  ggtitle("v2")

# %%
fig(6,5)
ggplot(feat.count.df, aes(x = feature_counts_overfive)) + 
  geom_histogram(binwidth = 500, fill = "turquoise4", color = "white", alpha = 1) +
  labs(x = "Feature Counts >5", y = "Number of Samples") +
  theme_minimal(base_size=16)+
  ggtitle("v2")

# %%
fig(6,5)
ggplot(feat.count.df, aes(x=lib_size, y=feature_counts_overfive, color=mito_fraction)) + geom_point()+ 
  scale_y_continuous(breaks = seq(0, 40000, 1000)) +
  theme_minimal(base_size=16) +
  geom_hline(yintercept = 5000, color = "red", linetype = "dashed") + geom_vline(xintercept = 500000, color = "red", linetype = "dashed")+
  theme(axis.text.x = element_text(angle = 45))+
  ggtitle("v2")
