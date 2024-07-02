# inputs:
## 1. Taiji main output direactory
## 2. meta data

# outputs:
# data:
## 1. pr.csv: processed PageRank scores
## 2. rna.csv: gene expression of all genes
## 3. rna_TF.csv: gene expression of all TFs
## 4. mycolors.rds: fixed color palette, will be used in all downstream analysis
# plots:
## 1. age distribution 
## 2. sex distribution
## 3. cluster purity distribution

# import libraries
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(dplyr)

# set up ----------------------------------------------------
## create a sub-dir 'post-analysis/' under Taiji main output directory
maindir <- "/home/jupyter/output_20230310/post-analysis/" 
setwd(maindir)
set.seed(42)
tmp <- lapply(list.files(path = "scripts/functions/", pattern = "*\\.R",full.names = T),source)

## ------------- load rna data ----------------------------------
readRNA <- function(dir){
  
  # read raw file
  fl <- paste0(dir,"/expression_profile.tsv")
  rna <- read.table(fl, check.names = F)
  row.names(rna) <- toupper(rownames(rna))
  row.names(rna) <- gsub("\\(.*","",toupper(rownames(rna))) # remove "(UNKNOWN_TAIJI)"
  rna <- rna[complete.cases(rna),]
  
  # column-wise TPM-normalized
  ## remove columns whose expression values are all zeroes
  rna <- rna[,colSums(rna)>0]
  rna <- 10^6 * as.data.frame(scale(as.matrix(rna), center = F, scale = colSums(as.matrix(rna))))
  
  # remove rows whose expression levels are all zeroes
  rna <- rna[rowSums(rna)>0,]
  
  return(rna)
}

rna <- readRNA("../RNASeq/")
genes <- rownames(rna)


##load pagerank data -------------------------
load_pr <- function(dir,samples=NULL){
  fl <- paste0(dir,"/GeneRanks.tsv")
  df <- read.table(fl, check.names = F)
  df <- df[!rownames(df) %in% c("DNMT1"),]
  
  # some samples have zero expression value for all genes and need to be removed (optional)
  # some samples have small ATAC-seq peaks so pagerank is abnormal
  if (!is.null(samples)){
      df <- df[,colnames(df) %in% samples]
  }
  return(df)
}

pr <- load_pr("../")

# # ------------- read group info ----------------------------------
meta <- read.csv("metadata_20230310_expand.csv") %>% dplyr::filter(id %in% colnames(pr)) %>%
        dplyr::inner_join(meta_s, by = "subject.subjectGuid", suffix = c("","")) %>%
        dplyr::filter(!preClust == "Doublet") %>%
        dplyr::mutate(subject.group2 = ifelse(subject.group=="HC", subject.group, "non-HC")) %>%
        dplyr::mutate(subject.group3 = ifelse(subject.group%in%c("Converter","Non-converter"), "pre-RA",subject.group))


# unify meta and pr and rna
## unify colomns and TFs
All_TFs <- intersect(rownames(pr), rownames(rna))
samples <- unique(meta$id)

pr <- pr[All_TFs,samples]
rna_TF <- rna[All_TFs,samples] # subset of rna with only TFs' expression
rna <- rna[,samples]



# save to file
write.table(pr, "pr.tsv", sep = "\t", quote=F)
write.table(rna_TF, "rna_TF.tsv", sep = "\t", quote=F)
write.table(meta, "meta.tsv", sep='\t', quote=F, row.names=F)


## fix the annotation colors-------------------------------------------------------------
n1 = length(unique(meta$subject.subjectGuid)) # check how many patients are included, if larger than 12, then need to create a customized color palette
n2 = length(unique(meta$preClust))
n3 = length(unique(meta$subject.group))
n4 = length(unique(meta$subject.group2))
n5 = length(unique(meta$subject.group3))


qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(22)
mycolors <- list(subject.id = sample(col_vector,n1),
                 cell.type= colorRampPalette(brewer.pal(11,"Spectral"))(n2),
                 subject.group = brewer.pal(n3,"Set1")[c(2,1,3,4)],
                 subject.group2 = brewer.pal(8,"Set1")[c(3,5)],
                 subject.group3 = brewer.pal(8,"Set1")[c(1,3,6)],
                 TF.group = brewer.pal(5,"Set1")[c(2,1,3,4,5)],
                 kmeans = brewer.pal(5,"Set2")[1:5]
)
names(mycolors$subject.id) <- sort(unique(meta$subject.subjectGuid))
names(mycolors$cell.type) <- sort(unique(meta$preClust))
names(mycolors$subject.group) <- sort(unique(meta$subject.group))
names(mycolors$subject.group2) <- sort(unique(meta$subject.group2))
names(mycolors$subject.group3) <- sort(unique(meta$subject.group3))
names(mycolors$TF.group) <- c(sort(unique(meta$subject.group)),"multitasker")
names(mycolors$kmeans) <- paste0("C",1:5)

saveRDS(mycolors, "mycolors.rds")

# age distribution
df <- meta %>% dplyr::distinct(subject.subjectGuid, .keep_all = TRUE)
p <- ggplot(data=df, aes(x=subject.age, fill=subject.group))+
        geom_histogram(position = "identity") + 
        theme_bw()+ 
        scale_fill_manual(values = mycolors$subject.group)+
        theme(text = element_text(size = 20))
pdf("age_dist.pdf", width = 6, height = 4)
print(p)
dev.off()

# sex distribution
df <- meta %>% dplyr::distinct(subject.subjectGuid, .keep_all = TRUE) %>% group_by(subject.biologicalSex, subject.group) %>%
        mutate(count = n()) %>% distinct(subject.group, subject.biologicalSex, count)
p <- ggplot(data=df, aes(x=subject.biologicalSex, y=count, fill=subject.group))+
        geom_bar(position = "dodge", stat = "identity") + 
        theme_bw()+ 
        scale_fill_manual(values = mycolors$subject.group)+
        theme(text = element_text(size = 20))
pdf("sex_dist.pdf", width = 6, height = 4)
print(p)
dev.off()

# cluster purity distribution
p1 <- ggplot(meta, aes(x=preClust, y=purity, fill=subject.group))+
    geom_boxplot()+
    scale_fill_manual(values = mycolors$subject.group)+
    theme(legend.position = "none", axis.text.x = element_text(angle = 45))+
    ggtitle("cluster purity distribution")+
    facet_wrap(~subject.group, ncol = 1)
pdf("purity_dist.pdf",width = 10, height = 5)
print(p1)
dev.off()
