# set working directory
maindir <- "../code/Fig_signature/" # needs to be modified
datadir <- "../../data/Fig_signature/"
resultdir <- "../../results/"
setwd(maindir)

# import functions
source('findOptimal.R') # find optimal hyper-parameters for Kmeans clustering
source('specificTFs_cell.R') # identify Kmeans group-specific TFs. Dependencies are unpairedtTest.R, wilcoxonTest.R
source('utils.R') # some basic functions
source('gsea.R') # GO analysis

# import libraries
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(data.table)
library(compositions)
library(UpSetR)

# set plot parameters, which is suitable for CNS publication
plot.format = theme(
    plot.background=element_blank(),
    panel.grid=element_blank(),
    panel.background=element_blank(),
    panel.border=element_rect(color="black", linewidth=0.5, fill=NA),
    axis.line=element_blank(),
    axis.ticks=element_line(color="black",linewidth=0.5),
    axis.text=element_text(color="black", size=7),
    axis.title=element_text(color="black", size=7),
    plot.title=element_text(color="black", size=7),
    legend.background=element_blank(),
    legend.key=element_blank(),
    legend.text=element_text(color="black", size=7),
    legend.title=element_text(color="black", size=7))

# set cell types of interest
selected_cell_states <- c("B intermediate","B memory", "B naive", "CD14 Mono", "CD16 Mono","CD4 Naive","CD8 Naive","CD4 TCM","CD8 TEM","MAIT","NK","NK_CD56bright","Treg")

# prepare input
pr <- read.csv(paste0(datadir, 'pr.csv'), row.names = 1)
meta <- read.csv(paste0(datadir, 'meta.csv'), row.names = 1)

# load customized color palette
mycolors <- readRDS(paste0(datadir, 'mycolors.rds'))
print(mycolors)

# calculate the total number of scRNA-seq cells and scATAC-seq cells
sum(meta$RNA_cells)
sum(meta$ATAC_cells)
total_cells <- sum(meta$RNA_cells) + sum(meta$ATAC_cells)
print(paste0("total cells: ", total_cells))

# median of cells per cluster
tmp <- meta |> mutate(total_cells = RNA_cells + ATAC_cells)
median(tmp$total_cells)
print(paste0("median of cells per cluster: ", median(tmp$total_cells)))

# median of cells per sample
tmp <- meta |> mutate(total_cells = RNA_cells + ATAC_cells) |> group_by(subject.subjectGuid) |> summarize(n = sum(total_cells))
median(tmp$n)
print(paste0("median of cells per sample: ", median(tmp$n)))

# percentage of each cell type
tmp <- meta |> mutate(total_cells = RNA_cells + ATAC_cells) |> group_by(preClust) |> summarize(n = sum(total_cells)) |> arrange(-n) |> mutate(pct=round(100*n/total_cells, 2))
sum(head(tmp, n= 13)$n)/total_cells # top 13 cell types total percentage
print(paste0("top 13 cell types total percentage: ", sum(head(tmp, n= 13)$n)/total_cells))


# across groups
tmp <- meta |> mutate(total_cells = RNA_cells + ATAC_cells) |> group_by(preClust, sample, subject.group3) |> summarize(n = sum(total_cells)) |> arrange(-n) 
# |> mutate(pct=round(100*n/total_cells, 2))
df4 <- tmp %>% filter(preClust %in% selected_cell_states) %>% dplyr::group_by(sample) %>% mutate(total=sum(n)) %>% dplyr::group_by(preClust) %>% mutate(per=round(100*n/total,2))

print('Plotting cell type distribution across groups...')


p <- ggplot(df4, aes(x = preClust, y = per, fill = subject.group3)) +
    geom_boxplot(outlier.size = 0.5) +
    # geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(1), dotsize = 0.1)+
    labs(x = "", y = "Percent of labelled cells") +
    scale_fill_manual(name="group",values = mycolors$subject.group3) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)
    )

pdf(paste0(resultdir, "13_cell_states_dist_3_groups.pdf"), width=4.5, height=3.5)
print(p)
dev.off()

#### We applied CLR (Centered Log-Ratio) transformation to handle compositional contraints of compositional data, such as percentages that sum to 100%.Then we performed standard statistical tests like Kruskal-Wallis test.
print('Testing if cell types are significantly different across groups...')
comp_data <- meta |> mutate(total_cells = RNA_cells + ATAC_cells) |> filter(preClust %in% selected_cell_states)|> group_by(preClust, sample) |> summarize(n = sum(total_cells)) |> arrange(-n) |> 
    tidyr::pivot_wider(names_from = 'preClust', values_from = 'n', values_fill = 0) |> tibble::column_to_rownames('sample')
comp_data <- comp_data/rowSums(comp_data)

# CLR transformation
clr_data <- clr(comp_data)

# get disease group
disease_groups <- meta |> distinct(sample, subject.group3) |> arrange(factor(sample, levels=rownames(clr_data))) |> pull(subject.group3)
disease_groups <- factor(disease_groups, labels = c('ERA','At-Risk', 'CON'))

# perform Kruskal-Wallis test for each cell type
kruskal_test <- function(data) {
  results <- sapply(1:ncol(data), function(i) {
    test <- kruskal.test(data[,i] ~ disease_groups)
    c(statistic = test$statistic, p.value = test$p.value)
  })
  rownames(results) <- c("statistic", "p.value")
  colnames(results) <- colnames(data)
  return(results)
}

clr_results <- kruskal_test(clr_data)
print(clr_results[,selected_cell_states])

# statistics of purity scores
print(paste0("mean of purity scores: ", mean(meta$purity)))
print(paste0("sd of purity scores: ", sd(meta$purity)))

p1 <- ggplot(meta|>filter(preClust%in%selected_cell_states), aes(x=preClust, y=purity, fill=preClust))+
    geom_boxplot()+
    scale_fill_manual(values = mycolors$cell.type)+
    ggtitle("cluster purity distribution")+
    theme(
      axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)
    )+plot.format

# cluster purity distribution across kmeans
pdf(paste0(resultdir, "purity_across_13_cell_types.pdf"),width = 5, height = 3.5)
print(p1)
dev.off()

# Kmeans clustering
## The following inputs are needed for Kmeans clustering:
### 1. PageRank scores matrix, which is the output of Taiji. Each row is one TF, each column is one sample, and each cell is PageRank score. 
### 2. meta data of each sample. Typical meta data includes cell types, participant id, clinical status, etc.

# feature filtering
##The first step is to prepare data for Kmeans clustering. Before PCA reduction, we need to filter the features since some features are outliers, i.e. some TFs have abnormally high scores in very few samples. We'll filter out these outliers using function `find_outlier_feature`.
### 1. `iqr_multiplier` determines the upper bound for outliers (default is 50)
### 2. `sample_threshold`: minimum number of samples that must exceed the upper bound for a feature to be considered an outlier (default is 2)

print('Removing outliers...')
outlier_f <- find_outlier_feature(pr, iqr_multiplier = 100, sample_threshold = 3)
print(paste0("number of outliers: ", length(outlier_f)))

outlier_tfs <- rownames(pr)[outlier_f]
writeLines(outlier_tfs, paste0(resultdir, 'outlier_TFs_iqr100_s3.txt'))

tag <- "" # can be customized for different versions of result
pr <- pr[setdiff(rownames(pr),outlier_tfs),]
pr_normed2 <- zscore(pr) ## zscore of pagerank row-wise

# load data first
data <- pr_normed2

# perform pca
pca <- data %>% as.matrix %>% t %>% prcomp 
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

# visualize first PCs
mydata <- pca %>% .$x %>% data.frame %>% mutate(Samples = rownames(.)) %>% 
    mutate(cell.type = meta[match(rownames(.), meta$id),"preClust"]) %>% 
    mutate(subject.group = meta[match(rownames(.), meta$id),"subject.group3"])

myPlot <- ggplot(mydata, aes(x = PC1, y = PC2, colour = cell.type, shape=subject.group, stroke=1.5)) + geom_point(size = 5) + scale_color_manual(values = mycolors$cell.type) + 
        plot.format + scale_shape(solid = FALSE) + ggtitle("PCA by PageRank") +
        xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
        ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))

pdf(paste0(resultdir, "PCA_",tag,"_PC1_2_group3.pdf"), width = 7, height = 5)
print(myPlot)
dev.off()

# plot pca cumulating variance
df2 <- data.frame(x=c(1:length(pca$sdev)),y=cumsum(percentVar))
p2 <- ggplot(df2)+aes(x=x,y=y)+geom_point()+geom_line()+
    labs(x="Principal Component",y="Cumulative Proportion of \nVariance Explained")+
    plot.format
print(p2)

## save to file
pdf(paste0(resultdir, "PCA_",tag,"_cumul_variance.pdf"), width=2, height=2)
print(p2)
dev.off()


# hyperparameter selection
## After manually checking the above cumulative PCA plot, we decide to use PC=500 to reduce the data based on "elbow" method, which explained 85% variance. 
## Then we need to determine the best distance metric and number of K for Kmeans clustering. To evaluate the clustering quality, we computed:
### 1. **average silhouette width** which combines both *cohesion* (how close data points in a cluster are to each other)
### 2. **gap statistic** to determine the number of clusters
### 3. **elbow method** to determine the number of principal components

## We used five common distance metrics:
### - Euclidean distance
### - Manhattan distance
### - Kendall correlation
### - Pearson correlation
### - Spearman correlation

## Number of k is selected from the range [3, 20]. The lower limit is considering the number of disease states while the upper limit is the balance of searching space and searching speed.

PCNo <- 500
df <- data %>% as.matrix %>% t
data_reduced <- as.data.frame(prcomp(df, rank. = PCNo)$x)
write.csv(data_reduced, paste0(resultdir, 'data_reduced_PC',PCNo,'_',tag,'.csv'))

# this step takes ~3hr
set.seed(3)
df2 <- findOptimal(data_reduced, max_k = 10, cores = detectCores())
saveRDS(df2,paste0(resultdir, 'hyper_opt_PC',PCNo,'_',tag,'.rds'))

# save to file
pdf(paste0(resultdir, "kmeans_param_PC",PCNo,"_",tag,".pdf"),width=3.5,height = 2)
df2$silhouette_plot+plot.format
df2$elbow_plot+plot.format
df2$gap_plot+plot.format
dev.off()

# Kmeans clustering
## From the distance metrics plot, we can select the best hyperparam combos: 
### - K=5
### - distance metric=Pearson Correlation.

# get final Kmeans clustering result
set.seed(3)
clusterNo <- 5
PCNo <- 500
df <- data %>% as.matrix %>% t
data_reduced <- as.data.frame(prcomp(df, rank. = PCNo)$x)
Cluster <- kmeans(data_reduced, centers = clusterNo, nstart = 25,iter.max = 50) 

df2 <- as.data.frame(Cluster$cluster) %>% tibble::rownames_to_column("id") %>% dplyr::rename(kmeans=`Cluster$cluster`) 
df2$kmeans <- paste0("G",df2$kmeans)

# get updated meta file
info <- meta %>% inner_join(df2, by = "id")
print(head(info))
print(dim(info))

# adjust the order of Kmeans: G3-->G4; G4-->G5; G5-->G3 (K=5)
info <- info |> mutate(kmeans = ifelse(kmeans=='G5', 'G2', ifelse(kmeans=='G2', 'G4', ifelse(kmeans=='G4', 'G1', ifelse(kmeans=='G1', 'G5', kmeans)))))

# save to file
write.csv(info, paste0(resultdir, 'meta_kmeans_',tag,'_',clusterNo,'.csv'))

# cluster distribution across Kmeans groups
## cell type distribution across Kmeans groups

file <- paste0(resultdir, 'meta_kmeans_',tag,'_',clusterNo,'.csv')
test_group <- 'preClust'
f1 <- paste0("id ~ kmeans + ", test_group)
f2 <- paste0("id ~ ", test_group)


df1 <- read.csv(file) %>% 
    stats::aggregate(as.formula(f1), data = ., FUN = function(x){NROW(x)})

# add null distribtution
df2 <- read.csv(file) %>%
    stats::aggregate(as.formula(f2), data = ., FUN = function(x){NROW(x)}) %>% mutate(kmeans="All")
df3 <- rbind(df1,df2) %>% arrange(kmeans)

# write to file
tmp <- df3 %>% tidyr::pivot_wider(names_from = preClust, values_from = id, values_fill = 0) 
write.csv(tmp,paste0(resultdir, 'cell_states_freq_across_kmeans_',tag,'_',clusterNo,'.csv'),quote=F, row.names=F)

# normalize
df3 <- read.csv(paste0(resultdir, 'cell_states_freq_across_kmeans_',tag,'_',clusterNo,'.csv')) %>% tibble::column_to_rownames('kmeans') 
df4 <- df3/rowSums(df3)

# visualize all cell types
p <- pheatmap(df4, cluster_rows=FALSE, cluster_cols=FALSE, gaps_row = 1, cellheight = 12, cellwidth = 12, angle_col = 90)

# select major cell states and visualize
print(selected_cell_states)
df5 <- df4[,sub(' ','.',selected_cell_states)]

p <- pheatmap(df5, cluster_rows=FALSE, cluster_cols=FALSE, gaps_row = 1, cellheight = 12, cellwidth = 12, angle_col=90)

## save to file
pdf(paste0(resultdir, 'cell_state_freq_across_kmeans_',tag,'_',clusterNo,'.pdf'), width = 4, height = 4)
print(p)
dev.off()

# At-Risk/ERA vs CON ratio across Kmeans groups
## we want to test if the ratio of At-Risk/ERA vs CON is significantly different across Kmeans groups
file <- paste0(resultdir, 'meta_kmeans_',tag,'_',clusterNo,'.csv')
test_group <- 'subject.group2'
# test_group <- 'subject.group'
f1 <- paste0("id ~ kmeans + ", test_group)
f2 <- paste0("id ~ ", test_group)


df1 <- read.csv(file) %>% 
    stats::aggregate(as.formula(f1), data = ., FUN = function(x){NROW(x)})

# add null distribtution
df2 <- read.csv(file) %>%
    stats::aggregate(as.formula(f2), data = ., FUN = function(x){NROW(x)}) %>% mutate(kmeans="All")
df3 <- rbind(df1,df2) %>% arrange(kmeans)


# write to file
tmp <- df3 %>% tidyr::pivot_wider(names_from = !!as.name(test_group), values_from = id, values_fill = 0)
write.csv(tmp,paste0(resultdir, 'disease_state_freq_across_kmeans_',test_group,'_',tag,'_',clusterNo,'.csv'),quote=F, row.names=F)

# df3[df3$kmeans==x & df3$subject.group!="Non-converter", "id"]
# ratios <- df3 |> tidyr::pivot_wider(names_from = test_group, values_from = 'id', values_fill = 0) |> as.data.frame() |> mutate(ratio=(Converter+ERA)/`Non-converter`)
ratios <- df3 |> tidyr::pivot_wider(names_from = test_group, values_from = 'id', values_fill = 0) |> as.data.frame() |> mutate(ratio=(`At-Risk/ERA`/CON))
null = ratios[ratios$kmeans=='All','ratio']
ratios$adjusted.ratio <- ratios$ratio/null
print(ratios)

p4 <- ggplot(ratios, aes(x = kmeans, y = adjusted.ratio, fill=kmeans)) +
    geom_bar(stat = "identity", colour = "black") +
    labs(x = "", y = "ratio") + 
    scale_fill_manual(values = c("#808080",brewer.pal(clusterNo,"Set2")[1:clusterNo])) +
    theme_bw() +
    theme(
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
    )


## save to file
pdf(paste0(resultdir, 'disease_state_ratio_across_kmeans_',test_group,'_',tag,'_',clusterNo,'.pdf'), width = 4, height = 4)
print(p4)
dev.off()

## statistical test
## From above plot, it's obvious that G2 is enriched in At-Risk/ERA compared to controls. 
## Then we used chi-squared goodness of fit test to see if the enrichment is significant or not compared to the null distribution. 
## R function `chisq.test` can be used as follow:
### `chisq.test(x,p)`, where x is a numeric vector representing the observed probability and p represents the expected probability.

# test for G2, subject.group2
test_k = 'G2'
test_group = 'subject.group2'
a <- df3 |> filter(kmeans=='All' & !!as.name(test_group)!='CON') |> pull(id) 
b <- df3 |> filter(kmeans=='All' & !!as.name(test_group)=='CON') |> pull(id) 
a2 <- df3 |> filter(kmeans==test_k & !!as.name(test_group)!='CON') |> pull(id) 
b2 <- df3 |> filter(kmeans==test_k & !!as.name(test_group)=='CON') |> pull(id) 

x <- c(a2,b2)
p <- c(a/(a+b), b/(a+b))
print(chisq.test(x,p=p))


## G2 clusters ratio per cell type in CON and At-Risk/ERA
## From the previous analysis, we already know:
### - G2 is significantly enriched in At-Risk/ERA
### - G2 is multi-lineage group
## We then want to see the distribution of the enrichment across cell types. In other words, what cell types in G2 that contributed more in the enrichment?

C_target = 'G2'
g <- "subject.group2"
meta3 <- read.csv(file, row.names = 1) %>% 
        mutate(preClust2 = ifelse(preClust %in% c("B memory", "B naive", "B intermediate"), "B cell", 
                                  ifelse(preClust %in% c("CD16 Mono", "CD14 Mono"), "Monocytes", 
                                         ifelse(preClust %in% c("NK", "NK_CD56bright"), "NK cell", preClust))),
                kmeans2 = ifelse(kmeans==C_target, kmeans, "other"))


# get the number of total clusters across cell types
total <- meta3 |> group_by(preClust2, !!as.name(g)) |> summarise(total=n())


# get the number of G2 clusters across cell types
sub <- meta3 |> filter(kmeans==C_target) |> group_by(preClust2, !!as.name(g)) |> summarise(sub=n())


# calculate the percentage of G2 clusters per cell type in CON and At-Risk/ERA respectively
df <- sub %>% right_join(total, by=c("preClust2",g)) %>% mutate_all(~replace(., is.na(.), 0)) %>% mutate(pct = 100*sub/total) %>%
            mutate(label=paste0(sub,"(", sprintf("%.0f", pct), "%)"))


# save to file
write.csv(df, paste0(resultdir, 'pct_G2_per_cell_type_two_groups_',tag,"_",clusterNo,'.csv'))

## statistical test
## We used chi-squared goodness of fit test to see if the percentage is significant different from CON and At-Risk/ERA. 
selected_cell_states2 <- c('B cell', 'Monocytes', 'CD4 Naive', 'CD8 Naive', 'CD4 TCM', 'CD8 TEM', 'MAIT', 'NK cell', 'Treg')


#### calculate p-values for all cell types
pvalues <- lapply(selected_cell_states2, function(x){
    print(x)
    tmp <- df %>% filter(preClust2==x) 
    a <- tmp |> filter(subject.group2=='At-Risk/ERA') |> pull(total)
    b <- tmp |> filter(subject.group2=='CON') |> pull(total)
    c <- tmp |> filter(subject.group2=='At-Risk/ERA') |> pull(sub)
    d <- tmp |> filter(subject.group2=='CON') |> pull(sub)
    
    null.prop <- c(a/(a + b), b/(a + b))
    res <- chisq.test(c(c, d), p = null.prop)
    print(res$p.value)    
    return(res$p.value)
})

names(pvalues) <- selected_cell_states2

## select cell lineages with significant enrichment
cells <- names(pvalues)[which(pvalues <= 0.3)]

## plot vertical histogram
df2 <- df %>% filter(preClust2 %in% cells)
df2[[g]] <- factor(df2[[g]], levels = sort(unique(df2[[g]]), decreasing = T))

p1 <- ggplot(df2, aes(x = preClust2, y = pct, fill = !!as.name(g))) + 
    geom_bar(position = position_dodge(), stat = "identity", width = .7) +
    scale_fill_manual(values = mycolors[[g]]) + 
    geom_text(aes(label = label), position = position_dodge(width = .7), vjust=0, size = 3) +
    plot.format + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5), axis.title.x = element_blank())


# no label version
p2 <- ggplot(df2, aes(x = preClust2, y = pct, fill = !!as.name(g))) + 
    geom_bar(position = position_dodge(), stat = "identity", width = .7) +
    scale_fill_manual(values = mycolors[[g]]) + 
    plot.format + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5), axis.title.x = element_blank())


df3 <- df |> filter(preClust2 %in% selected_cell_states2)


p3 <- ggplot(df3, aes(y = !!as.name(g), x = pct, fill = !!as.name(g))) + 
    geom_bar(position = position_dodge(), stat = "identity", width = .7) +
    scale_fill_manual(values = mycolors[[g]]) + 
    facet_wrap(~preClust2)+
    geom_text(aes(label = label), position = position_dodge(width = .7), vjust=0, size = 3) +
    plot.format + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5), axis.title.x = element_blank(),
                        legend.key = element_blank(), strip.background = element_rect(colour="black", fill=NA)) 



p4 <- ggplot(df3, aes(y = !!as.name(g), x = pct, fill = !!as.name(g))) + 
    geom_bar(position = position_dodge(), stat = "identity", width = .7) +
    scale_fill_manual(values = mycolors[[g]]) + 
    facet_wrap(~preClust2)+
    plot.format + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5), axis.title.x = element_blank(),
                        legend.key = element_blank(), strip.background = element_rect(colour="black", fill=NA)) 



## save to file
pdf(paste0(resultdir, 'pct_G2_per_cell_type_two_groups_',g,'_',tag,"_",clusterNo,'.pdf'), width = 5, height = 4)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()


# signature TFs, genes, and pathways
## Kmeans group-specific TFs
## The heatmap shows distinct TF profiles across Kmeans groups. We used Wilcoxon rank-sum test (a.k.a Mann-Whitney U test) to identify Kmeans group-specific TFs. 
## We selected Wilcoxon test because our data doesn't fit normal distribution. If data is normally or log-normally distributed, t-test can be applied here as well.
## To be more specific, we divided the clusters into two groups: target group and background group. 
### Target group included the clusters in the Kmeans group of interest, for example G2. 
### Background group included all the remaining clusters. 
## We experimented on several combinations of p-value and log2 fold change. We chose p-value of 0.01 and log2FC of 0.5 considering a reasonable number of TFs for all 5 groups.
tag <- "" 
clusterNo <- 5

# create dir for results
dir.create(file.path(resultdir, paste0('kmeans_',clusterNo)), showWarnings = FALSE)
dir.create(file.path(resultdir, paste0('kmeans_',clusterNo,'/GO')), showWarnings = FALSE)

setwd(paste0(resultdir,"kmeans_",clusterNo))

## identify kmeans cluster-specific TFs
main <- function(p = 0.001, d = 1, data=pr, meta=meta2, group_by='kmeans'){
    tryCatch({
        subject.types <- unique(meta[[group_by]])
        ls <- lapply(subject.types,function(y) 
          specificTFs_cell(y, p = p, d = d, df = data, group = meta, group_name = group_by, test = "wilcoxon"))

        ### GO analysis
        L <- list.files(path = "./", pattern = paste0(".*_p",p,"_d",d,".*.txt"), full.names = T)
        print(L)
        lapply(L, function(x) gsea(x, output_dir = "GO/", showCategory = 20, go.height = 10, go.width = 7))
    },error=function(e){})

}                     

# load data and meta data
meta2 <- read.csv(paste0(resultdir,'meta_kmeans_',tag,'_',clusterNo,'.csv'), row.names=1)


# this step takes several minutes
main(data=pr, meta=meta2, p = 0.01, d = 0.5)

L <- list.files(pattern = "G.*p0.01_d0.5.*txt", full.names = T)
TFs <- lapply(L, function(x) readLines(x))
names(TFs) <- paste0("G",1:clusterNo)

p <- upset(fromList(TFs), 
      nintersects = 40, 
      nsets = 6, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
      )

print(p)

pdf(paste0(resultdir,'kmeans_group_TFs_upset.pdf'), width = 4, height = 4)
print(p)
dev.off()

# total number of unique TFs
n_TFs <- unlist(TFs) |> unique() |> length()
print(paste0("Total number of unique TFs: ", n_TFs))

# read signature TFs
# pathways <- c('SUMOylation of intracellular receptors','Transcriptional regulation by RUNX2','YAP1- and WWTR1 (TAZ)-stimulated gene expression','Deactivation of the beta-catenin transactivating complex','NOTCH3 Intracellular Domain Regulates Transcription')
tfs <- read.csv('GO/G2_p0.01_d0.5_gene409_wilcoxon_GO.csv') |> filter(grepl('SUMO|RUNX2|YAP1|beta-catenin|NOTCH3', Description)) |> 
    select(Description, p.adjust,geneID) |> mutate(D2 = ifelse(grepl('SUMO', Description), 'SUMO', 
                           ifelse(grepl('RUNX2', Description), 'RUNX2', 
                                  ifelse(grepl('NOTCH3', Description), 'NOTCH3',                            
                                         ifelse(grepl('beta-catenin', Description), 'Wnt', 
                                                ifelse(grepl('YAP1', Description), 'YAP1', 'SUMO')))))) |>  group_by(D2) 


# signature TF list
l <- c('SUMO','RUNX2','YAP1','Wnt','NOTCH3')
regulators <- lapply(l, function(x) {
    r <- paste0(tfs |> filter(D2==x) |> pull(geneID), collapse='/')
    strsplit(r,'/')[[1]] |> unique()
})
names(regulators) <- l

saveRDS(regulators, paste0(resultdir,'signature_tfs.rds'))

p <- upset(fromList(regulators), 
      nintersects = 40, 
      nsets = 6, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1,set_size.show = TRUE
      )

print(p)

pdf(paste0(resultdir,'signature_TFs_upset.pdf'), width = 4, height = 3)
print(p)
dev.off()


L <- list.files(path = paste0(resultdir, '/kmeans_',clusterNo,'/GO/') , pattern = 'p0.01_d0.5.*csv', full.names = T)


## This plot shows the top pathways for each Kmeans group. We curated a list of G2-specific pathways associated with RA pathogenesis:
### 1. Transcriptional regulation by RUNX2
### 2. YAP1- and WWTR1 (TAZ)-stimulated gene expression
### 3. SUMOylation of intracellular receptors
### 4. NOTCH3 Intracellular Domain Regulates Transcription
### 5. Deactivation of the beta-catenin transactivating complex

# df <- do.call("rbind", lapply(L, function(x){read.csv(x) %>% filter(grepl('R-HSA', ID) & p.adjust<=0.05) %>% slice_min(order_by = p.adjust, n=30)%>%mutate(group=gsub('_.*','',gsub('.*/','',x)))}))
df <- do.call("rbind", lapply(L, function(x){read.csv(x) %>% filter(p.adjust<=0.05) |> mutate(ID_class = factor(gsub("[^a-zA-Z]", "", ID), levels = c('RHSA','GO','hsa'))) |>
                                             arrange(ID_class, p.adjust) |> head(n=35)|>
                                             mutate(group=gsub('_.*','',gsub('.*/','',x)))}))

df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text=x)))
df <- df %>% arrange(group, GeneRatio)
df$Description <- factor(df$Description, levels = unique(df$Description))
 

# selected representative pathways
# shortened version of plot
selected_pathways <- c('Signaling by NTRKs',
                       'SUMOylation of intracellular receptors','Transcriptional regulation by RUNX2','YAP1- and WWTR1 (TAZ)-stimulated gene expression','NOTCH3 Intracellular Domain Regulates Transcription','Deactivation of the beta-catenin transactivating complex',
                       'cartilage development', 'alpha-beta T cell activation',
                       'Formation of the beta-catenin:TCF transactivating complex','Repression of WNT target genes','Transcriptional regulation by RUNX3','Interleukin-4 and Interleukin-13 signaling',
                       'B cell activation')

df2 <- df %>% filter(Description %in% selected_pathways)
df2$Description <- factor(df2$Description, levels =selected_pathways)

p <- ggplot(df2, aes(x=group, y=Description, size=Count, color=p.adjust))+
    geom_point()+scale_color_gradient(low="red",high="blue", name="p.adjust")+theme_bw()+
    ylab("")+xlab("")+plot.format


# save to file
pdf(paste0(resultdir,'summary_representative_pathways.pdf'), width = 6, height = 3)
print(p)
dev.off()

## signature pathways across cell types
## Next, we identified G2-specific TFs along with pathways for each cell type. 
## Target group included the G2 clusters in the cell type of interest and the background included the remaining clusters of the same cell type. 
## p-value of 0.01 and log2 fold change of 0.5 were used for calling the specific TFs. 

# create dir for results
dir.create(file.path(resultdir, paste0('kmeans_cell_type_',clusterNo)), showWarnings = FALSE)
dir.create(file.path(resultdir, paste0('kmeans_cell_type_',clusterNo,'/GO')), showWarnings = FALSE)

setwd(paste0(resultdir,'kmeans_cell_type_',clusterNo))

meta3 <- meta2 |> mutate(preClust2 = ifelse(preClust %in% c("B memory", "B naive", "B intermediate"), "B cell", ifelse(preClust %in% c("CD16 Mono", "CD14 Mono"), "Monocytes", ifelse(preClust %in% c("NK", "NK_CD56bright"), "NK cell", preClust)))) |> mutate(cell_k=paste0(preClust2,'_',kmeans))


selected_cell_states3 <- c('B cell', 'Monocytes', 'CD4 Naive', 'CD4 TCM', 'CD8 Naive', 'CD8 TEM', 'NK cell', 'Treg')

# identify specific TFs for each cell type, which takes ~5 min
tmp <- lapply(selected_cell_states3,function(x){
    meta4 <- meta3 |> filter(preClust2==x)
    df <- pr[,meta4$id]
    ls <- specificTFs_cell(paste0(x,'_G2'), p = 0.01, d = 0.5, df = df, group = meta4, group_name = 'cell_k', test = "wilcoxon")
})

# GO analysis for each group of TFs, which takes ~10 min
L <- list.files(path = "./", pattern = paste0(".*_p",0.01,"_d",0.5,".*.txt"), full.names = T)
print(L)
lapply(L, function(x) gsea(x, output_dir = "GO/", showCategory = 20, go.height = 10, go.width = 7))

setwd("GO/")

L <- list.files(pattern = 'gene.*csv')
df <- do.call("rbind", lapply(L, function(x){read.csv(x) %>% filter(grepl('R-HSA', ID) & p.adjust<=0.05) %>% slice_min(order_by = p.adjust, n=50)%>%mutate(group=gsub('_.*','',x))}))


# write to file
write.table(df, 'Reactome_pathway_p0.05_G2_specific_TFs_across_cell_states_top50.tsv', quote=F, row.names=F, sep = '\t')

df <- read.table('Reactome_pathway_p0.05_G2_specific_TFs_across_cell_states_top50.tsv', sep = '\t', header = T)
df$GeneRatio <- sapply(df$GeneRatio, function(x) eval(parse(text=x)))
df <- df %>% arrange(group, GeneRatio)
df$Description <- factor(df$Description, levels = unique(df$Description))
df <- df %>% mutate(log10_p.adjust = ifelse(p.adjust<=1e-5, 5, -log10(p.adjust))) |> group_by(group) |> slice_min(order_by=p.adjust, n = 50)
tail(df)
dim(df)

p <- ggplot(df, aes(x=group, y=Description, size=Count, color=log10_p.adjust))+
    geom_point()+scale_color_gradient(low="blue",high="red", name="-log10(p.adjust)")+theme_bw()+
    ylab("")+xlab("")+plot.format+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))



# selected pathways
selected_pathways <- c('Activation of HOX genes during differentiation',
                       'NGF-stimulated transcription',
                       'Signaling by NOTCH1 in Cancer',
                       'Transcriptional regulation by RUNX3',
                       'Deactivation of the beta-catenin transactivating complex',
                       'NOTCH3 Intracellular Domain Regulates Transcription',
                       'SUMOylation of intracellular receptors', 
                       'YAP1- and WWTR1 (TAZ)-stimulated gene expression', 
                       'Transcriptional regulation by RUNX2')
selected_cell_states <- c('B cell','CD4 Naive','CD4 TCM','CD8 Naive','CD8 TEM','Monocytes','NK cell','Treg')
df3 <- df %>% filter(Description %in% selected_pathways & group %in% selected_cell_states)
df3$Description <- factor(df3$Description, levels=selected_pathways)

p2 <- ggplot(df3, aes(x=group, y=Description, size=Count, color=p.adjust))+
    geom_point()+
    scale_color_gradient(low="red",high="blue", name="p.adjust", limits=c(0, 0.05))+
    # scale_color_gradient2(
    #     low = "red",
    #     mid = "gray",
    #     high = "blue",
    #     midpoint = 0.05,
    #     limits = c(0, 1),
    #     name = "p.adjust"
    # ) +
    ylab("")+xlab("")+
    # guides(color='none')+
    plot.format+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save to file
setwd('../')
pdf(paste0(resultdir,'summary_representative_G2_pathways_across_cell_types.pdf'), width = 6, height = 2.5)
print(p2)
dev.off()

## Top regulated genes of signature TFs
## Now we have identified signature pathways along with involved TFs associated with RA pathogenesis. 
## Taiji generated the regulatory network showing the regulatory relationship between TF and its target genes (**regulatee**) with edge weight, which represents the regulatory strength. 
## Based on this, we can curate a list of representative regulatee genes by selecting the top 10 genes regulated by the signature TFs involved in each pathway ranked by the mean edge weight.
## Required input: network files "/taiji_output/Network/sample*/edges_combined.csv". Each sample/pseudo-bulk cluster has one network file.

setwd(resultdir)
# read signature TFs
pathways <- c('SUMOylation of intracellular receptors','Transcriptional regulation by RUNX2','YAP1- and WWTR1 (TAZ)-stimulated gene expression','Deactivation of the beta-catenin transactivating complex','NOTCH3 Intracellular Domain Regulates Transcription')
tfs <- read.csv(paste0(resultdir,'/kmeans_5/GO/G2_p0.01_d0.5_gene409_wilcoxon_GO.csv')) |> filter(Description %in% pathways) |> select(Description, p.adjust,geneID)


# signature TF list
regulators <- unique(unlist(lapply(pathways, function(x) {
    r <- tfs |> filter(Description==x) |> pull(geneID) 
    strsplit(r,'/')[[1]]
})))
print(paste0("number of signature TFs: ", length(regulators)))

# selected sample list
fg <- meta2 %>% filter(kmeans == 'G2' & subject.group2 == "At-Risk/ERA") %>% pull(id)
print(paste0("number of selected samples: ", length(fg)))

# # this step takes ~30 min
# # to save time, the output file is provided for demo purpose so that this step can be skipped
# get_csv_file(regulators, fg, '37_signature_TFs_vs_all_genes')

dt <- read.csv(paste0(datadir,'37_signature_TFs_vs_all_genes_regulatees.tsv'), sep=' ')
names(dt) <- c("regulator","regulatees", "weight","sample")
dt <- dt[!grepl("[0-9]{3}", dt$regulatees),]
dt <- dt |> group_by(regulator, regulatees) |> summarise(mean_edge_weight = mean(weight))
head(dt)
dim(dt)

tfs$gene <- unlist(lapply(pathways, function(x){
    regulators <- tfs |> filter(Description==x) |> pull(geneID) 
    regulators <- strsplit(regulators,'/')[[1]]
    regulatee <- dt |> filter(regulator %in% regulators) |> group_by(regulatees) |> summarise(mean_mean_ew=mean(mean_edge_weight)) |> slice_max(mean_mean_ew, n =10) |> pull(regulatees)    
    paste0(regulatee, collapse = '/')
}))


write.csv(tfs, paste0(resultdir,'Table_sig_pathways_TFs_regulatees.csv'))


# Heatmap
# change working directory to resultdir
setwd(resultdir)
print("Plotting heatmaps...")

# take the top ten TFs for each group and plot partial heatmap
L <- list.files(path = paste0('kmeans_',clusterNo,'/'), pattern = "G.*p0.01_d0.5.*txt", full.names = T)
TFs <- lapply(L, function(x) readLines(x)[1:10])
names(TFs) <- paste0("G",1:clusterNo)
annotation_row <- tibble("TF.group" = unlist(lapply(1:length(TFs), function(x) {rep(names(TFs[x]), length(TFs[[x]]))})),"TFs" = unlist(TFs))
annotation_row[annotation_row$TF.group=="G2","TFs"] <- c("ZNF304","HKR1","NFE4","CTCFL","ETV4","ZNF254","SOX9","FOXL2","VSX2","HEY1")
annotation_row <- annotation_row |> group_by(TFs) |> filter(n()==1)|> tibble::column_to_rownames("TFs")


get_heatmap <- function(df, info, row_ann, mycolors,...){
        
    info <- info %>% dplyr::filter(id %in% colnames(df)) %>% dplyr::arrange(kmeans, preClust)
   
    # sort df make sure the order of info and df matched
    df <- df[,match(info$id,names(df))]

    annotation_col = data.frame(
      # subject.group = info$subject.group,
      # cell.type = info$preClust,
      subject.group2 = info$subject.group2,
      subject.group3 = info$subject.group3,
      kmeans = info$kmeans
    )
    rownames(annotation_col) = names(df)
  
    df2 <- df[rownames(row_ann),rownames(annotation_col)]
    print(dim(df2))

    p1 <- pheatmap(df2, 
                 fontsize = 7, angle_col = 90, 
                 cellwidth = 0.2, cellheight = 7,
                 cluster_cols = F, cluster_rows = F, 
                 clustering_distance_cols = 'correlation', 
                 clustering_distance_rows = 'correlation', 
                 clustering_method = 'average',
                 annotation_row = row_ann,
                 annotation_col = annotation_col,
                 breaks = seq(-2, 4, by = 0.1),
                 show_colnames = F,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(-2, 4, by = 0.1))),
                 annotation_names_row = F, annotation_names_col = F,
                 annotation_colors = mycolors, border_color = NA,...)
    return(p1)
}

mycolors$TF.group <- mycolors$kmeans

p1 <- get_heatmap(pr_normed2, meta2, mycolors, row_ann=annotation_row)

# save to file
pdf(paste0("hp_",tag,"_",clusterNo,"_top10_TFs.pdf"))
print(p1)
dev.off()

# Heatmap of all TFs
get_heatmap <- function(df, info, mycolors,...){
        
    info <- info |> dplyr::filter(id %in% colnames(df)) |> arrange(kmeans,preClust)
   
    # sort df make sure the order of info and df matched
    df <- df[,match(info$id,names(df))]

    annotation_col = data.frame(
      subject.group = info$subject.group,
      cell.type = info$preClust,
      subject.subjectGuid = info$subject.subjectGuid,
      kmeans = info$kmeans
    )
    rownames(annotation_col) = names(df)
  
    p1 <- pheatmap(df, 
                 fontsize = 7, angle_col = 90, 
                 cellwidth = 0.3, cellheight = 0.3,
                 clustering_distance_cols = 'correlation', 
                 clustering_distance_rows = 'correlation', 
                 clustering_method = 'average',
                 annotation_col = annotation_col,
                 breaks = seq(-2, 4, by = 0.1),
                 show_colnames = F, 
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(-2, 4, by = 0.1))),
                 annotation_names_row = T, annotation_names_col = T,
                 annotation_colors = mycolors, 
                 border_color = NA,...)
    return(p1)
}

tag <- "" # can be customized for different versions of result
clusterNo <- 5
meta2 <- read.csv(paste0(resultdir,'meta_kmeans_',tag,'_',clusterNo,'.csv'), row.names=1)

p1 <- get_heatmap(pr_normed2, meta2, mycolors, cluster_cols=F)

pdf(paste0("hp_normed2_all_TFs_",tag,"_",clusterNo,".pdf"), height=6, width = 11)
print(p1)
dev.off()

# the total heatmap with partial labels
L <- list.files(path = 'kmeans_5/', pattern = "G.*p0.01_d0.5.*txt", full.names = T)
TFs <- unlist(lapply(L, function(x) readLines(x)[1:10]))

## heatmap sorted by cell type
hp_basic <- function(df,info,mycolors,...){
  
    info <- info %>% dplyr::filter(id %in% colnames(df)) %>% 
            dplyr::filter(preClust %in% selected_cell_states) %>%
            dplyr::arrange(preClust, kmeans, subject.group3)

    # sort df make sure the order of info and df matched
    df <- df[,match(info$id,names(df))]

    annotation_col = data.frame(
      cell.type = info$preClust,
      subject.group2 = info$subject.group2,
      subject.group3 = info$subject.group3,
      kmeans = info$kmeans
    )
    rownames(annotation_col) = names(df)

    p1 <- pheatmap::pheatmap(df, 
                    fontsize = 7, show_rownames = F,
                    angle_col = 90, show_colnames = F,
                    cluster_cols = F, cluster_rows = T, 
                    cellwidth = 0.3, cellheight = 0.1,
                    clustering_distance_cols = 'correlation', 
                    clustering_distance_rows = 'correlation', 
                    clustering_method = 'average',
                    annotation_col = annotation_col,
                    annotation_colors = mycolors, border_color = NA,
                    # breaks = seq(-2, 4, by = 0.1),
                    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(-2, 4, by = 0.1))),...)
    return(p1)
  
}

pr_normed2[pr_normed2 > 4] <- 4
pr_normed2[pr_normed2 < -2] <- -2

selected_cell_states <- c("B intermediate","B memory", "B naive", "CD14 Mono", "CD16 Mono","CD4 Naive","CD8 Naive","CD4 TCM","CD8 TEM","MAIT","NK","NK_CD56bright","Treg")
p2 <- hp_basic(pr_normed2, meta2, mycolors)

# save to file
pdf(paste0("hp_normed2_by_cell_type.pdf"), height=10, width = 10)
print(p2)
dev.off()

## heatmap of participants in G2 across cell types
C_target = 'G2'
g <- "subject.group2"
tag <- ""
clusterNo <- 5

file <- paste0('meta_kmeans_',tag,'_',clusterNo,'.csv')
meta3 <- read.csv(file, row.names = 1) %>% 
        mutate(preClust2 = ifelse(preClust %in% c("B memory", "B naive", "B intermediate"), "B cell", ifelse(preClust %in% c("CD16 Mono", "CD14 Mono"), "Monocytes", ifelse(preClust %in% c("NK", "NK_CD56bright"), "NK cell", preClust))),
                kmeans2 = ifelse(kmeans==C_target, kmeans, "other"))


selected_cell_states3 <- c('B cell', 'Monocytes', 'CD4 Naive', 'CD4 TCM', 'CD8 Naive', 'CD8 TEM', 'NK cell', 'Treg')

g2_cluster <- meta3 |> filter(RNA_cells>0) |> filter(preClust2 %in% selected_cell_states3) |> filter(kmeans=='G2') |> 
    group_by(subject.subjectGuid, preClust2) |> summarise(n=n())
data_p <- g2_cluster |> tidyr::pivot_wider(names_from = subject.subjectGuid, values_from = n, values_fill = 0) |> tibble::column_to_rownames('preClust2')

# manually add participant CU1044 for consistency, which doesn't have signature cluster
data_p <- cbind(data_p,data.frame('CU1044'=rep(0,nrow(data_p))))

# get the participant order
pt_orders <- read.csv(paste0(datadir,"heatmap_col_order.csv"), row.names = 1) |> mutate(new_sub=1:67) |> mutate(new_id=paste0('p', new_sub)) |> 
    tibble::rownames_to_column('id') |> tibble::column_to_rownames('new_id')

# change the participant name
colnames(data_p) <- rownames(pt_orders)[match(colnames(data_p), pt_orders$id)]

# reorganize it by participant order of `pt_orders`
data_p <- data_p[,intersect(rownames(pt_orders), colnames(data_p))]

# only non-CON groups with rows clustered
colors = list(subject.group3 = c('#E41A1C','#4DAF4A','#FFFF33'))
names(colors$subject.group3) = c("ERA","CON","At-Risk")

breaksList <- c(seq(-1,12,1),seq(17,22,5))
p1 <- pheatmap::pheatmap(data_p, filename = 'hp_cell_type_across_patients_row_clustered.pdf', 
                         cellheight = 7, fontsize = 7, cellwidth = 7, 
                         cluster_cols = F, cluster_rows = T,
                         breaks = breaksList,
                         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(length(breaksList)),
                         annotation_names_col = FALSE, 
                         annotation_col=pt_orders[,'subject.group3',drop=F], annotation_colors=colors
                        )

data_pm <- data_p
data_pm[data_pm > 10] <- 10

meta_ps <- pt_orders |> filter(!subject.group3=='CON')
data_ps <- data_pm[,intersect(names(data_pm),rownames(meta_ps))]
data_ps <- data_ps[rowSums(data_pm)>0,]

p2 <- pheatmap::pheatmap(data_ps, filename = 'hp_cell_type_across_patients_noCON_row_clustered.pdf', 
                         cellheight = 7, fontsize = 7, cellwidth = 7, 
                         cluster_cols = F, cluster_rows = T,
                         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(11),
                         # breaks = breaksList,
                         annotation_names_col = FALSE, 
                         annotation_col=meta_ps[,'subject.group3',drop=F], annotation_colors=colors
                        )



# save to file
pdf(paste0("hp_n_clusters_cell_type_across_patients_",tag,".pdf"), height=2, width = 10)
print(p1)
dev.off()

# save to file
pdf(paste0("hp_n_clusters_cell_type_across_patients_",tag,"_noCON.pdf"), height=2, width = 6)
print(p2)
dev.off()

## top signature TFs across cell types
## we selected top G2-specific TFs to visualize the mean PageRank scores
TFs <- readLines(paste0(resultdir,'/kmeans_5/G2_p0.01_d0.5_gene409_wilcoxon.txt'))[1:100]


info <- meta2 %>% mutate(kmeans2 = ifelse(kmeans=="G2", 'G2', "other"))
df <- pr_normed2[TFs,] |> tibble::rownames_to_column('TF') |> tidyr::pivot_longer(cols = !TF,names_to = 'id', values_to = 'value')
df1 <- info %>% dplyr::inner_join(df, by = "id") %>% dplyr::group_by(preClust, TF, kmeans2) %>% dplyr::summarise(mean.value=mean(value)) |> ungroup() |> 
    filter(preClust %in% selected_cell_states) |> dplyr::mutate(group=paste0(preClust,'_',kmeans2))


# second round filtering: diff(G2-other)>0.5
df2 <- df1 |> dplyr::group_by(TF, kmeans2) |> dplyr::summarise(mean2.value=mean(mean.value)) |> 
    tidyr::pivot_wider(names_from = kmeans2, values_from = mean2.value) |> mutate(diff=G2-other) 
TFs2 <- df2 |> filter(diff>0.5) |> pull(TF)

selected_cell_states4 <- c('B memory','CD14 Mono','CD4 Naive','CD8 Naive','CD4 TCM','CD8 TEM','MAIT','NK','Treg')

m <- df1 |> filter(TF %in% TFs2) |> select(TF, mean.value, group) |> tidyr::pivot_wider(names_from = group, values_from = mean.value) %>% tibble::column_to_rownames("TF")
m <- m[,c(paste0(selected_cell_states4,'_G2'), paste0(selected_cell_states4,'_other'))]


# set annotation 
annotation_col <- df1 |> distinct(kmeans2, group) |> tibble::column_to_rownames("group")

m[m > 2] <- 2
m[m < -2] <- -2

p <- pheatmap::pheatmap(m, 
                fontsize = 12, show_rownames = T,
                angle_col = 90, show_colnames = T,
                cluster_cols = F, cluster_rows = T, 
                cellwidth = 12, cellheight = 12,
                clustering_distance_cols = 'correlation', 
                clustering_distance_rows = 'correlation', 
                clustering_method = 'average',
                annotation_col = annotation_col, 
                annotation_colors = mycolors,
                border_color = NA,
                filename = NA)

# save to file
pdf(paste0("mean_pr_top_TFs_across_cell_types",tag,".pdf"), height=10, width = 10)
print(p)
dev.off()

print(sessionInfo())


