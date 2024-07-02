## create new sub-dir and set the new sub-dir as the working dir
createDir <- function(subdir){
    dir.create(file.path(maindir, subdir), showWarnings = FALSE)
    setwd(file.path(maindir, subdir))
    print(getwd())
}

## normalize data across rows to make sure the values is in range [0,1] after norm
scaleData <- function(df){
  return(as.data.frame(t(apply(as.matrix(df), MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X))))))
}
                               
## get z-scores of data 
zscore <- function(df){
  return(as.data.frame(t(scale(t(as.matrix(df))))))
}                      
                           
## calculate covariance, x is the matrix
co.var <- function(x) ( matrixStats::rowSds(x) / matrixStats::rowMeans2(x) )
    
## load gene expression data, need to specify the dir 
## clean the data by removing columns with colSum==0 and rows with rowSum==0
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
                               
##load pagerank data -------------------------
# remove DNMT1 since it's not a TF
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