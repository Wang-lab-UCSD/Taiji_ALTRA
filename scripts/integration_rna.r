library(Seurat)
library(SeuratDisk)
library(dplyr)
library(H5weaver)
library(hise)
library(BPCells)

options(future.globals.maxSize = 80000 * 1024^2)

read_h5_seurat <- function(h5_file,
                           target = "matrix",
                           feature_names = "id",
                           ...) {
  
  if(!requireNamespace("Seurat", versionCheck = list(op = ">=", version = "3.1.0"))) {
    stop("Can't find the Seurat package. Please install with install.packages(\"Seurat\")")
  }
  
  assertthat::assert_that(is.character(h5_file))
  assertthat::assert_that(length(h5_file) == 1)
  assertthat::assert_that(typeof(feature_names) == "character")
  assertthat::assert_that(length(feature_names) == 1)
  assertthat::assert_that(feature_names %in% c("id","name"))
  
  mat <- read_h5_dgCMatrix(h5_file,
                           target = target,
                           feature_names = feature_names)
  
  cell_meta <- read_h5_cell_meta(h5_file,
                                 target = target)
  rownames(cell_meta) <- cell_meta$barcodes
  
  rownames(cell_meta) <- cell_meta$barcodes
  
  feat_meta <- read_h5_feature_meta(h5_file,
                                    target = target)
  
  cite <- FALSE
  
  # Check for CITE-seq data
  if("feature_type" %in% names(feat_meta)) {
    if("Antibody Capture" %in% feat_meta$feature_type) {
      cite <- TRUE
    }
  }
  
  if(cite) {
    cite_feat <- feat_meta[feat_meta$feature_type == "Antibody Capture",]
    feat_meta <- feat_meta[feat_meta$feature_type != "Antibody Capture",]
    
    cite_mat <- mat[cite_feat$id,]
    mat <- mat[feat_meta$id,]
  }
    
  so <- Seurat::CreateAssayObject(counts = mat,...) # added this line to force there's no duplicated colnames()
    
  so <- Seurat::CreateSeuratObject(counts = so,
                                   meta.data = cell_meta,
                                   ...)
  if(cite) {
    so[["ADT"]] <- Seurat::CreateAssayObject(counts = cite_mat)
  }
  
  so
}


meta <- read.csv("../dataFreeze_20241127_cross_section.csv", row.names=1)
dirs <- meta |> pull(file.id) |> unique()
rna_files <- unlist(lapply(dirs, function(x) list.files(path = paste0('../cache/',x),pattern = "labeled.h5", full.names = T)))
                           
                                                                                    
# use BPCells with Seurat objects
# Loop through h5ad files and output BPCells matrices on-disk
data.list <- c()
metadata.list <- c()

for (i in 1:length(rna_files)) {
    
    seurat_obj <- read_h5_seurat(rna_files[i])
    # write_matrix_dir(mat = seurat_obj[["RNA"]]$counts, 
    #                  dir = paste0(gsub(".h5", "", rna_files[i]), "_BP"))     # create BP cell folder

    # Load in BP matrices
    mat <- open_matrix_dir(dir = paste0(gsub(".h5", "", rna_files[i]), "_BP"))
    # Get metadata
    metadata.list[[i]] <- seurat_obj@meta.data
    data.list[[i]] <- mat
}
# Name layers
names(data.list) <- sub('-0[0-9].*','',sub('.*P[0-9]_PB','KT',rna_files))

# Add Metadata
for (i in 1:length(metadata.list)) {
  metadata.list[[i]]$sample.sampleKitGuid <- names(data.list)[i]
}

metadata <- Reduce(rbind, metadata.list)         
                           
# create Seurat object
merged.object <- CreateSeuratObject(counts = data.list, meta.data = metadata)
print(merged.object)
                           
# sample representative cells from each dataset
merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object, verbose = FALSE)
merged.object <- SketchData(object = merged.object, ncells = 2000, method = "LeverageScore", sketched.assay = "sketch")
print(merged.object)
saveRDS(object = merged.object, file = "merged_rna_67samples.rds")                           
saveRDS(object = merged.object@meta.data, file = "merged_rna_67samples_metadata.rds")
print("object saved!! Now performing UMAP analysis")   
                           
                           
# perform UMAP -- without integration                           
DefaultAssay(merged.object) <- "sketch"
merged.object <- FindVariableFeatures(merged.object, verbose = F)
merged.object <- ScaleData(merged.object, verbose = F)
merged.object <- RunPCA(merged.object, verbose = F)   
                           
merged.object <- FindNeighbors(merged.object, dims = 1:30, reduction = "pca")
merged.object <- FindClusters(merged.object, resolution = 2, cluster.name = "unintegrated_clusters")

merged.object <- RunUMAP(merged.object, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
p1 <- DimPlot(merged.object, reduction = "umap.unintegrated", group.by = c("pbmc_sample_id", "seurat_clusters","pool_id"))

          
# perform UMAP -- with integration
merged.object <- IntegrateLayers(merged.object, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca", dims = 1:30, k.anchor = 20, reference = which(Layers(merged.object, search = "data") %in% c("data.KT00443")),
    verbose = F)
                    
# cluster the integrated data
merged.object <- FindNeighbors(merged.object, reduction = "integrated.rpca", dims = 1:30)
merged.object <- FindClusters(merged.object, resolution = 2)
merged.object <- RunUMAP(merged.object, dims = 1:30, reduction = "integrated.rpca", reduction.name = "umap.integrated")    

p2 <- DimPlot(merged.object, reduction = "umap.integrated", group.by = c("pbmc_sample_id", "seurat_clusters","seurat_pbmc_type"))
                           
                           
# Integrate the full datasets
merged.object <- ProjectIntegration(object = merged.object, sketched.assay = "sketch", assay = "RNA", reduction = "integrated.rpca")
merged.object <- ProjectData(object = merged.object, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "integrated.rpca.full",
    full.reduction = "integrated.rpca.full", dims = 1:30, refdata = list(celltype.full = "seurat_pbmc_type"))

merged.object <- RunUMAP(merged.object, reduction = "integrated.rpca.full", dims = 1:30, reduction.name = "umap.full", reduction.key = "UMAP_full_")

p3 <- DimPlot(merged.object, reduction = "umap.full", group.by = c("pbmc_sample_id","celltype.full","pool_id"), alpha = 0.1)
                           
# save umap plot
pdf("Integration_rna_67samples.pdf", width=16, height=4)
print(p1)
print(p2)
print(p3)
dev.off()


# save object
print(format(object.size(merged.object), units = "Mb"))                           
saveRDS(object = merged.object, file = "merged_rna_67samples.rds")  
                           
                           
                           
## part 2: QC
merged.object <- readRDS("merged_rna_67samples.rds")
                           
### calculate percent.mt and percent.ribo
# Connect to AnnotationHub
library(AnnotationHub)
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# get the AnnotationHub ID the most recent database
id <- ahDb |> mcols() |> rownames() |> tail(n=1)                           
# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, return.type = "data.frame")                                                       # Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)
                                       
# Extract IDs for mitochondrial genes
mt <- annotations %>% dplyr::filter(seq_name == "MT") %>% dplyr::pull(gene_id)
ribo_ensg <- annotations |> dplyr::filter(grepl('^RP[SL]', gene_name)) |> pull(gene_id)

merged.object <- PercentageFeatureSet(merged.object, features = intersect(Features(merged.object), ribo_ensg), col.name = "percent.ribo", assay = 'RNA')
merged.object <- PercentageFeatureSet(merged.object, features = intersect(Features(merged.object), mt), col.name = "percent.mt", assay = 'RNA')

meta <- read.csv("meta_kmeans_pr_1031_5.csv", row.names = 1) |> distinct(sample.sampleKitGuid, subject.group3, sample)
meta2 <- read.csv("../Guid_KitGuid_67samples.csv") |> mutate(sample.new.id=sub('.*_','',title)) |> distinct(sample.sampleKitGuid, sample.new.id)

df <- merged.object@meta.data |> dplyr::select(barcodes, sample.sampleKitGuid, nCount_RNA, nFeature_RNA, n_reads, percent.mt, percent.ribo) |> left_join(meta, by='sample.sampleKitGuid') |> left_join(meta2, by='sample.sampleKitGuid') |> mutate(log10GenesPerUMI=log10(nFeature_RNA)/log10(nCount_RNA))


# add co_clusters information
itg <- readRDS('cell_co_clusters_all_67samples.RDS')
df <- df |> left_join(itg, by=c('barcodes','sample'))

meta <- read.csv("meta_kmeans_pr_1031_5.csv", row.names = 1) |> distinct(sample.sampleKitGuid, subject.group3, sample)
meta2 <- read.csv("../Guid_KitGuid_67samples.csv") |> mutate(sample.new.id=sub('.*_','',title)) |> distinct(sample.sampleKitGuid, sample.new.id)
meta3 <- read.csv("meta_kmeans_pr_1031_5.csv", row.names = 1) |> distinct(sample.sampleKitGuid, cluster, kmeans) |> dplyr::rename(Co_clusters=cluster)

itg <- readRDS('cell_co_clusters_all_67samples.RDS')

df <- merged.object@meta.data |> left_join(meta, by='sample.sampleKitGuid') |> left_join(itg, by=c('barcodes','sample'))|> left_join(meta3, by=c('sample.sampleKitGuid','Co_clusters'))|> left_join(meta2, by='sample.sampleKitGuid') |> mutate(log10GenesPerUMI=log10(nFeature_RNA)/log10(nCount_RNA))

merged.object@meta.data <- df
print('updated meta data!')
saveRDS(object = merged.object, file = "merged_rna_67samples.rds")  

# filter seurat object
filtered_seurat <- subset(x=merged.object, subset=(kmeans %in% paste0('G',1:5)))
print(filtered_seurat)

# save to object
# saveRDS(object = df, file = "merged_rna_67samples_metadata.rds")
saveRDS(object = filtered_seurat, file = "merged_rna_67samples_703701.rds")  

mycolors <- readRDS('mycolors_1031.rds')
# alpha=0.1
merged.object <- readRDS("merged_rna_67samples_703701.rds")
merged.object@meta.data$sample.new.id <- factor(merged.object@meta.data$sample.new.id, levels = paste0('p',1:67))
names(mycolors$subject.id) <- paste0('p',1:67)

# p1 <- DimPlot(merged.object, reduction = "umap.full", group.by = c("sample.new.id","pool_id")) & NoAxes()
# p1[[1]]$layers[[1]]$aes_params$alpha = alpha
# p1[[2]]$layers[[1]]$aes_params$alpha = alpha
# p2 <- DimPlot(merged.object, reduction = "umap.full", group.by = c("kmeans","Group"), cols=c(mycolors$cell.type, mycolors$kmeans)) & NoAxes()
p3 <- DimPlot(merged.object, reduction = "umap.full", group.by = c("sample.new.id","Group"), cols=c(mycolors$cell.type, mycolors$subject.id)) & NoAxes()

pdf("Integration_rna_67samples_703701_umap.full_v2.pdf", width=15, height=4)
# print(p1)
# print(p2)
print(p3)
dev.off()


# p1 <- DimPlot(merged.object, reduction = "umap.unintegrated", group.by = c("sample.new.id","pool_id")) & NoLegend() & NoAxes()
# p1[[1]]$layers[[1]]$aes_params$alpha = alpha
# p1[[2]]$layers[[1]]$aes_params$alpha = alpha
# p2 <- DimPlot(merged.object, reduction = "umap.unintegrated", group.by = c("kmeans","Group"), cols=c(mycolors$cell.type, mycolors$kmeans))& NoLegend() & NoAxes()

# pdf("Integration_rna_67samples_all_703701_umap.unintegrated.pdf", width=6, height=3)
# print(p1)
# print(p2)
# dev.off()


# p1 <- DimPlot(merged.object, reduction = "umap.integrated", group.by = c("sample.new.id","pool_id")) & NoLegend() & NoAxes()
# p1[[1]]$layers[[1]]$aes_params$alpha = alpha
# p1[[2]]$layers[[1]]$aes_params$alpha = alpha
# p2 <- DimPlot(merged.object, reduction = "umap.integrated", group.by = c("kmeans","Group"), cols=c(mycolors$cell.type, mycolors$kmeans))& NoLegend() & NoAxes()

# pdf("Integration_rna_67samples_all_703701_umap.integrated.pdf", width=6, height=3)
# print(p1)
# print(p2)
# dev.off()










