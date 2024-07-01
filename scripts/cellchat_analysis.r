# basically follow the cellchat tutorial with minor customization
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(CellChat)
library(patchwork)

meta <- read.table("dataFreeze_20230310_filtered.csv", sep = ",", header = T) 
arrows <- list.files(path="cache/",pattern = ".arrow",recursive = T, full.names = T)
rnas <- list.files(path = 'cache/', pattern = "labeled.h5",recursive = T, full.names = T)
samples <- unique(meta$sample.sampleKitGuid)
print(samples)
print(length(samples))

meta2 <- read.csv('/home/jupyter/output_20230310/post-analysis/heatmap/meta_kmeans_5.csv') |> mutate(preClust2=ifelse(preClust %in% c("B memory", "B naive", "B intermediate"), "B cell", ifelse(preClust %in% c("CD16 Mono", "CD14 Mono"), "Monocytes", preClust))) |> filter(RNA_cells>0)
head(meta2)
dim(meta2)

get_seurat_object <- function(i){
	sample <- samples[i]
	arrow_file <- arrows[grepl(meta[which(meta$sample.sampleKitGuid==sample & grepl("arrow",meta$file.fileType)),"file.id"],arrows)]
	rna_file   <- rnas[grepl(meta[which(meta$sample.sampleKitGuid==sample & grepl("scRNA",meta$file.fileType)),"file.id"],rnas)]
	proj_name <- gsub("-.*","",gsub(".*/","",arrow_file))
	emb_file <- paste0(proj_name,"/RNAIntegration/GeneIntegrationMatrix/Save-Block1-JointCCA.rds")
	print("================================")
	print(paste0("sample: ", sample))
	print("================================")

	q <- LoadH5Seurat(sub("_labeled.h5","_labeled_v2.h5seurat",rna_file))

	itg <- readRDS(sub("-JointCCA","-JointCCA-cluster",emb_file))
	itg <- itg[grepl('reference',rownames(itg)),]
	rownames(itg) <- sub('_reference','',rownames(itg))

	# add co-embedded clusters label 
	q@meta.data$coemb.clusters <- paste0('C',itg[rownames(itg)[match(rownames(q@meta.data),rownames(itg))],'Co_clusters'])

	# only select clusters in final Taiji analysis
	cs <- meta2 |> filter(sample.sampleKitGuid==!!sample) |> pull(cluster)
	Idents(q) <- "coemb.clusters"  
	qs <- subset(q, cells = WhichCells(q, idents = cs))
	
	# save object
	file <- paste0("cellchat/",sample,"_rna_coembed.h5seurat")
	SaveH5Seurat(qs, file)
	print(paste0("Finished generating seurat object for cell chat: ", sample))

}


get_ccc_net <- function(i){
	sample <- samples[i]
	file <- paste0("cellchat/",sample,"_rna_coembed.h5seurat")

	qs <- LoadH5Seurat(file)

	cellchat <- createCellChat(object = qs, group.by = "ident", assay = "RNA")
	cellchat <- setIdent(cellchat, ident.use = "ident") 
	print(levels(cellchat@idents))
	groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
	print(groupSize)

	# set the ligand-receptor interaction database
	CellChatDB <- CellChatDB.human
	# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
	CellChatDB.use <- subsetDB(CellChatDB)
	cellchat@DB <- CellChatDB.use

	# subset the expression data of signaling genes for saving computation cost
	cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
	future::plan("multisession", workers = 4) # do parallel
	cellchat <- identifyOverExpressedGenes(cellchat)
	cellchat <- identifyOverExpressedInteractions(cellchat)

	# inference of cell-cell communication network
	cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 20)

	# filter out cell-cell communication if there're only few cells in certain cell groups
	cellchat <- filterCommunication(cellchat, min.cells = 10)

	# extract inferred cellular communication network as data frame
	df.net <- subsetCommunication(cellchat) # at the level of ligand/receptors

	# save to file
	write.csv(df.net, paste0("cellchat/",sample,"_all_ccc.csv"), row.names=F, quote=F)	

	# infer the cell-cell communication at a signaling pathway level
	cellchat <- computeCommunProbPathway(cellchat)

    # calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

    # save object 
	saveRDS(cellchat,paste0("cellchat/",sample,"_cellchat_obj.rds"))	
	print(paste0("Finished ccc analysis: ", sample))
}

o <- lapply(seq_along(samples), get_seurat_object)
# o <- lapply(c(44:67), get_seurat_object)

o <- lapply(seq_along(samples), get_ccc_net)



