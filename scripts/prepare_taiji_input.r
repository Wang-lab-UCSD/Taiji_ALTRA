.libPaths()
.libPaths("/home/jupyter/R/lib")

# devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.refGene")
# BiocManager::install('scDataviz')
# remotes::install_github("mojaveazure/seurat-disk")

library(ArchR)
library(tidyverse)
library(hise)
library(Seurat)
library(SeuratDisk)
library(H5weaver)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
library(scDataviz)

# set up hyper-params
minCell=20
minPeak=200

# These functions allow us to extract the fragment positions directly from the .arrows:
chr_fragments_datatable <- function(h5_con, chr_name) {
    rg <- h5read(h5_con, paste0("/Fragments/",chr_name,"/Ranges"))

    data.table(
        chr = chr_name,
        start = rg[,1],
        end = rg[,1] + rg[,2],
        bc = rep(h5read(h5_con, paste0("/Fragments/",chr_name,"/RGValues")),
                 h5read(h5_con, paste0("/Fragments/",chr_name,"/RGLengths"))),
        n_umi = 1
    )
}

extract_arrow_fragments <- function(arrow_file) {
    h5_con <- H5Fopen(arrow_file)
    contents <- h5ls(h5_con)
    chrs <- contents$name[contents$group == "/Fragments"]

    frags <- chr_fragments_datatable(h5_con, chrs[1])
    for(i in 2:length(chrs)) {
        frags <- rbind(frags, chr_fragments_datatable(h5_con, chrs[i]))
    }

    setorderv(frags, c("chr","start"))

    frags
}

convert_arrow <- function(arrow_file) {
    print(paste("Extracting fragments from",arrow_file))
    fragment_dt <- extract_arrow_fragments(arrow_file)
    fwrite(fragment_dt,
            sub("_archr.arrow","_fragments.tar.gz",arrow_file),
            col.names = FALSE,
            sep = "\t")
}

# main function
main <- function(i,min.cell=minCell){
    sample <- samples[i]
    arrow_file <- arrows2[grepl(meta[which(meta$sample.sampleKitGuid==sample & grepl("arrow",meta$file.fileType)),"file.id"],arrows2)]
    rna_file   <- rnas[grepl(meta[which(meta$sample.sampleKitGuid==sample & grepl("scRNA",meta$file.fileType)),"file.id"],rnas)]
    proj_name <- gsub("-.*","",gsub(".*/","",arrow_file))
    emb_file <- paste0(proj_name,"/RNAIntegration/GeneIntegrationMatrix/Save-Block1-JointCCA.rds")
    print("================================")
    print(paste0("sample: ", sample))
    print("================================")
   
    # 1. process scRNA-seq and scATAC-seq cells respectively----------------------------------------
    if (!file.exists(sub("_labeled.h5","_labeled_v2.h5seurat",rna_file))){
        
        # create ArchR project
        RefGeneAnnotations <- createGeneAnnotation(TxDb= TxDb.Hsapiens.UCSC.hg38.refGene, OrgDb = org.Hs.eg.db)
        addArchRGenome('hg38')
        proj <- ArchRProject(arrow_file, outputDirectory = proj_name, geneAnnotation =  RefGeneAnnotations, copyArrows = F)
        
        # filter doublets 
        proj <- filterDoublets(proj, filterRatio=0.5)
        saveArchRProject(proj)

        addArchRThreads(5)

        ## Needs ArchR 1.0.2 to run. 
        proj<- addIterativeLSI(
          ArchRProj = proj,
          useMatrix = "TileMatrix", 
          name = "IterativeLSI", 
          iterations = 4,
          varFeatures = 25000,
          sampleCellsPre = 10000,
          projectCellsPre = FALSE,
          sampleCellsFinal = 100000,
          force = TRUE
        )

        proj <- addClusters(
          input = proj,
          reducedDims = "IterativeLSI",
          method = "Seurat",
          name = "Clusters",
          resolution = 3,
          force = TRUE
        )

        proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI", force = TRUE)
        saveArchRProject(proj)
        #Up to here, should only take 4 minutes

        ## load seurat reference to label the scRNA-cells
        reference <- LoadH5Seurat("./pbmc_multimodal.h5seurat")
        ## remove platelets from reference
        Idents(reference) <- 'celltype.l2'
        reference <- subset(x=reference, idents = 'Platelet', invert=TRUE)

        q <- read_h5_seurat(rna_file, feature_names = "name")
        q <- SCTransform(q, verbose = FALSE)

        anchors <- FindTransferAnchors(
            reference = reference,
            query = q,
            normalization.method = "SCT",
            reference.assay = "SCT",
            dims = 1:50,
            reference.reduction = "spca",
            recompute.residuals = FALSE
        )

        q <- TransferData(
          anchorset = anchors, 
          reference = reference,
          query = q,
          refdata = list(
            celltype.l1 = "celltype.l1",
            celltype.l2 = "celltype.l2",
            celltype.l3 = "celltype.l3")
        )
        q <- IntegrateEmbeddings(
          anchorset = anchors,
          reference = reference,
          query = q, 
          new.reduction.name = "ref.spca"
        )

        DefaultAssay(q) <- "RNA"

        # DimPlot(object = q, reduction = "ref.spca", group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
        SaveH5Seurat(q, sub("_labeled.h5","_labeled_v2.h5seurat",rna_file))

        ## Cell type labeling
        ## make the sample size large enough to avoid sampling
        ## tune dimsToUse because the co-cluster quality seems bad, and it doesn't improve
        proj <- addGeneIntegrationMatrix(
          ArchRProj = proj,
          useMatrix = "GeneScoreMatrix", 
          matrixName = "GeneIntegrationMatrix", 
          reducedDims = "IterativeLSI", 
          seRNA = q,
          addToArrow = FALSE,
          sampleCellsATAC = 50000,
          sampleCellsRNA = 50000,
          # dimsToUse = 1:50,
          # dims = 1:50,
          groupRNA = "predicted.celltype.l2",
          nameCell = "predictedCell_Un", 
          nameGroup = "predictedGroup_Un", 
          nameScore = "predictedScore_Un", 
          force=TRUE
        )

        saveArchRProject(proj)
        print(paste0(length(unique(proj@cellColData$predictedCell_Un))," scRNA cells are used to match scATAC cells"))
        # up to here, should take 15 min
        
    }else{
        print(paste0("already processed scRNA-seq and scATAC-seq files for ", sample, ". Now load the files and get the co-embedded clusters!"))
        q <- LoadH5Seurat(sub("_labeled.h5","_labeled_v2.h5seurat",rna_file))
        proj <- loadArchRProject(proj_name)
    }
    
    # 2. coembed scRNA-seq and scATAC-seq cells together---------------------------
    if (!file.exists(paste0("cluster_meta/",proj_name,"_cluster_identity.csv"))){
        ## load integration file
        itg <- readRDS(emb_file)
        head(itg)
        
        set.seed(10)
        print("integration of all cells from scATAC and scRNA")
        clusters <- scDataviz::clusKNN(as.matrix(itg[,1:(ncol(itg)-5)])) # wrapper function for Seurat's `FindNeighbors` and `FindClusters`
        # the process would break down if setting the resolution below 1
        # takes 10 min

        ## save clustering result
        names(clusters) <- rownames(itg)
        head(clusters)
        itg$Co_clusters <- clusters
        saveRDS(itg,sub("-JointCCA","-JointCCA-cluster",emb_file))
        print("finished writing integration results")
        print(table(itg$Assay))

        ## create cluster meta file
        cM <- as.matrix(confusionMatrix(itg$Co_clusters, itg$Group))
        preClust <- colnames(cM)[apply(cM, 1 , which.max)]
        cluster <- paste0("C",rownames(cM))
        purity <- apply(cM, 1 , function(x) round(max(x)/sum(x),3))
        cells <- apply(cM, 1, sum)                
        RNA_cells <- unlist(lapply(unique(itg$Co_clusters), function(x) sum(itg$Co_clusters == x & itg$Assay == "RNA")))
        names(RNA_cells) <- unique(itg$Co_clusters)
        ATAC_cells <- unlist(lapply(unique(itg$Co_clusters), function(x) sum(itg$Co_clusters == x & itg$Assay == "ATAC")))
        names(ATAC_cells) <- unique(itg$Co_clusters)
        meta.cluster <- as.data.frame(cbind(preClust,cluster,purity,cells,RNA_cells,ATAC_cells)) %>%
                        dplyr::mutate_if(colnames(.) %in% c("purity","cells","RNA_cells","ATAC_cells"), as.numeric) %>%
                        dplyr::arrange(-cells)
        
        filter <- meta.cluster$RNA_cells>=min.cell & meta.cluster$ATAC_cells>=min.cell
        print(meta.cluster[filter,])
        print(sum(filter))
                                    
        ## write to cluster file
        dir.create(file.path("cluster_meta"),showWarnings = F)
        file.path <- paste0("cluster_meta/",proj_name,"_cluster_identity.csv")
        write.table(meta.cluster, file.path, quote = F, row.names=F, sep=",")

        ## add to archr project and get plots
        cluster_atac <- clusters[grepl("_query",names(clusters))]
        names(cluster_atac) <- gsub("_query","",names(cluster_atac))

        head(cluster_atac)

        proj <- addCellColData(ArchRProj = proj, data = paste0("C",cluster_atac),cells = names(cluster_atac), name = "Co_clusters", force = T)
        saveArchRProject(proj)

        head(proj@cellColData)

        p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedGroup_Un")
        p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters")
        p3 <- plotEmbedding(proj, colorBy = "cellColData", name = "DoubletEnrichment")
        p4 <- plotEmbedding(proj, colorBy = "cellColData", name = "Co_clusters")
        p5 <- plotEmbedding(proj, colorBy = "cellColData", name = "predictedScore_Un")

        plotPDF(p1,p2,p3,p4,p5, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
        
        ## visualize the co-embedded clusters
        ## umap colored by assay, cell type and clusters
        set.seed(1) # Always do this prior to UMAP
        UMAPParams <- list(n_neighbors = 40, min_dist = 0.2, metric="cosine", verbose=FALSE)
        UMAPParams$X <- as.data.frame(itg[, grep("CC_", colnames(itg))])
        UMAPParams$ret_nn <- FALSE
        UMAPParams$ret_model <- FALSE
        UMAPParams$n_threads <- 1
        uwotUmap <- do.call(uwot::umap, UMAPParams)
                                    
        prefix <- gsub("Save-","",sub("-Joint.*","",gsub(".*/","",emb_file)))
        p1 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = itg$Assay,
          randomize = TRUE, 
          size = 0.2,
          title = paste0(prefix, " colored by Assay"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())
        p2 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = itg$Co_clusters,
          randomize = TRUE, 
          size = 0.2,
          title = paste0(prefix, " colored by Assay"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                 legend.position = "none")
        p3 <- ggPoint(
          x = uwotUmap[,1], 
          y = uwotUmap[,2], 
          color = itg$Group,
          randomize = TRUE, 
          size = 0.2,
          title = paste0(prefix, " colored by Assay"),
          xlabel = "UMAP Dimension 1",
          ylabel = "UMAP Dimension 2",
          rastr = TRUE
        )+ theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                 axis.text.y = element_blank(), axis.ticks.y = element_blank())
                                    
        pdf(paste0(proj_name,"/Plots/Plot-UMAP-",prefix,"_RNA_ATAC.pdf"))
        print(p1)
        print(p2)
        print(p3)
        dev.off()
        print(paste0("finished getting co-embebed clusters for ", sample, "!"))
    
    }else{
        print(paste0("already got co-embedded clusters for ", sample, ". Now get the pseudobulk!"))
        itg <- readRDS(sub("-JointCCA","-JointCCA-cluster",emb_file))    
    }
                                    
    # 3. prepare pseudo-bulk input for Taiji 
    ## combine scRNA cells
    print("----------------------------------------------")
    print(paste0("Now combine scRNA cells ", sample, "!"))
    mat <- q@assays$RNA@counts

    getCounts <- function(i){
      c1 = gsub("_reference","",rownames(itg[which(itg$Co_clusters==i),]))
      df1 <- mat[,which(colnames(mat) %in% c1)]
      if (!is.null(ncol(df1))){
          if (ncol(df1)>=min.cell){
              df2 <- rowSums(df1) 
              df3 <- data.frame("geneID" = names(df2),"count" = unname(df2))
              dir.create(file.path("./pseudobulk/"),showWarnings = F)
              write.table(df3, file = paste0("./pseudobulk/",proj_name,"_C",i,"_rna.tsv"), 
                          quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
          }
          print(paste0("number of scRNA cells in cluster ",i,": ", ncol(df1)))
      }else{print(paste0("number of scRNA cells in cluster ",i,": 0!"))}
    }

    c <- lapply(unique(itg$Co_clusters),getCounts)


    # combine scATAC cells
    ## write to fragments.tar.gz in the same directories in the cache.
    ## Should take ~2 minutes per file
    convert_arrow(arrow_file)

    ## get original matrix 
    file.path <- gsub("_archr.arrow","_fragments.tar.gz",arrow_file)
    mat <- fread(file.path, header = FALSE)

    c1 = gsub("_query","",rownames(itg[which(itg$Co_clusters=="1"),]))

    str(mat)
    # proj@cellColData
    c1 = gsub(".*#","",c1)
    tail(c1)

    # gather peaks in the same cluster
    getPeaks <- function(i,min.peak=minPeak){
      c1 = gsub("_query","",rownames(itg[which(itg$Co_clusters==i),]))
      c1 = gsub(".*#","",c1)
      df1 <- mat[which(mat[[4]]%in%c1),1:4]

      if (!is.null(nrow(df1))){
          n_cell = nrow(unique(df1[,4]))
          print(paste0("number of scATAC cells in cluster ",i,": ", n_cell))
          if (nrow(df1)>=min.peak && n_cell>=min.cell){
              print(paste0("number of scATAC peaks in cluster ",i,": ", nrow(df1)))
              write.table(df1, file = paste0("./pseudobulk/",proj_name,"_C",i,"_atac.bed"),quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
          }
      }else{print(paste0("number of scATAC cells in cluster ",i," is zero!"))}
    }

    c <- lapply(unique(itg$Co_clusters),getPeaks)
    print("-----------------------------------")
    print(paste0("success for sample ",sample))
}
                                    
# read meta data
# meta <- read.csv("20x20-v5.csv")
meta <- read.table("dataFreeze_20230310_filtered.csv", sep = ",", header = T)                            
arrows <- list.files(pattern = ".arrow",recursive = T)
arrows2 <- arrows[grepl("cache",arrows)]
rnas <- list.files(pattern = "labeled.h5",recursive = T)
samples <- unique(meta$sample.sampleKitGuid)
print(samples)
print(length(samples))
                                    
# run all samples                                    
o <- lapply(seq_along(samples),main)
                                    
# only deal with additional CU controls
# o <- lapply(c(46:67),main)   
# o <- lapply(c(38,43,44,48,49,50,52,53,54),main)   
                               
                      
# main(38)
sessionInfo()


