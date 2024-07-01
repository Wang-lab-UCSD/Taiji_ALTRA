# set up
# args = commandArgs(trailingOnly=TRUE)
# p = as.numeric(args[1])
# cell_type = args[2]

library(tidyverse)
odir <- "/home/jupyter/inputs/"
dir.create(file.path(odir),showWarnings = F)
setwd(odir)

rna.fl <- readLines("rna.fl")
atac.fl <- readLines("atac.fl")
# time <- "year1_day7"
# meta.fl <- paste0("/home/jupyter/metadata_",time,".csv")
time <- "20230310"
meta.fl <- "/home/jupyter/dataFreeze_20230310_filtered.csv"


rna.fl.sample <- lapply(rna.fl, function(x) sub("_rna.tsv","",x))
atac.fl.sample <- lapply(atac.fl, function(x) sub("_atac.bed.gz","",x))
common.sample <- lapply(intersect(atac.fl.sample, rna.fl.sample), function(x) sub(".*/","",x))
                        

# create more comprehensive meta data by incorporating subject information
meta_files <- list.files(path = "/home/jupyter/cluster_meta",pattern = "*csv",full.names = T)
read_plus <- function(x){
    df <- read.csv(x)
    df$sample <- gsub("_cluster.*","",gsub(".*/","",x))
    return(df)
}
meta.cluster <- do.call(rbind,lapply(meta_files,read_plus))
meta.sample <- read.csv(meta.fl) %>% 
                dplyr::mutate(sample=gsub("-.*","",gsub(".*/","",file.name))) %>%
                dplyr::select(sample.visitName, subject.biologicalSex, subject.partnerCode, subject.subjectGuid,cohort.cohortGuid,subject.age,sample)
meta <- meta.cluster %>% dplyr::mutate(id = paste0(sample,"_",cluster)) %>%
                        dplyr::inner_join(meta.sample, by = "sample") %>%
                        dplyr::filter(id %in% common.sample)
write.table(meta, paste0("../metadata_",time,"_expand.csv"), quote = F, row.names = F, sep = ",")
                        
# main function                        
main <- function(filter.type=NULL,cell.iden=NULL,min.purity=NULL){
    # filter.type = c("ident","purity","both")
    # cell.iden = unique(meta$preClust)
    # min.purity belongs [0.1]
    
    rnaseq = "RNA-seq:"
    atacseq = "ATAC-seq:"
    hicseq = "HiC:"
    
    # get meta based on different filters
    if (is.null(filter.type)){
        meta2 <- meta
        output_file <- paste0("input_",time,".yml")
    }else if(filter.type=="ident"){
        meta2 <- meta[grepl(cell.iden,meta$preClust),]
        output_file <- paste0("input_",time,"_",gsub("-|\\s+","",cell.iden),".yml")
    }else if (filter.type=="purity"){
        meta2 <- meta[which(meta$purity>=min.purity),]
        output_file <- paste0("input_",time,"_",gsub("\\.","",min.purity),".yml")
    }else if (filter.type=="both"){
        meta2 <- meta[meta$purity>=min.purity & grepl(cell.iden,meta$preClust),]
        output_file <- paste0("input_",time,"_",gsub("-|\\s+","",cell.iden),"_",gsub("\\.","",min.purity),".yml")
    }
    meta2$id <- paste0(meta2$sample,"_",meta2$cluster)
    
    # loop and add
    for(i in seq_along(rna.fl)){

      i_rna = rna.fl[i]
      tis = gsub("_rna.tsv","",gsub(".*/","",i_rna))

      # check if tis passed the filter
      if (tis %in% meta2$id){
          i_atac = paste0("pseudobulk/",tis,"_atac.bed.gz")

          if(i_atac %in% atac.fl){
            line_tis1 = paste0("  - id: ", tis, "_RNA")
            line_tis1_a = paste0("  - id: ", tis, "_ATAC")
            line_tis1_b = paste0("  - id: ", tis, "_HiC")
            line_tis2 = paste0("    group: ", tis)
            line_tis3 = "    replicates:"
            line_tis4 = "      - rep: 1"
            line_tis5 = "        files:"
            line_tis6 = paste0("          - path: /home/jupyter/",i_rna)
            line_tis6_a = paste0("          - path: /home/jupyter/",i_atac)
            line_tis6_b = "          - path: /home/jupyter/epitensor_loop_top10p_87311.txt"
            line_tis7 = "            tags: ['GeneQuant']"
            line_tis7_a = "            tags: ['PairedEnd']"
            line_tis7_b = "            tags: ['ChromosomeLoop']"

            rnaseq = c(rnaseq, line_tis1, line_tis2, line_tis3, line_tis4, line_tis5, line_tis6, line_tis7)
            atacseq = c(atacseq, line_tis1_a, line_tis2, line_tis3, line_tis4, line_tis5, line_tis6_a, line_tis7_a)
            hicseq = c(hicseq, line_tis1_b, line_tis2, line_tis3, line_tis4, line_tis5, line_tis6_b, line_tis7_b)
          }
      }
    }

    combine_line = c(rnaseq, atacseq, hicseq)
    write.table(combine_line,output_file,row.names = F,col.names=F,quote=F)
}

main()
# main(filter.type="ident",cell.iden=cell_type)
# main(filter.type="purity",min.purity=p)
# main(filter.type="both",cell.iden=cell_type,min.purity=p)


