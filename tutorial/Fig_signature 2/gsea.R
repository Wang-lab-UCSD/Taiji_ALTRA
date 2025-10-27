library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)

gsea <- function(x,key=NULL,is.table=T,is.gene.list=T,output_dir="./", showCategory=10, go.width=5, go.height=5){
  
  # get genes and gene id -------------------
  prefix <- gsub(".*/","",gsub(".txt","",x))
  gene <- readLines(x)
  if (length(gene) >= 1){
      df <- AnnotationDbi::select(org.Hs.eg.db, gene, "ENTREZID", "SYMBOL") %>% tidyr::drop_na()
      gene.id <- df[,2]
      if (is.gene.list){
          df[,1] %>% writeLines(paste0(output_dir,prefix,"_geneList.txt"))
      }

      # GO term over-representation test-------------------
      ego <- enrichGO(gene = gene, 
                      OrgDb = org.Hs.eg.db, ont = "BP",keyType = "SYMBOL",
                      pvalueCutoff = 0.1,qvalueCutoff = 0.5
      )
      ego2 <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)

      ## select specific GO terms----------------
      if (!is.null(key)){
        key.term <- which(ego$ID%in%key)
        if (length(key.term)>=1){
            o <- ego[key.term]
            write.table(o, paste0(output_dir,prefix,"_specificGOterms.txt"),quote=F, row.names = F)
        }
      }

      if(!is.null(ego2) && nrow(ego2)>2){
        p1 <- dotplot(ego2, showCategory=showCategory)+ggtitle("GO")
        pdf(file = paste0(output_dir,prefix,"_go.pdf"),width = go.width,height = go.height)
        print(p1)
        dev.off()
      }

      # KEGG over-representation test ---------------------------
      kk <- enrichKEGG(gene = gene.id, organism = "hsa", 
                       minGSSize = 2, pvalueCutoff = 0.1, qvalueCutoff = 0.5)
      if (!is.null(kk) && nrow(kk)>2){
        p2 <- dotplot(kk, showCategory=showCategory)+ggtitle("KEGG")
        pdf(file = paste0(output_dir,prefix,"_kegg.pdf"),width = go.width,height = go.height)
        # p <- cowplot::plot_grid(p1, p2, ncol = 2)
        print(p2)
        dev.off()
      }

      # Reactome over-representation analysis ----------------------------
      x <- enrichPathway(gene=gene.id, pvalueCutoff = 0.1, readable=TRUE)
      if(!is.null(x) && nrow(x)>2){
        p3 <- dotplot(x, showCategory=showCategory)+ggtitle("Reactome")
        pdf(file = paste0(output_dir,prefix,"_reactome.pdf"),width = go.width,height = go.height)
        print(p3)
        dev.off()
      }


      # output tables or not -------------------------------------
      if(is.table){
        kk2 <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
        cols <- intersect(intersect(names(ego2@result), names(kk2@result)), names(x@result))
        output <- rbind(ego2@result[,cols],kk2@result[,cols],x@result[,cols])
        write.table(output,paste0(output_dir,prefix,"_GO.csv"),row.names = F, sep = ",")
      }
  }
}
