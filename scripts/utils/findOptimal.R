#### findOptimal function is to determine the best distance method and cluster numbers ####
#### arguments: ####
##>>>>>>> 
library(factoextra)
#### main body ####
findOptimal <- function(data_reduced, output_file="distanceMetrics.pdf"){
  
  # elbow method
  fviz_nbclust(data_reduced, kmeans, method = "wss", k.max = 60)
  
  # silhouette method
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",k.max = 60)
  df2 <- g$data
  
  # use pearson distance
  res.dist <- get_dist(data_reduced, method = "pearson")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = 60)
  df2 <- rbind(df2,g$data)
  
  # use spearman distance
  res.dist <- get_dist(data_reduced, method = "spearman")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = 60)
  df2 <- rbind(df2,g$data)
  
  # use kendall distance
  res.dist <- get_dist(data_reduced, method = "kendall")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = 60)
  df2 <- rbind(df2,g$data)
  
  # use manhattan distance
  res.dist <- get_dist(data_reduced, method = "manhattan")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = 60)
  df2 <- rbind(df2,g$data)
  
  # visualization
  df2$method <- c(rep("euclidean",60),rep('pearson',60),rep('spearman',60),rep('kendall',60),rep('manhattan',60))
  df2[,'clusters'] <- as.numeric(as.character(df2[,'clusters']))
  df2 <- df2[df2$clusters>=1,]
  p2 <- ggplot(df2)+aes(x=clusters,y=y,color=method,group=method)+
    geom_line()+geom_point()+
    labs(x='number of clusters K',y='average silhouette width',color='distance metrics')+
    scale_x_continuous(breaks =seq(0,60,by=5))+
    # theme(axis.title = element_text(size = 15),
    #       axis.text = element_text(size = 15),
    #       legend.title = element_text(size = 15),
    #       legend.text = element_text(size = 15)
    #       )+
    theme_light()
  # png("distanceMetrics.png",units="in", width=6, height=5, res=300)
  pdf(output_file,width=8.5,height = 5)
  print(p2)
  dev.off()
}

