library(factoextra)
library(parallel)

findOptimal <- function(data_reduced, max_k = 60, cores = 1) {
  if (!is.numeric(max_k) || max_k <= 0 || max_k != round(max_k)) {
    stop("max_k must be a positive integer")
  }
  
  # Elbow method
  elbow_plot <- fviz_nbclust(data_reduced, kmeans, method = "wss", k.max = max_k)
  
  # Gap statistic method
  gap_plot <- fviz_nbclust(data_reduced, kmeans, method = "gap_stat", k.max = max_k, nboot=500)

  # Function to compute silhouette for a given distance method
  compute_silhouette <- function(method) {
    res.dist <- get_dist(data_reduced, method = method)
    g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette", diss = res.dist, k.max = max_k)
    g$data$method <- method
    return(g$data)
  }
  
  # Parallel computation for different distance methods
  methods <- c("euclidean", "pearson", "spearman", "kendall", "manhattan")
  results <- mclapply(methods, compute_silhouette, mc.cores = cores)
  
  df2 <- do.call(rbind, results)
  df2$clusters <- as.numeric(as.character(df2$clusters))
  df2 <- df2[df2$clusters >= 1, ]
  
  # Create combined plot
  combined_plot <- ggplot(df2, aes(x = clusters, y = y, color = method, group = method)) +
    geom_line() +
    geom_point() +
    labs(title = "Silhouette Scores for Different Distance Metrics",
         x = "Number of groups K", y = "Silhouette Score",color='distance metrics')
    
  return(list(data = df2, elbow_plot = elbow_plot, silhouette_plot = combined_plot, gap_plot = gap_plot))
}
