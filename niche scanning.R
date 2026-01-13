#Cluster scanning using dbscan##
xenium.obj$x <- centroid_coords$x[match(rownames(xenium.obj@meta.data), rownames(centroid_coords))]
xenium.obj$y <- centroid_coords$y[match(rownames(xenium.obj@meta.data), rownames(centroid_coords))]
Idents(xenium.obj) <- "SCT_snn_res.0.1"
s.sub <- subset(xenium.obj, idents=c("Proximal_Tubules", "Reg_Proximal", "STC"))

coords_matrix <- cbind(
  x = s.sub$x,
  y = s.sub$y,
  cluster = s.sub$clusters
)

coords_xy <- as.matrix(coords_matrix[, c("x", "y")])
coords_xy <- matrix(as.numeric(coords_xy), ncol = 2)
coords_xy <- na.omit(coords_xy)

dbscan_result <- hdbscan(coords_xy, minPts = 10)
clusters <- dbscan_result$cluster

cluster_sizes <- table(clusters[clusters > 0])

small_clusters <- names(cluster_sizes[cluster_sizes <= 10])
clusters[clusters %in% small_clusters] <- -1 
cluster_names <- ifelse(clusters == -1, "NS", paste0("C", clusters))
cluster_names[clusters == 0] <- "Noise"

s.sub$cluster_name <- factor(cluster_names)
s.sub$cluster_name <- as.character(s.sub$cluster_name)

split_cluster <- function(data, cluster_name, k) {
  Idents(data) <- data$cluster_name
  cluster_cells <- WhichCells(data, idents = cluster_name)
  cluster_data <- data@meta.data[cluster_cells, c("x", "y")]
  kmeans_result <- kmeans(cluster_data, centers = k)
  new_clusters <- paste0(cluster_name, letters[1:k])[kmeans_result$cluster]

  data$new_cluster_name <- data$cluster_name
  data$new_cluster_name[cluster_cells] <- new_clusters

  data$cluster_name <- data$new_cluster_name
  
  return(data)
}

s.sub <- split_cluster(s.sub, "C3", 2)
s.sub <- split_cluster(s.sub, "C8", 3)
