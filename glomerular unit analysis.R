# Glomerulr niche-based analysis###
Idents(s.sub) <- "cluster_name"
s.sub <- subset(s.sub, idents=c("NS"), invert=T)
DefaultAssay(s.sub) <- "Xenium"
s.sub <- NormalizeData(s.sub)
s.sub <- ScaleData(s.sub)
avg_expression <- AverageExpression(s.sub, group.by = "cluster_name", assays = "Xenium")
avg_expression1 <- AverageExpression(s.sub, group.by = c("cluster_name","sample"), assays = "Xenium")
avg_expression2 <- AverageExpression(s.sub, group.by = c("cluster_name","disease"), assays = "Xenium")
sample_info <- sub(".*_", "", colnames(avg_expression1$Xenium))
disease_info <- sub(".*_", "",  colnames(avg_expression2$Xenium))

gene_vars <- apply(avg_expression$Xenium, 1, var)
genes_to_keep <- names(sort(gene_vars, decreasing = TRUE)[1:300])
expr_matrix <- t(avg_expression$Xenium[genes_to_keep,])

#UMAP clustering
umap_result <- umap(expr_matrix, n_neighbors = 20, min_dist = 0.2)
umap_df <- data.frame(UMAP1 = umap_result$layout[,1],
                      UMAP2 = umap_result$layout[,2],
                      cluster = rownames(expr_matrix),
                      sample = sample_info,
                      disease = disease_info,
                      name=expr_matrix@Dimnames[[1]])

kmeans_result <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = 6, iter.max=5)
umap_df$cluster <- kmeans_result$cluster
name_map <- data.frame(old = c(1,2,3,4,5,6),new = c("Normal","Crescent", "Normal", "Crescent", "Normal", "Normal"))
umap_df$cluster <- name_map$new[match(umap_df$cluster, name_map$old)]
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "UMAP Clustering Result")

sds <- slingshot(data = umap_df[, c("UMAP1", "UMAP2")], 
                 clusterLabels = umap_df$cluster,
                 start.clus = "1", extend = "y", stretch = 0.8)

# Get all curves from the SlingshotDataSet
curves <- slingCurves(sds)
num_curves <- length(curves)

# Create the base plot
p <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("UMAP with Slingshot Trajectory")
p

# Add each curve to the plot
for(i in 1:num_curves) {
  curve_data <- curves[[i]]$s
  curve_df <- data.frame(UMAP1 = curve_data[,1], UMAP2 = curve_data[,2])
  p <- p + geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2), 
                     inherit.aes = FALSE, color = "black", linewidth = 0.5)
}

# Display the plot
p

curve_data <- curves[[1]]$s
curve_order <- curves[[1]]$ord
if (!is.null(curve_order)) {
  ordered_curve <- curve_data[curve_order, ]
} else {
  ordered_curve <- curve_data
}
curve_df <- data.frame(UMAP1 = ordered_curve[,1], UMAP2 = ordered_curve[,2])

p1 <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
  geom_point(size = 3) +
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2), 
            inherit.aes = FALSE, color = "red", linewidth = 1.5) +
  theme_minimal() +
  ggtitle("UMAP with Ordered Trajectory") + theme(panel.grid = element_blank())

p1


umap_df$pseudotime <- slingPseudotime(sds)[,1]
pt <- slingPseudotime(sds)
pt <- pt[rownames(umap_df), , drop = FALSE]

umap_df$pseudotime <- ifelse(is.na(pt[,1]), pt[,2], pt[,1])

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  ggtitle("Pseudotime trajectory from upper right to lower left") +
  geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2), 
            inherit.aes = FALSE, color = "red", linewidth = 1.5) +
  theme(panel.grid = element_blank())

curve_data <- slingCurves(lines)[[1]]$s 
curve_df <- data.frame(UMAP1 = curve_data[,1], UMAP2 = curve_data[,2])

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = pseudotime)) +
  geom_point(size = 3) +
  scale_color_viridis_c(option = "plasma") +
  theme_minimal() +
  ggtitle("Pseudotime trajectory from upper right to lower left") +
  #geom_path(data = curve_df, aes(x = UMAP1, y = UMAP2), 
  #          inherit.aes = FALSE, color = "red", linewidth = 1.5)+
  theme(panel.grid = element_blank())

#umap_df$pseudotime <-16 -umap_df$pseudotime#for glomerular
umap_df$name_clean <- sub("^(C\\d+[a-z]*).*", "\\1", umap_df$name)
umap_df_sorted <- umap_df[order(umap_df$name_clean),]
#umap_df_sorted$pseudotime <- -umap_df_sorted$pseudotime

library(splines)
library(pheatmap)

library(tidyr)
cluster_proportions <- s.sub@meta.data %>%
  dplyr::group_by(cluster_name, clusters) %>%
  dplyr::summarise(count = n(), .groups = 'drop') %>%
  dplyr::group_by(cluster_name) %>%
  dplyr::mutate(proportion = count / sum(count)) %>%
  dplyr::ungroup()

cluster_pseudotime <- data.frame(
  cluster_name = umap_df_sorted$name_clean,
  pseudotime = umap_df_sorted$pseudotime
)


cluster_proportions <- cluster_proportions %>%
  left_join(cluster_pseudotime, by = "cluster_name")

data_long <- cluster_proportions %>%
  pivot_longer(cols = c(proportion),
               names_to = "variable", 
               values_to = "value")

target_clusters <- c("Podocyte","Mesangium","EC_Vascular","Damaged", "Normal") # Select cells

cluster_colors <- c(
  "CD4T" = "#1f77b4",             # 青
  "CD8T" = "#ff7f0e",             # オレンジ
  "CD4Eff" = "#1f77b4",             # 青
  "CD8Eff" = "#ff7f0e",             # オレンジ
  "Proliferating_Lymp" = "#2ca02c",   # 緑
  "Normal" = "#2ca02c",   # 緑
  "Eosino/Baso" = "#2ca02c",   # 緑
  "Tcelldim" = "#d62728",         # 赤
  "Damaged" = "#d62728",         # 赤
  "TREM2_Macrophage" = "#9467bd", # 紫
  "Resident_Macrophage" = "#8c564b",  # 茶
  "Inflammatory_Macrophage" = "#e377c2",       # ピンク
  "Neutrophil" = "#e967c2",       # ピンク
  "Fibroblast_FAP" = "#7f7f7f",   # グレー
  "Fibroblast" = "#bcbd22",       # 黄緑
  "Podocyte" = "#17becf",         # 水色
  "Mesangium" = "#aec7e8",        # 薄青
  "PEC" = "#ffbb78",              # 薄オレンジ
  "EC_Vascular" = "#98df8a",      # 薄緑
  "Glomerular_EC" = "#98df8a",      # 薄緑
  "Proximal_Tubules" = "#c5b0d5", # 薄紫
  "Proximal_PCT" = "#c5b0d5", # 薄紫
  "Proximal_Damaged" = "#c49c94",  # 薄茶
  "Reg_PT" = "#c49c94",  # 薄茶
  "Bcell" = "#c45c32"  # 薄茶
)

cluster_proportions_filtered <- cluster_proportions %>%
  filter(clusters %in% target_clusters)

data_long_filtered <- cluster_proportions_filtered %>%
  pivot_longer(cols = c(proportion),
               names_to = "variable", 
               values_to = "value")


max_pseudotime <- max(data_long_filtered$pseudotime, na.rm = TRUE)
x_start <- 1
arrow_y <- 0
y_end=0.15

plot.png <- ggplot(data_long_filtered, aes(x = pseudotime, y = value, color = clusters, fill = clusters)) +
  geom_smooth(method = "loess", se = TRUE, span=1, alpha=0.15) +
  theme_classic() +
  coord_cartesian(xlim = c(x_start, max_pseudotime), ylim = c(0, y_end)) +
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    plot.title = element_text(size = 15)
  ) +
  scale_color_manual(values = cluster_colors) +
  scale_fill_manual(values = cluster_colors) +
  annotate("segment", x = x_start, xend = 7.9, y = arrow_y, yend = arrow_y,
           colour = "blue", size = 1,
           arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) +
  annotate("segment", x = 8.1, xend = 11.9, y = arrow_y, yend = arrow_y,
           colour = "red", size = 1,
           arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) +
  annotate("segment", x = 12.1, xend = max_pseudotime, y = arrow_y, yend = arrow_y,
           colour = "green4", size = 1,
           arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) +

  annotate("segment", x = 8, xend = 8, y = 0, yend = y_end, colour = "blue", linetype = "dashed", size = 0.5) +
  annotate("segment", x = 12, xend = 12, y = 0, yend = y_end, colour = "red", linetype = "dashed", size = 0.5) +
  annotate("segment", x = max_pseudotime, xend = max_pseudotime, y = 0, yend = y_end, colour = "green4", linetype = "dashed", size = 0.5)
plot.png
