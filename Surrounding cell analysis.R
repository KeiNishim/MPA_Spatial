# Niche-surrounding cell analysis###
xenium.obj@meta.data$x <- suppressWarnings(as.numeric(centroid_coords$x))
xenium.obj@meta.data$y <- suppressWarnings(as.numeric(centroid_coords$y))
xenium.obj$cluster_name <- s.sub$cluster_name
 
calculate_cluster_info_safe <- function(xenium.obj, min_cells = 3) {
  df <- xenium.obj@meta.data %>%
    filter(!is.na(cluster_name), !is.na(x), !is.na(y))
  if (nrow(df) == 0) return(tibble())
  
  df %>%
    group_by(cluster_name) %>%
    summarise(
      n = n(),
      centroid_x = mean(x, na.rm = TRUE),
      centroid_y = mean(y, na.rm = TRUE),
      radius = {
        cx <- mean(x, na.rm = TRUE)
        cy <- mean(y, na.rm = TRUE)
        d <- sqrt((x - cx)^2 + (y - cy)^2)
        stats::quantile(d, 0.95, na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    filter(n >= min_cells, !is.na(radius), !is.na(centroid_x), !is.na(centroid_y))
}

find_surrounding_cells_safe <- function(cluster_info_row, obj) {
  if (any(is.na(c(cluster_info_row$centroid_x, cluster_info_row$centroid_y, cluster_info_row$radius)))) {
    return(tibble())
  }
  df <- obj@meta.data %>% filter(!is.na(x), !is.na(y))
  
  mutate(df, distance = sqrt((x - cluster_info_row$centroid_x)^2 + (y - cluster_info_row$centroid_y)^2)) %>%
    filter(distance > cluster_info_row$radius & distance <= 2 * cluster_info_row$radius)
}

calculate_surrounding_ratios_clusters <- function(xenium.obj) {
  ci <- calculate_cluster_info_safe(xenium.obj)
  if (nrow(ci) == 0) return(tibble())
  
  if (!"clusters" %in% colnames(xenium.obj@meta.data)) {
    stop("meta.data に 'clusters' 列がありません。")
  }
  
  results <- lapply(seq_len(nrow(ci)), function(i) {
    info <- ci[i, ]
    sc <- find_surrounding_cells_safe(info, xenium.obj)
    
    if (nrow(sc) == 0) {
      return(tibble(
        target_cluster = as.character(info$cluster_name),
        clusters = NA_character_,
        count = 0L,
        ratio = NA_real_
      ))
    }
    
    counts <- sc %>%
      group_by(clusters) %>%
      summarise(count = n(), .groups = "drop") %>%
      mutate(ratio = count / sum(count))
    
    tibble(
      target_cluster = as.character(info$cluster_name),
      clusters = counts$clusters,
      count = counts$count,
      ratio = counts$ratio
    )
  })
  
  bind_rows(results)
}

surrounding_ratios <- calculate_surrounding_ratios_clusters(xenium.obj)
surrounding_ratios_wide <- surrounding_ratios %>%
  mutate(target_cluster = as.character(target_cluster),
         clusters = as.character(clusters)) %>%
  replace_na(list(ratio = 0)) %>%
  pivot_wider(
    id_cols = target_cluster,
    names_from = clusters,
    values_from = ratio,
    values_fill = 0
  )

rownames(umap_df) <- sub("_.*", "", rownames(umap_df))

surrounding_ratios_wide <- surrounding_ratios %>%
  mutate(target_cluster = replace_na(target_cluster, "Unknown")) %>%
  pivot_wider(id_cols = target_cluster, 
              names_from = clusters, 
              values_from = ratio, 
              values_fill = 0)

surrounding_ratios_wide <- surrounding_ratios_wide %>%
  rename_with(~ "cell_id", matches("target_cluster"))
umap_df
umap_df$cell_id <- umap_df$name_clean
umap_df_with_ratios <- umap_df %>%
  left_join(surrounding_ratios_wide, by = "cell_id")

umap_df$cell_id <- rownames(umap_df)
head(umap_df_with_ratios)

num_colnames <- names(umap_df_with_ratios)[
  sapply(umap_df_with_ratios, is.numeric) & 
    !(names(umap_df_with_ratios) %in% c("UMAP1", "UMAP2", "cluster", "sample", "disease", "name", "cell_id", "pseudotime"))
]

data_long <- umap_df_with_ratios %>%
  pivot_longer(
    cols = all_of(num_colnames),
    names_to = "variable",
    values_to = "value"
  )

target_clusters <- c("Proximal_Tubules", "Proximal_Damaged", "LoH", "LoH_Damaged")   # 例
target_clusters <- c("Fibroblast_FAP","Fibroblast", "TREM2_Macrophage", "Neutrophil")   # 例
target_clusters <- c("CD4T","CD8T","Bcell")   # 例
target_clusters <- c("Fibroblast", "Fibroblast_FAP", "TREM2_Macrophage", "Resident_Macrophage", "Inflammatory_Macrophage")
target_clusters <- c("Neutrophil")

data_long_filtered <- data_long %>%
  filter(variable %in% target_clusters)

max_pseudotime <- max(data_long_filtered$pseudotime, na.rm = TRUE)
x_start <- 1
arrow_y <- 0
y_end <- 0.5 

cluster_colors <- c(
  "CD4T" = "#1f77b4",
  "CD8T" = "#ff7f0e",
  "CD4Eff" = "#1f77b4",
  "CD8Eff" = "#ff7f0e",
  "LoH" = "#2ca02c",
  "Proliferating Cells" = "#2ca02c",
  "Eosino/Baso" = "#2ca02c",
  "Tcelldim" = "#d62728",
  "Juxtaglomerular" = "#d62728",
  "TREM2_Macrophage" = "#9467bd",
  "Resident_Macrophage" = "#8c564b",
  "Inflammatory_Macrophage" = "#e377c2",
  "Neutrophil" = "#e377c2",
  "Fibroblast_FAP" = "#7f7f7f",
  "Fibroblast" = "#bcbd22",
  "Podocyte" = "#17becf",
  "Mesangium" = "#aec7e8",
  "PEC" = "#ffbb78",
  "EC_GC" = "#98df8a",
  "Glomerular_EC" = "#98df8a",
  "Proximal_Tubules" = "#c5b0d5",
  "Proximal_PCT" = "#c5b0d5",
  "Proximal_Damaged" = "#c49c94",
  "Reg_PT" = "#c49c94",
  "Bcell" = "#c45c32",
  "Proximal_PST" = "#c45c32",
  "LoH_Damaged" = "#17becf"
)

p <- ggplot(data_long_filtered, aes(x = pseudotime, y = value, color = variable, fill = variable)) +
  geom_smooth(method = "loess", se = TRUE, span = 1, alpha = 0.15) +
  theme_classic() +#theme(legend.position = "none")+
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
           arrow = arrow(ends = "last", type = "closed", length = unit(0.1, "inches"))) 
