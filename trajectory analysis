#trajectory analysis using monocle3
DefaultAssay(s.int) <- "Xenium"
data <- GetAssayData(object = s.int, slot = "counts")
celldata <- as.data.frame(s.int@meta.data)
cds <- new_cell_data_set(data, cell_metadata = celldata)
#cds <- as.cell_data_set(s.int)
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds, reduction_method = "UMAP",umap.n_neighbors=15L, umap.min_dist=0.2,
                        max_components=2)
plot_cells(cds=cds, color_cells_by = "seurat_clusters")
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(s.int[["Xenium"]])
plot_cells(cds, genes="PROM1")

cds <- cluster_cells(cds, reduction_method = "UMAP", cluster_method = "louvain")
cds <- learn_graph(cds, use_partition=F, close_loop=F, verbose=T, learn_graph_control=list(ncenter=200))
cds <- order_cells(cds)

colData(cds)$cluster_custom <- as.character(cds$clusters)


pp <- plot_cells(cds=cds, 
           color_cells_by="pseudotime",
           label_cell_group=F,
           label_groups_by_cluster=F,
           show_trajectory_graph=T,
           label_leaves=F,
           label_branch_points=F,
           label_roots =F,
           cell_size=0.5,
           graph_label_size = 4)
pp

#pseudotime-correlated population analysis
pseudotime_df <- data.frame(
  ind = names(s.sub$Cluster21),
  pseudotime = as.numeric(s.sub$Cluster21)
)

plot_data <- nn_count %>%
  left_join(pseudotime_df, by = "ind")

cell_types <- c("Neutrophil", "CD4T", "CD8T", "Bcell")  # select cells
plot_data_long <- plot_data %>%
  pivot_longer(
    cols = all_of(cell_types),
    names_to = "cell_type",
    values_to = "neighbor_count"
  )
plot_data_long <- plot_data_long %>%
  filter(pseudotime >= -1, pseudotime <= 40)
ggplot(plot_data_long, aes(x = pseudotime, y = neighbor_count, color = cell_type)) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(x = "Pseudotime", y = "Number of Neighboring Cells") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed")

#pseudotime-correlated expression analysis
expr_matrix <- as.matrix(GetAssayData(s.sub, slot = "data")[, non_na_cells])

min_cells <- 0.05 * ncol(expr_matrix)
expressed_genes <- rowSums(expr_matrix > 0) >= min_cells
expr_matrix <- expr_matrix[expressed_genes, ]
cors <- apply(expr_matrix, 1, function(x) cor(x, pseudotime_values, method="spearman"))
pvals <- apply(expr_matrix, 1, function(x) cor.test(x, pseudotime_values, method="spearman")$p.value)

gene_cors_df <- data.frame(
  gene = names(cors),
  correlation = cors, 
  pvalue = pvals,
  stringsAsFactors = FALSE
)

gene_cors_df$padj <- p.adjust(gene_cors_df$pvalue, method = "BH")

up_genes <- gene_cors_df %>%
  filter(correlation > 0 & padj < 0.05) %>%
  arrange(desc(correlation))
down_genes <- gene_cors_df %>%
  filter(correlation < 0 & padj < 0.05) %>%
  arrange(desc(correlation))

Y <- c(up_genes$gene[1:20],down_genes$gene[1:10]) 
Y <- c(up_genes$gene[1:20],down_genes$gene[1:10])  
Y <- Y[Y != c("TNC","CDH6")]

pt.matrix <- exprs(cds)[match(Y,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
#rownames(pt.matrix) <- genes;
#K means with 6 groups
library(ComplexHeatmap)
library(colorRamp2)
library(RColorBrewer)
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 1,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)
htkm
