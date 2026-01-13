LoadXenium2 <- function(data.dir, fov = 'fov', assay = 'Xenium') {
  data <- ReadXenium(
    data.dir = data.dir,
    type = c("centroids", "segmentations"),
  )
  
  segmentations.data <- list(
    "centroids" = CreateCentroids(data$centroids),
    "segmentation" = CreateSegmentation(data$segmentations)
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("segmentation", "centroids"),
    molecules = data$microns,
    assay = assay
  )
  
  xenium.obj <- CreateSeuratObject(counts = data$matrix[["Gene Expression"]], assay = assay)
  xenium.obj[["ControlCodeword"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Codeword"]])
  xenium.obj[["ControlProbe"]] <- CreateAssayObject(counts = data$matrix[["Negative Control Probe"]])
  FeaturePlot(xenium.obj, "PROM1")
  xenium.obj[[fov]] <- coords
  return(xenium.obj)
}

xenium.obj.list <- list()
for (i in 1:length(folders)){
xenium_data <- LoadXenium2(folders[[i]], fov="fov")
xenium_data$sample <- paste0("sample",i)
xenium.obj.list[[i]] <- xenium_data
}
rm(xenium_data)
options(future.globals.maxSize = 8 * 1024^3) 
xenium.obj <- merge(xenium.obj.list[[1]],c(xenium.obj.list[[2]],xenium.obj.list[[3]],xenium.obj.list[[4]],
                                           xenium.obj.list[[5]],xenium.obj.list[[6]],xenium.obj.list[[7]],xenium.obj.list[[8]],
                                           xenium.obj.list[[9]],xenium.obj.list[[10]],xenium.obj.list[[11]],xenium.obj.list[[12]]))
rm(xenium.obj.list)
gc()
xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 10 & nFeature_Xenium > 0)

#xenium.obj <- SCTransform(xenium.obj, assay = "Xenium", conserve.memory = T, ncells=1000, variable.features.n = 750)
##############
xenium.obj <- FindVariableFeatures(xenium.obj)
xenium.obj <- NormalizeData(xenium.obj)
xenium.obj <- ScaleData(xenium.obj)
gc()
##############

xenium.obj <- RunPCA(xenium.obj, npcs = 20)
gc()
xenium.obj <- IntegrateLayers(xenium.obj, method=HarmonyIntegration,
                              orig.reduction="pca",new.reduction="integrated.cca",verbose=T)

xenium.obj <- RunUMAP(xenium.obj, reduction="integrated.cca", dims=1:15, n.neighbors=15)
gc()
xenium.obj <- FindNeighbors(xenium.obj, reduction = "integrated.cca", dims = 1:20, k.param=15)
gc()
xenium.obj <- FindClusters(xenium.obj, resolution = 3, n.start=7, n.iter=7)


xenium.obj <- RenameIdents(xenium.obj, "")
Idents(xenium.obj) <- "Xenium_snn_res.1"
DimPlot(xenium.obj, label=T, repel=T, group.by="clusters", split.by="disease")
