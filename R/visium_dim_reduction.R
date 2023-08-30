#' Perform dimensionality reduction on Visium Object
#'
#'
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat RunUMAP
#' @importFrom Seurat DimPlot
#' @importFrom Seurat SpatialDimPlot
#'
#' @return
#' @export
#'
#' @examples
visiumDimReduction <- function(object){

  # Perform dimensionality reduction
  object <- RunPCA(object, assay = "SCT", verbose = FALSE)
  object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
  object <- FindClusters(object, verbose = FALSE)

  # Run UMAP
  object <- RunUMAP(object, reduction = "pca", dims= 1:30)

  # Create the UMAP plot with labeled clusters
  umap.plot <- DimPlot(object, reduction = "umap", label = TRUE)

  # Spatially plot the UMAP clusters on the tissue image
  spatial.cluster.plot <- SpatialDimPlot(object, lable = TRUE, label.size = 3)

  return(list("object" = object,
              "umap.plot" = umap.plot,
              "spatial.cluster.plot" <- spatial.cluster.plot))

  }
