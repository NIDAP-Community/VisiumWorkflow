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
visiumDimReduction <- function(object,
                               clusters.to.plot = c(-1)) {

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

  # Define specific clusters for spatial plot

  # Default to all clusters if no specific cluster was chosen
  if (-1 %in% clusters.to.plot) {
    cell.identity <- CellsByIdentities(object = object)
  } else {
    cell.identity <- CellsByIdentities(object = object,
                                       idents = clusters.to.plot)
  }

  # Spatial plot for chosen clusters
  identify.spatial.plot <- SpatialDimPlot(object,
                                             cells.highlight = cell.identity,
                                             facet.highlight = TRUE,
                                             ncol = 3)

  return(list("object" = object,
              "umap.plot" = umap.plot,
              "spatial.cluster.plot" = spatial.cluster.plot,
              "identify.spatial.plot" = identify.spatial.plot))

  }
