#' Apply SCTransform normalization to Seurat object counts
#'
#' @importFrom Seurat SCTransform
#'
#'
#' @return
#' @export
#'
#' @examples
normalizeVisiumObject <- function(object,
                                  genes,
                                  point.size = 1.3,
                                  min.transparency = 0.1,
                                  max.transparency = 1
                                  ) {
  object <- SCTransform(object, assay = 'Spatial', verbose = FALSE)

  gene.plot <- SpatialFeaturePlot(object,
                                  features = genes,
                                  pt.size.factor = point.size,
                                  alpha = c(min.transparency,max.transparency))

  return(list("object" = object, "gene.plot" = gene.plot))
}
