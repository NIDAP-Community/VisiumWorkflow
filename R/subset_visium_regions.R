#' Subset the Visium object by selected clusters to isolate a region
#'
#' @param object
#' @param subset.clusters
#'
#' @importFrom Seurat SpatialDimPlot
#'
#' @return
#' @export
#'
#' @examples
subsetVisiumRegions <- function(object,
                                subset.clusters) {


  region <- subset(object, idents = subset.clusters)

  subset.spatial.plot <- SpatialDimPlot(object,
                                        crop = FALSE,
                                        label = TRUE,
                                        pt.size.factor = 1,
                                        label.size = 3)

  return("object" = object, "subset.spatial.plot" = subset.spatial.plot)


}
