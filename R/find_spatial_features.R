#' Find spatially variable features from non-annotated groups
#'
#' @param object
#'
#' @importFrom Seurat VariableFeatures
#' @importFrom Seurat FindSpatiallyVariableFeatures
#' @importFrom Seurat SpatiallyVariableFeatures
#' @importFrom Seurat SpatialFeaturePlot
#'
#' @return
#' @export
#'
#' @examples
findSpatialFeatures <- function(object,
                                min.transparency = 0.1,
                                max.transparency = 1){

  # A list of variable features to evaluate
  var.features <- VariableFeatures(object)[1:1000]

  # Assess variable features in spatial context
  object <- FindSpatiallyVariableFeatures(object,
                                          assay = "SCT",
                                          features = var.features,
                                          selection.method = "moransi")

  # List all spatially variable features
  all.top.features <- SpatiallyVariableFeatures(object,
                                                selection.method = "moransi")

  # Choose the top 6 spatially variable features and plot
  select.top.features <- head(all.top.features, 6)
  select.feature.plot <- SpatialFeaturePlot(object,
                                            features = top.features,
                                            ncol = 3,
                                            alpha = c(min.transparency,
                                                      max.transparency))

  return(list("object"= object, "select.feature.plot" = select.feature.plot))

}
