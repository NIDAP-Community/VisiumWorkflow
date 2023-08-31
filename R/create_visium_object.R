#' Create the initial Seurat ibject from the Visium inputs
#'
#' @importFrom Seurat Load10X_Spatial
#' @importFrom Seurat VlnPlot
#' @importFrom Seurat SpatialFeaturePlot
#' @importFrom ggplot theme
#'
#' @return
#' @export
#'
#' @examples
createVisiumObject <- function(visium.files.directory,
                               tissue.slice.name = "slice1",
                               feature.file.name
                               ) {

  # Load the Visium experiment data as a Seurat object

  object <- Load10X_Spatial(data.dir = visium.files.directory,
                  slice = tissue.slice.name,
                  filename = feature.file.name
                  )

  # Show the distribution of counts in all spots using a violin plot
  violin.plot <- VlnPlot(object, features = "nCount_Spatial", pt.size = 0.1) +
    NoLegend()

  # Show the spatial distribution of counts in each spot on the tissue image
  spatial.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial") +
                    theme(legend.position = "right")

  return(list("object" = object,
              "violin.plot" = violin.plot,
              "spatial.plot" = spatial.plot))

}
