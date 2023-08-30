#' Create the initial Seurat ibject from the Visium inputs
#'
#' @importFrom Seurat Load10X_Spatial
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

  Load10X_Spatial(data.dir = visium.files.directory,
                  slice = tissue.slice.name,
                  filename = feature.file.name
                  )

  return()

}
