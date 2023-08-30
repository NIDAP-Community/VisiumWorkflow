#' Differential Expression in Visium object between annotated regions
#'
#' @return
#' @export
#'
#' @examples
visiumDifferentialExpression <- function(object,
                                         region.1,
                                         region.2,
                                         genes,
                                         min.transparency = 0.1,
                                         max.transparency = 1) {

  # Create a list of DE markers between the chosen regions
  de.markers <- FindMarkers(object, ident.1 = region.1, ident.2 = region.2)


  # Create a list of all selected genes that are DEGS
  de.gene.list <- c()

  de.gene.names <- rownames(de.markers)

  for (gene in genes) {
    if (gene %in% de.gene.names) {
      de.gene.list <- append(de.gene.list, gene)
    } else {
      print(paste0(gene, " was not found in the DEG list"))
    }
  }

  # Print the plots of all selected genes that are DEGs
  if (length(de.gene.list) == 0) {
    print("None of the chosen genes were found in the DEG list")
  } else {
    de.spatial.plot <- SpatialFeaturePlot(object = object,
                                          features = de.gene.list,
                                          alpha = c(min.transparency,max.transparency),
                                          ncol = 3)
  }

}
