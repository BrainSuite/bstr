#' Defines an S4 class for colormaps
#' @export

BssColormap <- setClass(
  "BssColormap",
  slots = list(
    cmap_type = "character",
    cmap_name = "character",
    cex = "numeric",
    values = "numeric",
    rgbcolors = "matrix"
  )
)

setMethod("initialize", valueClass = "BssColormap", signature = "BssColormap",
          function(.Object, cmap_type, cmap_name, values) {
            .Object@cmap_type <- cmap_type
            .Object@cmap_name <- cmap_name
            .Object@cex <- max(abs(values))
            .Object@values <- values
            return(.Object)
          })

setGeneric("get_colors", valueClass = "matrix",function(bss_cmap) {
  standardGeneric("get_colors")
})
