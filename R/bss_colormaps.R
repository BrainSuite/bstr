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

get_logpvalue_colormap <- function(cmap_name, values) {
  hexcolrs <- RColorBrewer::brewer.pal(11, cmap_name)
  hexcolrs[6] <- "#FFFFFF"
  whitecolr <-  "#FFFFFF"
  colfunc <- colorRampPalette(c("white"))
  colfunc(10)
  # First 5 colors are negative
  # Last 5 colors are positive
  negcolors <- hexcolrs[1:5]
  poscolors <- hexcolrs[7:11]

  fnmap <- colorRamp(hexcolrs)

  values_0_to_1 <- ( values - min(values) ) /
    (max(values) - min(values))

  rgbcolors <- fnmap(values_0_to_1)/255
  return (rgbcolors)

}

setMethod("get_colors", valueClass = "matrix", signature = "BssColormap", function(bss_cmap) {

  switch(bss_cmap@cmap_type,
         logpvalue = { bss_cmap@rgbcolors <-
           get_logpvalue_colormap(bss_cmap@cmap_name, bss_cmap@values)
         },
         tvalue = { bss_cmap@rgbcolors <-
           get_tvalue_colormap(bss_cmap@cmap_name, bss_cmap@values)
         }
  )

  return(bss_cmap@rgbcolors)

})
