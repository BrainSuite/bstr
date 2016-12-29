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

get_tvalue_colormap <- function(cmap_name, values) {
  hexcolrs <- rev(RColorBrewer::brewer.pal(11, cmap_name))
  hexcolrs[6] <- "#FFFFFF"
  N <- 256
  if ( all(values >= 0) ){
    if (min(values) == 0)
      tposmin <- min(values[values > 0])
    else
      tposmin <- min(values)
    tposmax <- max(values)
    tnegmax <- 0
    fnmap <-
      colorRamp(colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'YlOrRd')))(1:256))
  }
  else {
    tnegmin <- -1*min(abs(values[values < 0]))
    tnegmax <- -1*max(abs(values[values < 0]))
    tposmin <- min(abs(values[values > 0]))
    tposmax <- max(abs(values[values > 0]))
    totlen <- tposmax - tnegmax
    if ( tnegmin < totlen/256) {
      tnegmin <- sign(tnegmin)*totlen/256
    }

    if ( tposmin < totlen/256) {
      tposmin <- sign(tposmin)*totlen/256
    }

    neglen <- tnegmin - tnegmax
    poslen <- tposmax - tposmin
    midzero_len <- tposmin - tnegmin


    negcolors <- rev(hexcolrs[1:5]) # First 5 colors are negative
    poscolors <- rev(hexcolrs[7:11]) # Last 5 colors are positive

    negcolor_range <- colorRampPalette(negcolors)(round(neglen/(1.001*totlen)*256))
    midzero_color_range <- colorRampPalette(hexcolrs[6])(round(midzero_len/totlen*256))
    poscolor_range <- colorRampPalette(poscolors)(round(poslen/(1.001*totlen)*256))

    lut <- c(negcolor_range, midzero_color_range, poscolor_range)
    fnmap <- colorRamp(lut)
  }

  values_0_to_1 <- (values - tnegmax ) / (tposmax - tnegmax)
  rgbcolors <- fnmap(values_0_to_1)/255
  return (rgbcolors)

}
