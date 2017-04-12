#' An S4 class for representing colormaps
#' @slot cmap_type A character string for the type of colormap. Valid values are "corr_values", "tvalues", "log_pvalues", "log_pvalues_adjusted"
#' @slot cmap_name A character string for the title of the colormap. This string will be displayed on the colorbar.
#' @slot values A numeric vector containing the values or measures to be mapped.
#' @slot cex Maximum absolute value of the measure to be mapped.
#' @slot rgbcolors A matrix of RGB colors whose size is same as that of the values.
#' @slot lut Color look up table
#' @slot vmin,vmax Minimum and maximum values
#'
#' @export
BssColormap <- setClass(
  "BssColormap",
  slots = list(
    cmap_type = "character",
    cmap_name = "character",
    cex = "numeric",
    values = "numeric",
    rgbcolors = "matrix",
    lut = "character",
    vmin = "numeric",
    vmax = "numeric"
  )
)

setMethod("initialize", valueClass = "BssColormap", signature = "BssColormap",
          function(.Object, cmap_type, cmap_name, values) {
            .Object@cmap_type <- cmap_type
            .Object@cmap_name <- cmap_name
            .Object@cex <- max(abs(values))
            .Object@values <- values
            .Object@lut <- ""
            .Object@vmin <- 0
            .Object@vmax <- 0
            switch(.Object@cmap_type,
                   corr_values = { cmap <- get_tvalue_colors(cmap_name, values)}, # Use tvalue cmap for correlations
                   tvalues = { cmap <- get_tvalue_colors(cmap_name, values)},
                   log_pvalues = { cmap <- get_logpvalue_colors(cmap_name, values)},
                   log_pvalues_adjusted = { cmap <- get_logpvalue_colors(cmap_name, values)}
            )
            .Object@lut <- cmap$lut
            .Object@rgbcolors <- cmap$rgbcolors
            .Object@vmin <- cmap$vmin
            .Object@vmax <- cmap$vmax
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
         log_pvalues = { bss_cmap@rgbcolors <-
           get_logpvalue_colormap(bss_cmap@cmap_name, bss_cmap@values)
         },
         tvalues = { bss_cmap@rgbcolors <-
           get_tvalue_colors(bss_cmap@cmap_name, bss_cmap@values)
         }
  )

  return(bss_cmap@rgbcolors)

})

get_logpvalue_colors <- function(cmap_name, values) {

  #                   log10(0.05)/1.0001   -log10(0.05)/1.0001
  # |---------------|---|---------|----------|---|-------------------|
  # -pex         log10(0.05)      0             -log10(0.05)        pex

  # **important** Assume values is already log-transformed
  pex <- max(abs(values))
  pFDRneglog <-  1*log10(0.05)
  pFDRposlog <- -1*log10(0.05)
  N <- 256
  if (pex < -1*log10(0.05)) {
    pex <- -1*log10(0.05)*1.001
    lut <- get_color_palette('white', N)
  }
  else {
    totlen <- 2*pex
    neglen <- pex - abs(pFDRneglog)
    zerolen <- pFDRposlog - pFDRneglog
    poslen <- pex - pFDRposlog

    negcolors <- get_color_palette('rev_winter', round(neglen*N/(1.001*totlen)))
    zerocolors <- get_color_palette('white', round(zerolen*N/(1.001*totlen)))
    poscolors <- get_color_palette('spring', round(poslen*N/(1.001*totlen)))
    lut <- c(negcolors, zerocolors, poscolors)
  }
  lut <- colorRampPalette(lut)(256) # Set the length of the lut to 256
  fnmap <- colorRamp(lut)
  values_0_to_1 <- (values + pex ) / (2*pex + .Machine$double.eps)
  rgbcolors <- fnmap(values_0_to_1)/255

  return(list("rgbcolors"=rgbcolors, "lut"=lut, "vmin"=-1*pex, "vmax"=pex))
}

get_tvalue_colors <- function(cmap_name, values) {

  # |---------------|-----------|------------|-------------------|
  # tnegmax     tnegmin         0         tposmin             tposmax

  N <- 256
  tnegmax <- min( c(values[values < 0], 0) )
  tnegmin <- max( c(values[values < 0], tnegmax) )
  tposmax <- max( c(values[values > 0], 0) )
  tposmin <- min( c(values[values > 0], tposmax) )

  totlen <- tposmax - tnegmax
  poslen <- tposmax - tposmin
  zerolen <- tposmin - tnegmin
  neglen <- abs(tnegmax) - abs(tnegmin)

  if ( all(totlen == 0) ) {
    lut <- get_color_palette('white', N)
  }
  else {
    negcolors <- get_color_palette('rev_winter', round(neglen/(1.001*totlen)*N))
    zerocolors <- get_color_palette('white', round(zerolen/(1.001*totlen)*N))
    poscolors <- get_color_palette('spring', round(poslen/(1.001*totlen)*N))
    lut <- c(negcolors, zerocolors, poscolors)
  }
  lut <- colorRampPalette(lut)(256) # Set the length of the lut to 256
  fnmap <- colorRamp(lut)
  values_0_to_1 <- (values - tnegmax ) / (tposmax - tnegmax + .Machine$double.eps)
  rgbcolors <- fnmap(values_0_to_1)/255

  return(list("rgbcolors"=rgbcolors, "lut"=lut, "vmin"=tnegmax, "vmax"=tposmax))
}

colorbar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)

  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

save_colorbar <- function(filename, lut, vmin, vmax, labeltxt) {
  df <- data.frame(
    y=seq(vmin, vmax, length=256)
  )

  ggplot2::ggplot(df) +
    ggplot2::geom_raster(ggplot2::aes(x = 0.5, y=df$y, fill = df$y)) +
    ggplot2::scale_fill_gradientn(colours = lut)   +  ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank()) +
    ggplot2::guides(fill=FALSE) +
    ggplot2::labs(y=labeltxt) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=14)) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 16, vjust=-1)) +
    ggplot2::theme(axis.ticks.length=ggplot2::unit(0.25, "cm"),
          axis.text.y = ggplot2::element_text(margin=ggplot2::unit(c(1.5,1.5,1.5,1.5), "cm")) ) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks= scales::pretty_breaks(n=10), position='right') +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(plot.background = ggplot2::element_blank()) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1)) +
    ggplot2::ggsave(filename, device = "pdf", width = 1.3, height = 3.5, dpi = 600)

}

#' Get a color palette (Hex color code list) for a colormap
#' @param cmap_name name of the colormap
#' @param N number of colors
#' @export
get_color_palette <- function(cmap_name, N) {

  switch(cmap_name,
         rev_spring = {
           return ( rev(colorRampPalette(c("#FF00FFFF", "#FFFF00FF"))(ceiling(N))) )
         },
         spring = {
           return ( colorRampPalette(c("#FF00FFFF", "#FFFF00FF"))(ceiling(N)) )
         },
         rev_winter = {
           return ( rev(colorRampPalette(c("#0000FF", "#00FF80"))(ceiling(N))) )
         },
         winter = {
           return ( colorRampPalette(c("#0000FF", "#00FF80"))(ceiling(N)) )
         },
         white = {
           return ( colorRampPalette(c("#FFFFFF"))(ceiling(N)) )
         }
         )
}

save_BrainSuiteLUT <- function(filename, lut) {
  write(col2rgb(lut)/255, file=filename, sep = " ", ncolumns = 3)
}
