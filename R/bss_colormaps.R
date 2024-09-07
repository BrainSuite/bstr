# BrainSuite Statistics Toolbox in R (bstr)
# Copyright (C) 2017 The Regents of the University of California
# Creator: Shantanu H. Joshi, Department of Neurology, Ahmanson Lovelace Brain Mapping Center, UCLA
#
# This program is free software; you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation; version 2.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License version 2 for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

#' An S4 class for representing colormaps
#' @slot cmap_type A character string for the type of colormap. Valid values are "corr_values", "corr_values_masked_adjusted", "tvalues", "tvalues_adjusted", "log_pvalues", "log_pvalues_adjusted"
#' @slot cmap_name A character string for the title of the colormap. This string will be displayed on the colorbar.
#' @slot values A numeric vector containing the values or measures to be mapped.
#' @slot cex Maximum absolute value of the measure to be mapped.
#' @slot rgbcolors A matrix of RGB colors whose size is same as that of the values.
#' @slot lut Color look up table
#' @slot vmin,vmax Minimum and maximum values
#' @slot cnegmin Minimum negative value
#' @slot cnegmax Maximum negative value
#' @slot cposmin Minimum positive value
#' @slot cposmax Maximum positive value
#'
#' @export
BstrColormap <- setClass(
  "BstrColormap",
  slots = list(
    cmap_type = "character",
    cmap_name = "character",
    cex = "numeric",
    values = "numeric",
    rgbcolors = "matrix",
    lut = "character",
    vmin = "numeric",
    vmax = "numeric",
    cnegmin = "numeric",
    cnegmax = "numeric",
    cposmin = "numeric",
    cposmax = "numeric"
  )
)

setMethod("initialize", valueClass = "BstrColormap", signature = "BstrColormap",
          function(.Object, cmap_type, cmap_name, values) {
            .Object@cmap_type <- cmap_type
            .Object@cmap_name <- cmap_name
            .Object@cex <- max(abs(values))
            .Object@values <- values
            .Object@lut <- ""
            .Object@vmin <- 0
            .Object@vmax <- 0
            .Object@cnegmin <- 0
            .Object@cnegmax <- 0
            .Object@cposmin <- 0
            .Object@cposmax <- 0

            switch(.Object@cmap_type,
                   corr_values = { cmap <- get_tvalue_colors(cmap_name, values)}, # Use tvalue cmap for correlations
                   corr_values_masked_adjusted = { cmap <- get_tvalue_colors(cmap_name, values)},
                   tvalues = { cmap <- get_tvalue_colors(cmap_name, values)},
                   tvalues_adjusted = { cmap <- get_tvalue_colors(cmap_name, values)},
                   log_pvalues = { cmap <- get_logpvalue_colors(cmap_name, values)},
                   log_pvalues_adjusted = { cmap <- get_logpvalue_colors(cmap_name, values)}
            )
            .Object@lut <- cmap$lut
            .Object@rgbcolors <- cmap$rgbcolors
            .Object@vmin <- cmap$vmin
            .Object@vmax <- cmap$vmax
            .Object@cnegmin <- cmap$cnegmin
            .Object@cnegmax <- cmap$cnegmax
            .Object@cposmin <- cmap$cposmin
            .Object@cposmax <- cmap$cposmax

            return(.Object)
          })

setGeneric("get_colors", valueClass = "matrix",function(bstr_cmap) {
  standardGeneric("get_colors")
})
#' Get log pvalue colormap
#' @param cmap_name name of the colormap
#' @param values log-transformed p-values
#' @export
get_logpvalue_colormap <- function(cmap_name, values) {
  hexcolrs <- RColorBrewer::brewer.pal(11, cmap_name)
  hexcolrs[6] <- hexcolrstr$gray
  colfunc <- colorRampPalette(c("gray"))
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

setMethod("get_colors", valueClass = "matrix", signature = "BstrColormap", function(bstr_cmap) {

  switch(bstr_cmap@cmap_type,
         log_pvalues = { bstr_cmap@rgbcolors <-
           get_logpvalue_colormap(bstr_cmap@cmap_name, bstr_cmap@values)
         },
         tvalues = { bstr_cmap@rgbcolors <-
           get_tvalue_colors(bstr_cmap@cmap_name, bstr_cmap@values)
         }
  )

  return(bstr_cmap@rgbcolors)

})
#' Get log pvalue colors, generate lut, and calculate min/max values
#' @param cmap_name name of the colormap
#' @param values log-transformed p-values
#' @export
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
    lut <- get_color_palette('gray', N)
  }
  else {
    totlen <- 2*pex
    neglen <- pex - abs(pFDRneglog)
    zerolen <- pFDRposlog - pFDRneglog
    poslen <- pex - pFDRposlog

    negcolors <- get_color_palette('rev_winter', round(neglen*N/(1.001*totlen)))
    zerocolors <- get_color_palette('gray', round(zerolen*N/(1.001*totlen)))
    poscolors <- get_color_palette('spring', round(poslen*N/(1.001*totlen)))
    lut <- c(negcolors, zerocolors, poscolors)
  }
  lut <- colorRampPalette(lut)(256) # Set the length of the lut to 256
  fnmap <- colorRamp(lut)
  values_0_to_1 <- (values + pex ) / (2*pex + .Machine$double.eps)
  rgbcolors <- fnmap(values_0_to_1)/255

  cnegmin <- pFDRneglog
  cnegmax <- -pex
  cposmin <- pFDRposlog
  cposmax <- pex

  return(list("rgbcolors"=rgbcolors, "lut"=lut, "vmin"=-1*pex, "vmax"=pex,
              "cnegmin"=cnegmin, "cnegmax"=cnegmax, "cposmin"=cposmin, "cposmax"=cposmax))
}
#' Get tvalue colors, generate lut, and calculate min/max values
#' @param cmap_name name of the colormap
#' @param values t-values
#' @export
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
    lut <- get_color_palette('gray', N)
  }
  else {
    negcolors <- get_color_palette('rev_winter', round(neglen*N/(1.001*totlen)))
    zerocolors <- get_color_palette('gray', round(zerolen*N/(1.001*totlen)))
    poscolors <- get_color_palette('spring', round(poslen*N/(1.001*totlen)))
    lut <- c(negcolors, zerocolors, poscolors)
  }
  lut <- colorRampPalette(lut)(256) # Set the length of the lut to 256
  fnmap <- colorRamp(lut)
  values_0_to_1 <- (values - tnegmax ) / (tposmax - tnegmax + .Machine$double.eps)
  rgbcolors <- fnmap(values_0_to_1)/255


  return(list("rgbcolors"=rgbcolors, "lut"=lut, "vmin"=tnegmax, "vmax"=tposmax,
              "cnegmin"=tnegmin, "cnegmax"=tnegmax, "cposmin"=tposmin, "cposmax"=tposmax))
}
#' Generates the components of a colorbar for an individual voxelcoordinate and overlay
#' @param lut previously generated lut file
#' @param min minimum value for the colorbar
#' @param max maximum value for the colorbar (default is the negative of the minimum)
#' @param ticks vector of ticks
#' @param nticks number of desired tick marks
#' @param title title of the colorbar
#' @export
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
#' Creates the display of a colorbar for an individual voxelcoordinate and overlay
#' @param filename name of the file
#' @param lut previously generated lut file
#' @param vmin minimum value for the colorbar
#' @param vmax maximum value for the colorbar
#' @param labeltxt the label for the colorbar
#' @export
save_colorbar <- function(filename, lut, vmin, vmax, labeltxt) {
  df <- data.frame(
    y=seq(vmin, vmax, length=256)
  )

  ggplot2::ggplot(df) +
    ggplot2::geom_raster(ggplot2::aes_string(x = 0.5, y='y', fill = 'y')) +
    ggplot2::scale_fill_gradientn(colours = lut)   +  ggplot2::xlab('') + ggplot2::ylab('') +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank()) +
    ggplot2::guides(fill="none") +
    ggplot2::labs(y=labeltxt) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=16)) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(size = 18, vjust=0.0)) +
    ggplot2::theme(axis.ticks.length=ggplot2::unit(0.25, "cm"),
          axis.text.y = ggplot2::element_text(margin=ggplot2::unit(c(2,2,2,2), "cm")) ) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks= scales::pretty_breaks(n=10), position='right') +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::theme(plot.background = ggplot2::element_blank()) +
    ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1))
    ggplot2::ggsave(filename, width = 1.5, height = 4, dpi = 600)

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
           return ( colorRampPalette(c(hexcolrstr$white))(ceiling(N)) )
         },
         gray = {
           return ( colorRampPalette(c(hexcolrstr$gray))(ceiling(N)) )
         }
        )
}
#' Saves the colormap values for min and max to an ini file
#' @param filename name of the file
#' @param bstr_cmap BstrColormap object
#' @export
save_colormap_to_ini <- function(filename, bstr_cmap) {
  cmap_to_save <- list()
  cmap_to_save[["colormap"]] <- list(cmap_type=bstr_cmap@cmap_type,
                                     cmap_name = bstr_cmap@cmap_name,
                                     cnegmin = bstr_cmap@cnegmin,
                                     cnegmax = bstr_cmap@cnegmax,
                                     cposmin = bstr_cmap@cposmin,
                                     cposmax = bstr_cmap@cposmax)
  ini::write.ini(cmap_to_save, filename)
}
#' Saves and writes the lut for BrainSuite use
#' @param filename name of the file
#' @param lut previously generated lut
#' @export
save_BrainSuiteLUT <- function(filename, lut) {
  write(col2rgb(lut)/255, file=filename, sep = " ", ncolumns = 3)
}

hexcolrstr <- list(
  white = "#FFFFFF",
  gray = "#CCCCCC"
)
