#' Bss output class
#' Defines an S4 class for statistical outputs
#' @export

check_files <- function(object){
  if (!dir.exists(object@subjdir)) {
    stop(sprintf("Subjects directory %s does not exist.\n", object@subjdir), call. = FALSE)
  }

  if (!file.exists(object@csv)) {
    stop(sprintf("Demographics csv file %s does not exist.\n", object@csv), call. = FALSE)
  }
}

BssOutput <- setClass(
  "BssOutput",
  slots = list(
    outdir = "character"
  ),
  validity = check_files
)

setMethod("initialize", valueClass = "BssOutput", signature = "BssOutput", function(.Object, outdir) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
    .Object@outdir <- outdir
  }
  else {
    .Object@outdir <- outdir
    warning(sprintf("The output directory %s already exists. Will overwrite it's contents.", outdir),
            immediate. = TRUE, call.=FALSE)
  }
  return(.Object)
})

setGeneric("save_out", valueClass = "BssOutput", function(bss_out, bss_data, bss_model) {
  standardGeneric("save_out")
})

setMethod("save_out", valueClass = "BssOutput", signature = "BssOutput", function(bss_out, bss_data, bss_model) {
  print("in save")
})

BssCBMOutput <- setClass(
  "BssCBMOutput",
  contains = "BssOutput"
)

BssTBMOutput <- setClass(
  "BssTBMOutput",
  contains = "BssOutput"
)

BssROIOutput <- setClass(
  "BssROIOutput",
  contains = "BssOutput"
)

setMethod("save_out", valueClass = "BssCBMOutput", signature = "BssCBMOutput", function(bss_out, bss_data, bss_model) {

  log_pvalues <- log10_transform(bss_model@pvalues)
  s1 <- bss_data@atlas_surface

  s1$attributes <- log_pvalues
  bss_cmap <- new("BssColormap", "logpvalue", "Spectral", log_pvalues)
  s1$vColor <- get_colors(bss_cmap)
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  writedfs(file.path(bss_out@outdir, 'log_pvalues.dfs'), s1)

  log_pvalues <- log10_transform(p.adjust(bss_model@pvalues, method = 'BH'))
  s1$attributes <- log_pvalues
  bss_cmap <- new("BssColormap", "logpvalue", "Spectral", log_pvalues)
  s1$vColor <- get_colors(bss_cmap)
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  writedfs(file.path(bss_out@outdir, 'log_pvalues_adjusted.dfs'), s1)


  bss_model@tvalues[abs(log_pvalues) <= -1*log10(0.05)] <- 0
  bss_cmap <- new("BssColormap", "tvalue", "Spectral", bss_model@tvalues)
  s1$attributes <- bss_model@tvalues
  s1$vColor <- get_colors(bss_cmap)
  s1$vColor <- matrix(s1$vColor, nrow=3, ncol=s1$hdr$nVertices, byrow = TRUE)
  writedfs(file.path(bss_out@outdir, 'tvalues.dfs'), s1)

  return(bss_out)
  }
)

