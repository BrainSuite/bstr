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
