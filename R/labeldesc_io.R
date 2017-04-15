#' #' Read BrainSuite label description file
#' #' @export
#' read_label_desc_xml2 <- function() {
#'   label_desc_filename <- system.file("extdata", "brainsuite_labeldescriptions_14May2014.xml", package = 'bssr')
#'   x1 <- xml2::read_xml(label_desc_filename)
#'   labels <-xml2::xml_find_all(x1, './/label')
#'   return (data.frame(roiid = xml2::xml_attr(labels, 'id'), roiname = xml2::xml_attr(labels, 'fullname')))
#' }

#' Read BrainSuite label description file
#' @export
read_label_desc <- function() {
  label_desc_filename <- system.file("extdata", "brainsuite_labeldescriptions_14May2014.xml", package = 'bssr')
  fid <- file(label_desc_filename, open="r")
  alllines <- readLines(fid)
  alllines <- alllines[3:(length(alllines)-1)]
  roiid <- vector(mode = "numeric", length = length(alllines))
  roiname <- vector(mode = "character", length = length(alllines))
  for (ii in seq(alllines)) {
    tempstr <- gsub("<|/>", "", gsub("\"", "", alllines[ii]))
    roiid[ii] <- unlist(strsplit(unlist(strsplit(tempstr, " "))[2], '='))[2]
    roiname[ii] <- unlist(strsplit(tempstr, '='))[5]
  }
  close(fid)
  return (data.frame(roiid = roiid, roiname = roiname))
}


#' Get the ROI name from the ROI label
#' @param label_desc_df \code{\link{data.frame}} containing fields from the label description file
#' @param roiid ROI label identifier
#' @export
get_roi_name <- function(label_desc_df, roiid) {
  return(label_desc_df[label_desc_df$roiid == roiid,]['roiname'])
}
