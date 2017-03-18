#' LabelDesc I/O functions
#' @export
read_label_desc <- function() {
  label_desc_filename <- system.file("extdata", "brainsuite_labeldescriptions_14May2014.xml", package = 'bssr')
  x1 <- xml2::read_xml(label_desc_filename)
  labels <-xml2:::xml_find_all(x1, './/label')
  return (data.frame(roiid = xml2::xml_attr(labels, 'id'), roiname = xml2::xml_attr(labels, 'fullname')))
}

#' @export
get_roi_name <- function(label_desc_df, roiid) {
  return(label_desc_df[label_desc_df$roiid == roiid,]['roiname'])
}
