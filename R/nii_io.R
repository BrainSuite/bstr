read_nii_images_for_all_subjects <- function(nii_filelist, attrib_siz, mask_idx = NULL) {
  data_matrix <- vapply(nii_filelist, function(i) {
    cat(sprintf('%s\n', i))
    as.vector(RNifti::readNifti(i))
  }, FUN.VALUE = numeric(attrib_siz))
  colnames(data_matrix) <- NULL
  if ( !is.null(mask_idx) )
    data_matrix <- data_matrix[mask_idx,]
  return(t(data_matrix))
}
