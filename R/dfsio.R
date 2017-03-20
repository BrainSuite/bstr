#' I/O for dfs files
#' @export

readdfs <- function(filename) {

  if (! file.exists(filename))
    stop(sprintf('Surface file %s does not exist', filename), call. = FALSE)

  dfs_fid = file(filename, "rb")
  dfs_magic <- readChar(dfs_fid, 12)
  hdr <- data.frame(dfs_magic)
  if (!grepl("DFS", dfs_magic)) {
    stop(sprintf("File %s is not a dfs file.\n", filename), call. = FALSE)
  }
  hdr$hdrsize <- readBin(dfs_fid, integer(), 1)
  hdr$mdoffset <- readBin(dfs_fid, integer(), 1)
  hdr$pdoffset <- readBin(dfs_fid, integer(), 1)
  hdr$nTriangles <- readBin(dfs_fid, integer(), 1)
  hdr$nVertices <- readBin(dfs_fid, integer(), 1)
  hdr$nStrips <- readBin(dfs_fid, integer(), 1)
  hdr$stripSize <- readBin(dfs_fid, integer(), 1)
  hdr$normals = readBin(dfs_fid, integer(), 1)
  hdr$uvStart = readBin(dfs_fid, integer(), 1)
  hdr$vcoffset = readBin(dfs_fid, integer(), 1)
  hdr$labelOffset = readBin(dfs_fid, integer(), 1)
  hdr$vertexAttributes = readBin(dfs_fid, integer(), 1)
  seek(dfs_fid, hdr$hdrsize)
  dfs <- list()
  dfs$hdr <- hdr
  dfs$faces <- readBin(dfs_fid, what = integer(), n = 3 * hdr$nTriangles)
  dfs$faces <- matrix(dfs$faces, nrow=3, ncol=hdr$nTriangles, byrow = FALSE)
  dfs$vertices <- readBin(dfs_fid, double(), size = 4, n = 3 * hdr$nVertices)
  dfs$vertices <- matrix(dfs$vertices, nrow=3, ncol=hdr$nVertices, byrow = FALSE)
  if (hdr$normals > 0) {
    seek(dfs_fid, hdr$normals)
    dfs$normals <- readBin(dfs_fid, double(), size = 4, n = 3 * hdr$nVertices)
    dfs$normals <- matrix(dfs$normals, nrow=3, ncol=hdr$nVertices, byrow = FALSE)
  }
  if (hdr$vcoffset > 0) {
    seek(dfs_fid, hdr$vcoffset)
    dfs$vColor <- readBin(dfs_fid, double(), size = 4, n = 3 * hdr$nVertices)
    dfs$vColor <- matrix(dfs$vColor, nrow=3, ncol=hdr$nVertices, byrow = FALSE)
  }

  if (hdr$uvStart > 0) {
    seek(dfs_fid, hdr$uvStart)
    uv <- readBin(dfs_fid, double(), size = 4, n = 2 * hdr$nVertices)
    uv <- matrix(uv, nrow=2, ncol=hdr$nVertices, byrow = FALSE)
    dfs$u <- uv[1, ]
    dfs$v <- uv[2, ]
  }

  if (hdr$labelOffset > 0) {
    seek(dfs_fid, hdr$labelOffset)
    dfs$labels <- readBin(dfs_fid, integer(), size = 2, n = hdr$nVertices)
  }

  if (hdr$vertexAttributes > 0) {
    seek(dfs_fid, hdr$vertexAttributes)
    dfs$attributes <- readBin(dfs_fid, double(), size = 4, n = hdr$nVertices)
  }

  close(dfs_fid)
  return(dfs)
}

#' @export
readdfsattributes <- function(filename) {
  dfs_fid = file(filename, "rb")
  dfs_magic <- readBin(dfs_fid, character(12))
  hdr <- data.frame(dfs_magic)
  if (!grepl("DFS", dfs_magic)) {
    stop(sprintf("File %s is not a dfs file.\n", filename), call. = FALSE)
  }
  hdr$hdrsize <- readBin(dfs_fid, integer(), 1)
  hdr$mdoffset <- readBin(dfs_fid, integer(), 1)
  hdr$pdoffset <- readBin(dfs_fid, integer(), 1)
  hdr$nTriangles <- readBin(dfs_fid, integer(), 1)
  hdr$nVertices <- readBin(dfs_fid, integer(), 1)
  hdr$nStrips <- readBin(dfs_fid, integer(), 1)
  hdr$stripSize <- readBin(dfs_fid, integer(), 1)
  hdr$normals = readBin(dfs_fid, integer(), 1)
  hdr$uvStart = readBin(dfs_fid, integer(), 1)
  hdr$vcoffset = readBin(dfs_fid, integer(), 1)
  hdr$labelOffset = readBin(dfs_fid, integer(), 1)
  hdr$vertexAttributes = readBin(dfs_fid, integer(), 1)
  seek(dfs_fid, hdr$hdrsize)
  dfs <- list()
  if (hdr$vertexAttributes > 0) {
    seek(dfs_fid, hdr$vertexAttributes)
    dfs$attributes <- readBin(dfs_fid, double(), size = 4, n = hdr$nVertices)
  }

  close(dfs_fid)
  return(dfs$attributes)
}

read_dfs_attributes_for_all_subjects <- function(dfs_filelist, attrib_siz) {
  data_matrix <- vapply(dfs_filelist, function(i) {
    cat(sprintf('%s\n', i))
    readdfsattributes(i)
  }, FUN.VALUE = numeric(attrib_siz))
  colnames(data_matrix) <- NULL
  return(t(data_matrix))
}

#' @export
writedfs <- function(filename, s1) {

  ftype_header <- 'DFS_LE v2.00'
  hdrsize <- 184L
  mdoffset <- 0L        # Start of metadata.
  pdoffset <- 0L       # Start of patient data header.
  nTriangles <- dim(s1$faces)[2]
  nVertices <- dim(s1$vertices)[2]
  nStrips <- 0L
  stripSize <- 0L
  normals <- 0L
  uvoffset <- 0L
  vcoffset <- 0L
  precision <- 0L
  labelOffset <- 0L
  attributes <- 0L
  nextarraypos <- hdrsize + 12L * (nTriangles + nVertices)  # Start fields after the header

  if (!identical(s1$normals, NULL)) {
    normals <- nextarraypos
    nextarraypos <- nextarraypos + nVertices * 12L #12 bytes per normal vector (3 x float32)
  }
  if (!identical(s1$vColor, NULL)) {
    vcoffset <- nextarraypos
    nextarraypos <- nextarraypos + nVertices * 12L # 12 bytes per color coordinate (3 x float32)
  }
  if (!identical(s1$u, NULL) && !identical(s1$v, NULL)) {
    uvoffset <- nextarraypos
    nextarraypos <- nextarraypos + nVertices * 8L # 8 bytes per uv coordinate (2 x float32)
  }
  if (!identical(s1$labels, NULL)) {
    #print 'has labels'
    labelOffset <- nextarraypos
    nextarraypos <- nextarraypos + nVertices * 2L  # 4 bytes per label (int16)
  }
  if (!identical(s1$attributes, NULL)) {
    #print 'has attr'
    attributes <- nextarraypos
    nextarraypos <- nextarraypos + nVertices * 4L  #  4 bytes per attribute (float32)
  }

  dfs_fid <- file(filename, "wb")
  writeChar(ftype_header, dfs_fid, 12, eos = NULL)
  writeBin(hdrsize, dfs_fid)
  writeBin(mdoffset, dfs_fid)
  writeBin(pdoffset, dfs_fid)
  writeBin(nTriangles, dfs_fid)
  writeBin(nVertices, dfs_fid)
  writeBin(nStrips, dfs_fid)
  writeBin(stripSize, dfs_fid)
  writeBin(normals, dfs_fid)
  writeBin(uvoffset, dfs_fid)
  writeBin(vcoffset, dfs_fid)
  writeBin(labelOffset, dfs_fid)
  writeBin(attributes, dfs_fid)

  writeBin(rep(0L, 4L + 15L * 8L), dfs_fid, size = 1L)

  writeBin(as.vector(s1$faces), dfs_fid)
  writeBin(as.vector(s1$vertices), dfs_fid, size = 4)

  if (normals > 0)
    writeBin(as.vector(s1$normals), dfs_fid, size = 4)
  if (vcoffset > 0)
    writeBin(as.vector(s1$vColor), dfs_fid, size = 4)
  if (uvoffset > 0) {
    uv <- as.vector(as.matrix(rbind(s1$u, s1$v), nrow=2, ncol=s1$hdr$nVertices, byrow = FALSE))
    writeBin(uv, dfs_fid, size = 4)
  }
  if (labelOffset > 0)
    writeBin(s1$labels, dfs_fid, size = 2L)
  if (attributes > 0)
    writeBin(s1$attributes, dfs_fid, size = 4)

  close(dfs_fid)
}
