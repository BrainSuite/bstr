% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bssdata_class.R
\name{package_data}
\alias{package_data}
\title{Package data for reproducible statistical analysis.}
\usage{
package_data(
  type = "sba",
  subjdir = NULL,
  csv = "",
  hemi = "left",
  sbasmooth = 0,
  tbmsmooth = 0,
  dbasmooth = 0,
  measure = "FA",
  atlas = "",
  eddy = TRUE,
  outdir = NULL,
  exclude_col = ""
)
}
\arguments{
\item{type}{character string denoting type of analysis. Should be sba, tbm, or roi.}

\item{subjdir}{subject directory containing BrainSuite processed data.}

\item{csv}{filename of a comma separated (csv) file containing the subject demographic information.
The first column of this csv file
should be "subjID" and should have subject identifiers you wish to analyze. subjID can be alphanumeric
and should be exactly equal to the individual subject directory name.}

\item{hemi}{chaaracter string denoting the brain hemisphere. Should either be "left" or "right".}

\item{sbasmooth}{numeric value denoting the smoothing level for sba.}

\item{tbmsmooth}{numeric value denoting the smoothing level for tbm.}

\item{dbasmooth}{numeric value denoting the smoothing level for dba.}

\item{measure}{character specifying the brain imaging measure. If analyzing diffusion data, should be "FA".}

\item{atlas}{path name to the atlas}

\item{eddy}{boolean for specifying if the diffusion images were eddy-current corrected or not.}

\item{outdir}{output directory that will contain the copied data}

\item{exclude_col}{character string for the column in demographics csv (contains 1 or 0 for each row) specifying the subjects to exclude. 1 denotes include, 0 denotes exclude.}
}
\description{
Takes same parameters as \code{\link[=load_bstr_data]{load_bstr_data()}} and copies the data to a new directory specified by outdir. You can repeatedly call this function to copy data of different types (tbm -- nii.gz, sba -- .dfs files etc.) to the same output directory.
Prior to using this function, BrainSuite and svreg should be run on all subjects.
If required, smoothing should be performed on cortical surface or volumetric image based measures.
A csv file containing subject demographic information should exist. The first column of this csv file
should have the subject identifiers. Subject identifiers can be alphanumeric
and should be exactly equal to the individual subject directory names.
}
