% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bss_stat_ops.R
\name{bstr_corr}
\alias{bstr_corr}
\title{Test for Correlation between a variable \code{corr_var} and a brain imaging phenotype.}
\usage{
bstr_corr(corr_var, bstr_data, group_var = "", mult_comp = "fdr")
}
\arguments{
\item{corr_var}{Character variable name. This should be present in the demographics csv file associated
with \code{bstr_data}.}

\item{bstr_data}{Object of type \code{\link[=BstrData]{BstrData()}}.}

\item{group_var}{Character variable for groups in the data. Currently this argument is only used in plotting for ROI analysis.}

\item{mult_comp}{method for multiple comparisons correction. The default method is "fdr". See \code{\link[=bstr_p_adjust]{bstr_p_adjust()}} for valid values.}
}
\description{
Test for correlation between a brain imaging phenotype (cortical thickness, determinant
of the jacobian matrix) and \code{corr_var} using the Pearson's product moment correlation
coefficient. The brain imaging phenotype is automatically selected from the type of \code{bstr_data}.
}
\details{
\code{bstr_data} can be of the type "sba", "tbm", or "roi".
}
