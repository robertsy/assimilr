% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_plotting.R
\name{image_mat}
\alias{image_mat}
\title{Plot covariance}
\usage{
image_mat(A, xlab = "columns", ylab = "rows", plotlegend = TRUE, ...)
}
\arguments{
\item{A}{the matrix to plot}

\item{xlab, ylab, }{plotlegend plotting options}

\item{...}{additional parameters passed to image.plot=>image}
}
\description{
plot a squared matrix as a heatmap
take care of putting things in the right direction and add a legend
}
