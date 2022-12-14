% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComDim_PCA.R
\name{RowsPartition}
\alias{RowsPartition}
\title{Calculate horizontal partitions in the multi-block.}
\usage{
RowsPartition(Xn = Xn, Partitions = Partitions)
}
\arguments{
\item{Xn}{A block.}

\item{Partitions}{The number of partitions.}
}
\value{
The partitioned block.
}
\description{
Internal function of ComDim_PCA().
}
\examples{
b1 <- matrix(rnorm(500),10,50)
b1 <- RowsPartition(b1,2)
}
