% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComDim_PCA_MB.R, R/ComDim_PLS_MB.R
\name{Compress_Data_2020}
\alias{Compress_Data_2020}
\title{Compress large multi-block objects.}
\usage{
Compress_Data_2020(s_n = s_n, CompMethod = CompMethod, Partitions = Partitions)

Compress_Data_2020(s_n = s_n, CompMethod = CompMethod, Partitions = Partitions)
}
\arguments{
\item{s_n}{The multi-block object.}

\item{CompMethod}{It can be 'Normal' (default), 'Kernel', 'PCT', 'Tall' or 'Wide'.}

\item{Partitions}{The number of partitions.}
}
\value{
The compressed multi-block.

The compressed multi-block.
}
\description{
Internal function of ComDim_PCA().

Internal function of ComDim_PCA().
}
\examples{
b1 <- matrix(rnorm(5000),10,500)
b2 <- matrix(rnorm(2000),10,200)
blist <- list(b1 = b1, b2 = b2)
blist <- Compress_Data_2020(blist, 'Normal', 1)
b1 <- matrix(rnorm(5000),10,500)
b2 <- matrix(rnorm(2000),10,200)
blist <- list(b1 = b1, b2 = b2)
blist <- Compress_Data_2020(blist, 'Normal', 1)
}
