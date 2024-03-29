% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-gcms.R
\docType{data}
\name{gcms}
\alias{gcms}
\title{GCMS data on cell differentiation}
\format{
An object of class \code{matrix} (inherits from \code{array}) with 36 rows and 15 columns.
}
\source{
\href{https://www.ebi.ac.uk/metabolights/MTBLS283}{Metabolights}
}
\usage{
data(dataset3)
}
\description{
GCMS data on cell differentiation
}
\examples{
data(dataset3)
mb_d3 <- BuildMultiBlock(gcms)
}
\author{
Cabrero et al. (2019) Scientific data 6:256
(\href{https://doi.org/10.1038/s41597-019-0202-7}{DOI})
}
\keyword{dataset}
