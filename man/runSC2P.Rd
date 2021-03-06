% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runDE.R
\name{runSC2P}
\alias{runSC2P}
\title{Run DE analysis by using SC2P. Here we output two result tables corresponding to two forms of DE genes.}
\usage{
runSC2P(sce)
}
\arguments{
\item{sce}{is a simulated scRNA-seq dataset with two-group conditions, e.g., treatment vs control.}
}
\value{
a list of three tables: the first table summaries the DE result for both forms of DE genes. cont table represents the result for continous case. disc table shows the result for discontinous case.
}
\description{
Run DE analysis by using SC2P. Here we output two result tables corresponding to two forms of DE genes.
}
\examples{
de = runSC2P(simData$sce)
}
