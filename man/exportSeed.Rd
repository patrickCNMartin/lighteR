% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExportAndPlot.R
\name{exportSeed}
\alias{exportSeed}
\title{Export data from seed object}
\usage{
exportSeed(
  seed,
  file = NULL,
  dataType = c("retain", "dropped", "traits", "measures"),
  extension = ".csv",
  sep = ","
)
}
\arguments{
\item{seed}{a seed object}

\item{file}{filename to be used - note that this should only be the first part of the file name as extenssion will be added seperately}

\item{dataType}{character string describing which data to export. "retain","dropped","traits","measures"}

\item{extension}{file extension to be used}

\item{sep}{data seperator}
}
\value{
Create files in specified directory containing extracted data.
}
\description{
Export data from seed object
}
