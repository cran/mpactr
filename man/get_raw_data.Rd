% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mpactr-accessors.R
\name{get_raw_data}
\alias{get_raw_data}
\title{Return the input peak table from mpactr object.}
\usage{
get_raw_data(mpactr_object)
}
\arguments{
\item{mpactr_object}{The mpactr object that is created by calling the
import_data() function.}
}
\value{
a \code{data.table}.
}
\description{
\code{get_raw_data} a wrapper function to return the meta data object of the
given mpactr object.
}
\examples{
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

raw_data <- get_raw_data(data)

}
