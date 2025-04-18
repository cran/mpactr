% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quality_control.R
\name{plot_qc_tree}
\alias{plot_qc_tree}
\title{Visualize Filtering Summary as Tree Map}
\usage{
plot_qc_tree(mpactr_object)
}
\arguments{
\item{mpactr_object}{an \code{mpactr_object}.}
}
\value{
a tree map plot of class \code{ggplot}.
}
\description{
\code{plot_qc_tree()} visualizes the filtering summary as a treemap. Ion
status (see \code{\link[=qc_summary]{qc_summary()}}) is reported here as percentage of all
pre-filtered ions.
}
\examples{
data <- import_data(
  example_path("coculture_peak_table.csv"),
  example_path("metadata.csv"),
  format = "Progenesis"
)

data_filter <- filter_mispicked_ions(data,
  ringwin = 0.5,
  isowin = 0.01,
  trwin = 0.005,
  max_iso_shift = 3,
  merge_peaks = TRUE
)

plot_qc_tree(data_filter)
}
