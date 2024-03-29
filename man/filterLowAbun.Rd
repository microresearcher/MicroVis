% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featureFiltering.R
\name{filterLowAbun}
\alias{filterLowAbun}
\title{Filter Low Abundance Features}
\usage{
filterLowAbun(
  dataset = NULL,
  min_abun = NULL,
  min_proportion = NULL,
  silent = F
)
}
\arguments{
\item{dataset}{MicroVis dataset (mvdata object)}

\item{min_abun}{(Optional) (Default = 1) Abundance threshold. Only features
with at least this much abundance in min_proportion percent of samples
will be KEPT.}

\item{min_proportion}{(Optional) (Default = 20) Low abundance proportion
threshold. Only features with at least min_abun abundance in this percentage
of samples will be KEPT.}

\item{silent}{Argument that is ultimately passed onto runSampleFilter(),
runNormalization(), and runFeatureFilter(), telling them not to output
any messages.}
}
\value{
MicroVis dataset (mvdata object) with updated low abundance filtering
    parameters.
}
\description{
Filter Low Abundance Features
}
