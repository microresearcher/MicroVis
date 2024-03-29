% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sampleFiltering.R
\name{addSamples}
\alias{addSamples}
\title{Add Samples}
\usage{
addSamples(dataset = NULL, samples)
}
\arguments{
\item{dataset}{MicroVis dataset (mvdata object)}

\item{samples}{Vector of names of the samples to be added back in. These must
be samples previously removed by removeSamples(). Low-quality samples will
not be added back in using this function.}
}
\value{
MicroVis dataset (mvdata object) with an updated "ignored" lists of
    samples for runSampleFilter() to remove
}
\description{
Add Samples
}
