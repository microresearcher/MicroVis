% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featureFiltering.R
\name{filterNAs}
\alias{filterNAs}
\title{Filter Unassigned Taxa}
\usage{
filterNAs(dataset = NULL, keepNAs = F, ranks = NULL, temp = F, silent = F)
}
\arguments{
\item{dataset}{MicroVis dataset (mvdata object)}

\item{keepNAs}{Logical telling whether or not to keep taxa that aren't assigned
at the ranks specified}

\item{ranks}{Any taxa that have an unknown assignment at this/these rank(s)
will be filtered out.}

\item{temp}{If set to TRUE, it tells processDataset() to NOT update the active
dataset.}

\item{silent}{Argument that is ultimately passed onto runSampleFilter(),
runNormalization(), and runFeatureFilter(), telling them not to output
any messages.}
}
\value{
MicroVis dataset (mvdata object) with an updated list of taxa to be
    removed because they were unassigned at specified ranks
}
\description{
Filter out taxa that are unassigned at certain ranks
}
