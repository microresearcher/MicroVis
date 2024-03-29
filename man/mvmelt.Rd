% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meltData.R
\name{mvmelt}
\alias{mvmelt}
\title{Melt the metadata and abundance table at a given rank for the desired features}
\usage{
mvmelt(dataset = NULL, rank = NULL, features = NULL, min_n = 3)
}
\arguments{
\item{dataset}{MicroVis dataset (mvdata object). If not specified,
defaults to active dataset.}

\item{rank}{Get abundance table of features at this rank. Defaults to active
rank.}

\item{features}{Get the abundance data of certain features only. Defaults to
all features.}

\item{min_n}{Passed to cleanData() to limit to samples that belong to groups
with at least 3 samples}
}
\value{
Table of merged metadata and abundance table
}
\description{
Melt the metadata and abundance table at a given rank for the desired features
}
