% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analysisDeseq.R
\name{mvdeseq}
\alias{mvdeseq}
\title{DESeq2 Analysis}
\usage{
mvdeseq(dataset = NULL, factor = NULL, rank = NULL, dataset_name = NULL)
}
\arguments{
\item{dataset}{MicroVis dataset. Defaults to the active dataset}

\item{factor}{Factor along which to perform analysis. Defaults to active factor}

\item{rank}{Rank at which to select features to analyze. Defaults to active rank}

\item{dataset_name}{(Not recommended) Name of the dataset to save results
to. This should not need to be used by users since the function can
determine the name of the dataset directly passed to it, but not when
it is called within another function.}
}
\value{
Deseq analysis results
}
\description{
DESeq2 Analysis
}