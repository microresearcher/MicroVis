% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDeseq.R
\name{plotDeseq}
\alias{plotDeseq}
\title{Plot DESeq Results}
\usage{
plotDeseq(
  dataset = NULL,
  factor = NULL,
  rank = NULL,
  byGroup = T,
  alpha = 0.05
)
}
\arguments{
\item{dataset}{MicroVis dataset. Defaults to the active dataset}

\item{factor}{Factor along which to analyze data}

\item{rank}{Rank at which to select features}

\item{byGroup}{Group features by group? Defaults to TRUE}

\item{alpha}{Significance threshold. Defaults to 0.05}
}
\value{
MicroVis dataset containing DESeq results
}
\description{
Plot DESeq Results
}
