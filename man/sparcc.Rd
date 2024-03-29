% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparcc.R
\name{sparcc}
\alias{sparcc}
\title{Sparcc wrapper}
\usage{
sparcc(data, iter = 20, inner_iter = 10, th = 0.1)
}
\arguments{
\item{data}{Community count data matrix}

\item{iter}{Number of iterations in the outer loop}

\item{inner_iter}{Number of iterations in the inner loop}

\item{th}{absolute value of correlations below this threshold are considered zero by the inner SparCC loop.}
}
\description{
A reimplementation of SparCC algorithm (Friedman et Alm 2012, PLoS Comp Bio, 2012).
}
