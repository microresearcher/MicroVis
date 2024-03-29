% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cclasso.R
\name{mvCClasso}
\alias{mvCClasso}
\title{CClasso Network Analysis}
\usage{
mvCClasso(dataset = NULL, r_cutoff = 0, features1 = NULL, features2 = NULL)
}
\arguments{
\item{dataset}{MicroVis dataset. Defaults to the active dataset}

\item{r_cutoff}{Coefficient cutoff. Defaults to 0 (none)}

\item{features1}{(Optional) Features to consider for set 1}

\item{features2}{(Optional) Features to consider for set 2}
}
\value{
Correlation matrices
}
\description{
CClasso Network Analysis
}
