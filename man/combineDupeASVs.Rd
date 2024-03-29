% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeASVtab.R
\name{combineDupeASVs}
\alias{combineDupeASVs}
\title{Combine Identical ASVs}
\usage{
combineDupeASVs(taxa_data_list, combineDupes = T)
}
\arguments{
\item{taxa_data_list}{Taxa data with taxa names and un-processed abundance table.}

\item{combineDupes}{If set to TRUE (default), MicroVis will try to combine
duplicate features.}
}
\value{
Taxa data with identical taxa combined
}
\description{
Combine Identical ASVs
}
