% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chooseFactors.R
\name{rangetotext}
\alias{rangetotext}
\title{Translate Range Factor Levels into Interpretable Text}
\usage{
rangetotext(factor)
}
\arguments{
\item{factor}{Factor (list object) with its groups to be translated}
}
\value{
New factor (list object) with name appended with "_Range" and group
    names changed to the the formats as appropriate: "# or less", "# to #",
    or "# or more"
}
\description{
Translate Range Factor Levels into Interpretable Text
}
