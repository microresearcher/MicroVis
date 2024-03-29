% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featureFiltering.R
\name{clearFeatureFilt}
\alias{clearFeatureFilt}
\title{Clear Feature Filtering}
\usage{
clearFeatureFilt(dataset = NULL, temp = F, silent = F)
}
\arguments{
\item{dataset}{MicroVis dataset (mvdata object)}

\item{temp}{If set to TRUE, it tells processDataset() to NOT update the active
dataset.}

\item{silent}{Argument that is ultimately passed onto runSampleFilter(),
runNormalization(), and runFeatureFilter(), telling them not to output
any messages.}
}
\value{
MicroVis dataset (mvdata object) with all feature filtering parameters
    cleared. This function then calls processDataset() which will call the
    3 "run___" functions and therefore add all the features back to the dataset
    via runFeatureFilter()
}
\description{
Clear Feature Filtering
}
