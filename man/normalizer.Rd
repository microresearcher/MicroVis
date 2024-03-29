% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalization.R
\name{normalizer}
\alias{normalizer}
\title{Normalizing Function}
\usage{
normalizer(
  data,
  norm_type,
  sum_scale = NA,
  log_base = NA,
  impute_method = c("GBM", "SQ", "BL", "CZM")
)
}
\arguments{
\item{data}{Abundance table that will be normalized by in columns. In other
words, if scaling by samples, then samples should be in columns. Does not
matter for data transformation.}

\item{norm_type}{Normalization type: "sum", "none", "relative", "standardized",
"centered", "median", "range", "glog" (generalized log), "pseudolog", "log",}

\item{sum_scale}{Number to scale samples/features to.}

\item{log_base}{Base for log transformations.}

\item{impute_method}{Method to impute zeros by if using zCompositions approach}
}
\value{
Normalized abundance table.
}
\description{
Actual function that normalizes data, one method at a time. In
    other words, it is called three times by runNormalization() if the data
    is to be scaled by samples, scaled by features, and transformed.
}
