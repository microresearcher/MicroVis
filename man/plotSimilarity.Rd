% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSimilarity.R
\name{plotSimilarity}
\alias{plotSimilarity}
\title{Plot similarity matrix between samples}
\usage{
plotSimilarity(
  dataset = NULL,
  rank = NULL,
  dist_method = "bray",
  weighted = F,
  clust_method = "ward.D2",
  clust_num = 2,
  r_cutoff = 0
)
}
\arguments{
\item{dataset}{MicroVis dataset. Defaults to the active dataset}

\item{rank}{Rank of features to use for similarity calculation}

\item{dist_method}{Dissimilarity calculation method. One of either "pearson",
"spearman", bray", "euclidean", "jaccard", "unifrac", "manhattan", "canberra",
"clark", "kulczynski", "gower", "altGower", "morisita", "horn", "mountford",
"raup", "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord".
Defaults to "bray"}

\item{weighted}{If performing unifract distance, whether to use weighted or
unweighted unifrac. Defults to FALSE (unweighted)}

\item{clust_method}{Clustering method. One of either "ward.D", "ward.D2",
"single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
"median" (= WPGMC) or "centroid" (= UPGMC). Defaults to "ward.D2"}

\item{clust_num}{Number of clusters to try to make}

\item{r_cutoff}{R-values with absolute value below this cutoff will be shaded
white. Defaults to 0}
}
\value{
Similarity matrix heatmap
}
\description{
Plot similarity matrix between samples
}
