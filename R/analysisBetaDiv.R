#' Beta Diversity Analysis
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param factor Factor along which to conduct beta-diversity analysis
#' @param stratifier (Optional) A second factor along which to stratify groups
#' @param method Method for dissimilarity calculation. One of either "bray",
#'     "euclidean", "jaccard", "unifrac", "manhattan", "canberra", "clark",
#'     "kulczynski", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#'     "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord"
#' @param weighted (Optional) If method is set to "unifrac", whether to perform
#'     weighted or unweighted unifrac. Defaults to FALSE
#'
#' @return List of data tables containing various beta-diversity results
#' @export
#'
bdiv <- function(dataset=NULL, factor=NULL, stratifier=NULL, method='bray', weighted=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- factor[factor %in% names(dataset$factors)]
  if(!length(factor)) factor <- dataset$active_factor

  stratifier <- stratifier[stratifier %in% names(dataset$factors)[!(names(dataset$factors) %in% factor)]]

  results <- pnova(dataset,
                   dist=method, weighted=weighted,
                   factors=c(factor, stratifier))

  return(results)
}

#' Unifrac Distance Calculation
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param weighted Whether to perform weighted or unweighted unifrac.
#'     Defaults to FALSE
#' @param normalized Normalize dissimilarities to a range of 0 to 1? Defaults to
#'     TRUE
#'
#' @return "dist" object of unifrac dissimilarity values between samples
#' @export
#'
mvunifrac <- function(dataset=NULL,weighted=F,normalized=T) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  ps <- makePS(dataset)

  unifrac_dist <- phyloseq::UniFrac(ps,weighted=weighted,normalized=normalized)

  return(unifrac_dist)
}
