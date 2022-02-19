#' @title Feature Remover
#'
#' @description This function removes features prior to any processing and is
#'     one of the processing functions called by processDataset(). It should be
#'     used WITH CAUTION as it can drastically modify a dataset, making it
#'     completely unrepresentative if used incorrectly.
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp This parameter has no use in this function and can be removed
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with features filtered out based on
#'     parameters set by the filter___ feature filtering functions
#'
runFeatureRemover <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  # The first time this is run on a dataset, save the original ASV table in a
  #   table called "raw_input"
  if(is.null(dataset$data$raw_input)) dataset$data$raw_input <- dataset$data$orig

  excluded <- dataset$data$proc$excluded_features

  if(length(excluded)) {
    if(!silent) cat(paste0('\n\n|~~~~~~~~~~~~~  EXCLUDING FEATURES  ~~~~~~~~~~~~~|\n'))

    new_orig <- dataset$data$orig

    new_orig <- new_orig[!(colnames(new_orig) %in% unname(unlist(excluded)))]

    dataset$data$orig <- new_orig

    if(!silent) for(ft in names(excluded)) cat(paste0('\n  Excluding ',
                                                      length(excluded[[ft]]),' ASVs from ',
                                                      ft))
  } else dataset$data$orig <- dataset$data$raw_input

  return(dataset)
}

#' Exclude specific features at a specified taxonomic rank. Only for taxonomic
#'     datasets with multiple taxonomic ranks.
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param features Vector of features (at any rank) to select
#' @param rank (Optional) Specify a rank from which to select features. Searches
#'     through all ranks by default
#' @param temp This parameter has no use in this function and can be removed
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#' @param reset Reset the feature remover to include all features prior to processing
#'
#' @return Dataset with list of excluded features that is passed to runFeatureFilter
#' @export
#'
excludeFeatures <- function(dataset=NULL, features, reset=F, rank=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(reset) {
    dataset$data$proc$excluded_features <- NULL
    dataset <- processDataset(dataset, temp=temp, silent=silent)
    return(dataset)
  }

  current_excluded <- dataset$data$proc$excluded_features

  allranks <- get('taxaRanks',envir = mvEnv)
  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- allranks[(allranks %in% getRanks(dataset))][[1]]

  if(!(rank %in% c('single_rank','functional'))) {
    new_excluded <- sapply(features,
                           function(ft) list(TaxatoASV(dataset$data,ft,taxa_rank=rank)))
    names(new_excluded) <- paste(names(new_excluded),rank)
  }

  if(length(new_excluded)) {
    for(ft in names(new_excluded)) if(!length(new_excluded[[ft]])) {
      message('\nNo ASVs found for a ',ft)
      new_excluded[[ft]] <- NULL
    }
    new_excluded <- new_excluded[!(names(new_excluded) %in% names(current_excluded))]
    dataset$data$proc$excluded_features <- append(current_excluded, new_excluded)
  } else stop('None of the features were found in "',
            dataset$name,
            '" at the ',rank,' rank')

  if(length(dataset$data$proc$selected)) {
    message('\nRemoving selected feature list')
    dataset$data$proc$selected <- NULL
  }

  dataset <- processDataset(dataset, temp=temp, silent=silent)

  return(dataset)
}
