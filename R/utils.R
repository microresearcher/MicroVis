#' Save Dataset in its Current State
#'
#' @param dataset_name (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#'
#' @return MicroVis dataset (mvdata object)
#' @export
#'
mvsave <- function(dataset_name=NULL) {
  if(is.null(dataset_name)) dataset_name <- 'active_dataset'

  dataset <- get(dataset_name,envir = mvEnv)
  if(class(dataset)!='mvdata') {
    message('"',dataset_name,'" is not a MicroVis dataset object')
  }

  print(dataset)

  return(dataset)
}

#' Clean up a metadata factor column in a merged data table
#'
#' @param merged_data Datatframe of merged metadata and abundance table
#' @param factor Factor to "clean" by ensuring its column is a factor with only
#'     levels that the filtered dataset contains
#'
#' @return Merged metadata+abundance table with the metadata column of the
#'     specified factor re-leveled to only include levels that exist in the
#'     column
#' @export
#'
cleanData <- function(merged_data, factor) {
  # Make sure only groups in the factor's subset are included
  clean_data <- merged_data[merged_data[[factor$name]] %in% factor$subset,]

  # Reset the metadata column for the factor in this merged table so that its
  #   levels are only the groups in the factor's subset.
  #   This is so that ggplot doesn't get confused looking for levels that don't
  #   don't actually exist in that column
  clean_data[[factor$name]] <- factor(clean_data[[factor$name]],
                                      levels=factor$subset)

  return(clean_data)
}

#' Clean Metadata Factors
#'
#' @param melted_data Dataframe with metadata and abundance data.
#' @param factor_name Factor to refactor in metadata.
#' @param verbose If TRUE (default), print out the count of samples in each group
#'     of the factor and warn the user if fewer than 3 samples are in a group.
#'
#' @return A nested list containing the sample count for each group in the factor.
#'
cleanGroups <- function(melted_data, factor_name, verbose=T) {
  # List of the groups currently in the factor of the melted data
  included_grps <- unique(melted_data[[factor_name]])

  for(grp in included_grps) {
    # If computing statistics, need to make sure that only groups with
    #   at least 3 samples are kept
    sample_size <- nrow(melted_data[melted_data[[factor_name]]==grp,])
    if(sample_size<3) {
      included_grps <- included_grps[!(included_grps %in% grp)]
      if(verbose) {
        message('\n WARNING: ',grp,' (',factor_name,')',' has only ',sample_size,' sample(s), and will be excluded from this analysis and figure\n')
      }
    }
  }
  if(length(included_grps)<2) {
    message('\nERROR: ',length(included_grps),' group(s) remaining with samples size of at least 3 -- analysis and plotting will now terminate\n')
  }

  clean_data <- melted_data[melted_data[[factor_name]] %in% included_grps,]

  clean_data[[factor_name]] <- factor(clean_data[[factor_name]],
                                      levels=included_grps)

  return(clean_data)
}

#' Get Lowest Rank in a Dataset
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#'
#' @return Lowest rank of the data in a dataset.
#' @export
#'
getLowestRank <- function(dataset=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  ds_ranks <- getRanks(dataset)

  lowest_rank <- rev(ds_ranks)[[1]]

  return(lowest_rank)
}

#' Get All Ranks of a Dataset
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#'
#' @return List of ranks in a dataset.
#' @export
#'
getRanks <- function(dataset=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  if(dataset$features!='taxa') return('functional')
  else if(dataset$features=='merged') return('single_rank')

  ds_ranks <- colnames(dataset$data$taxa_names)

  return(ds_ranks)
}

#' Get list of features in a dataset
#'
#' @param dataset Dataset of interest. Default is the currently active dataset
#' @param ranks Vector of ranks at which to get the list of features
#' @param allRanks If set to TRUE, will get the list of all features at all ranks
#'
#' @return List of features in a dataset at specified ranks
#' @export
#'
getFeatures <- function(dataset=NULL,ranks=NULL,allRanks=F) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  if(allRanks) ranks <- getRanks(dataset)
  else {
    ranks <- ranks[ranks %in% getRanks(dataset)]
    if(is.null(ranks)) ranks <- dataset$data$proc$active_rank
  }

  features <- c()
  for(rank in ranks) features <- c(features,colnames(dataset$data$proc[[rank]]))

  features <- features[!(features %in% 'Other')]

  return(features)
}

#' List Significant Features of a Dataset in Current State
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#' @param factor (Optional) Factor along which to perform statistical analysis.
#'     Default is the active factor of the dataset
#' @param ranks (Optional) Ranks at which to get the significant features. Default
#'     is the active rank of the dataset
#' @param allRanks If set to TRUE, will get significant features at all ranks
#'     of the dataset
#' @param alpha (Optional) Significance threshold. Default is 0.05.
#' @param silent If set to TRUE, output text will not be printed
#' @param dataset_name (Optional) Name of the dataset being passed so that results
#'     can be saved to that dataset
#'
#' @return List of significant features.
#' @export
#'
listsigs <- function(dataset=NULL,factor=NULL,
                     ranks=NULL, allRanks=F,
                     alpha=0.05,
                     silent=F,
                     dataset_name=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else if(is.null(dataset_name)) {
    dataset_name <- deparse(substitute(dataset))
  }

  print(dataset)

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  if(allRanks) {
    ranks <- getRanks(dataset)
    # No domain abundance table since analyses don't usually need to be done at
    #   the domain level
    ranks <- ranks[ranks!='domain']
  } else {
    ranks <- ranks[ranks %in% getRanks(dataset)]
    if(is.null(ranks)) ranks <- dataset$data$proc$active_rank
  }

  sigs <- c()
  for(rank in ranks) {
    fts <- colnames(dataset$data$proc[[rank]])
    merged_abun <- mvmelt(dataset)

    if(!length(dataset$stats[[factor]][[rank]]$stats)) dataset <- univar(dataset=dataset,
                                                                         factor=factor,
                                                                         rank=rank,
                                                                         features=fts,
                                                                         dataset_name=dataset_name)

    stats <- dataset$stats[[factor]][[rank]]$stats
    sigs <- c(sigs,stats$.y.[stats$p.adj<=alpha])

    if(!silent) cat('\n\nThese are the significant features at the',rank,'rank:\n',
                    paste(stats$.y.[stats$p.adj<=alpha],collapse = '\n '))
  }

  if(!silent) cat('\n\n')
  return(sigs)
}

#' Standardized naming of the stratification of an analysis
#'
#' @param stratifiers Factors by which an analysis is stratified
#'
#' @return Stratified name as a string
#' @export
#'
nameStratification <- function(stratifiers) {
  return(paste0('by_',paste0(stratifiers,collapse = '_and_')))
}
