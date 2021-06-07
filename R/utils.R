#' Save Dataset in its Current State
#'
#' @param name New name to save the currently active dataset
#'     to. If none is specified, simply returns the dataset
#'
#' @return MicroVis dataset (mvdata object)
#' @export
#'
mvsave <- function(name=NULL) {
  dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(name)) {
    dataset$name <- NA
    return(dataset)
  }

  dataset$name <- name
  name <- str_replace_all(name, '[- +=!@#$%^&*()]', '_')

  assign(name, dataset, pos=1)

  print(dataset)
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
#' @param alpha (Optional) Significance threshold. Defaults to 0.05
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

#' List features uniquely over/under-expressed in a specified group
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor along which to perform univariate analysis. Defaults to
#'     the active factor
#' @param groups Group(s) to identify unique features in. Cycles through all
#'     groups by default
#' @param rank Rank at which to analyze features
#' @param alpha Significance threshold. Defaults to 0.05
#' @param param Perform parametrized or nonparametrized univariate analysis?
#'     Defaults to FALSE (nonparametrized)
#' @param dataset_name (Not recommended) Name of the dataset to save statistics
#'     to. This should not need to be used by users since the function can
#'     determine the name of the dataset directly passed to it, but not when
#'     it is called within another function.
#'
#' @return A grouped list of features unique to each group of a factor in a dataset
#' @export
#'
listUniques <- function(dataset=NULL,factor=NULL,groups=NULL,
                        rank=NULL,
                        alpha=0.05,
                        param=F,
                        dataset_name=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else if(is.null(dataset_name)) {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  if(length(dataset$factors[[factor]]$subset)<3) {
    stop('There must be 3 or more groups')
  }

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  if(is.null(dataset$stats[[factor]][[rank]])) dataset <- univar(dataset=dataset,
                                                                 factor=factor,
                                                                 rank=rank,
                                                                 param=param,
                                                                 dataset_name=dataset_name)

  pw_stats <- dataset$stats[[factor]][[rank]]$pw_stats

  discriminators.bygroup <- list()
  for(grp in dataset$factors[[factor]]$subset) {
    discriminators.bygroup[[grp]] <- getDiscriminatingFeatures(pw_stats,grp,alpha=alpha)
  }

  return(discriminators.bygroup)
}

#' Get features that uniquely identify one group from other groups
#'
#' @param pairwise_statistics Data frame of pairwise statistics for a comparison
#'     of 3 or more groups
#' @param unique_group The group of interest -- features that are uniquely
#'     differentially expressed in this group compared to the comparison groups
#'     will be identified
#' @param comparison_groups The groups to compare the group of interest to.
#'     Defaults to all other groups in the subsetted dataset
#' @param alpha Significance threshold. Defaults to 0.05
#'
#' @return List of features uniquely differentially expressed in "unique_group"
#' @export
#'
getDiscriminatingFeatures <- function(pairwise_statistics,
                                      unique_group,
                                      comparison_groups=NULL,
                                      alpha=0.05) {
  if(!all(c('.y.','group1','group2','p.adj') %in% colnames(pairwise_statistics))) {
    stop('"pairwise_statistics" argument must be a table with .y., group1, group2, and p.adj columns')
  }

  all_groups <- unique(unlist(pairwise_statistics[c('group1','group2')]))

  unique_group <- unique_group[unique_group %in% all_groups]
  if(is.null(unique_group)) unique_group <- all_groups[1]

  comparison_groups <- comparison_groups[comparison_groups %in%
                                           all_groups[!(all_groups %in% unique_group)]]
  if(is.null(comparison_groups)) comparison_groups <- all_groups[!(all_groups %in% unique_group)]

  # Filter the pairwise_statistics for only those that:
  #   1) Meet the alpha threshold
  #   2) Have the unique_group (group of interest) in either the "group1" or
  #     "group2" column
  #   3) Have one the comparison_groups in either the "group1" or "group2" column
  pw_sigs <- pairwise_statistics[pairwise_statistics$p.adj<=alpha &
                                   (pairwise_statistics$group1==unique_group |
                                      pairwise_statistics$group2==unique_group) &
                                   (pairwise_statistics$group1 %in% comparison_groups |
                                      pairwise_statistics$group2 %in% comparison_groups),]

  # For each feature, if the number of rows equals the number of comparison groups,
  #   then that means the group of interest (unique_group) is different from all
  #   of the comparison groups, and is thus unique
  discriminating_fts <- sapply(unique(pw_sigs$.y.),
                               function(x) {
                                 ifelse(nrow(pw_sigs[pw_sigs$.y.==x,])==length(comparison_groups),
                                        T,F)
                               })

  return(names(discriminating_fts)[discriminating_fts])
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
