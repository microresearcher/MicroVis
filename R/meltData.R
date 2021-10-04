#' Melt the metadata and abundance table at a given rank for the desired features
#'
#' @param dataset MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#' @param rank Get abundance table of features at this rank. Defaults to active
#'     rank.
#' @param features Get the abundance data of certain features only. Defaults to
#'     all features.
#' @param min_n Passed to cleanData() to limit to samples that belong to groups
#'     with at least 3 samples
#'
#' @return Table of merged metadata and abundance table
#' @export
#'
mvmelt <- function(dataset=NULL,rank=NULL,features=NULL,min_n=3) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  md <- dataset$metadata

  abun <- dataset$data$proc[[rank]]

  # Any taxa names that start with a number get prepended with "X", so this "X"
  #   character needs to be removed
  colnames(abun) <- lapply(colnames(abun), function(x) {
    if(startsWith(x,'X') & is.numeric(type.convert(substr(x, 2, 2), as.is=T))) {
      sub('X','',x)
    } else {
      x
    }})

  if(!is.null(features)) {
    if(!all(features %in% colnames(abun))) {
      if(!any(features %in% colnames(abun))) {
        message('\nWARNING: None of the specified features were found. Defaulting to all features at ',rank,' rank')
        features <- NULL
      }
      else {
        message('\nWARNING: Some of the specified features were not found.')
        features <- features[features %in% colnames(abun)]
      }
    }
    abun <- abun[features]
  }

  abun$sample <- rownames(abun)
  merged <- merge(md,abun,by='sample')
  for(f in names(dataset$factors)) merged <- cleanData(merged,dataset$factors[[f]])

  return(merged)
}

#' Make a melted data table grouped by strata
#'
#' @param dataset MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#' @param factor Main factor whose groups will be in each grouped data table
#' @param stratifiers Stratifier(s) to group the melted data table by
#' @param rank Rank of dataset to choose features from
#' @param min_n Passed to cleanGroups()
#' @param verbose If set to TRUE, will print processing text
#'
#' @return Grouped, melted data table of metadata and abundance at a specified
#'     rank
#' @export
#'
mvstratify <- function(dataset=NULL, factor=NULL, stratifiers=NULL, rank=NULL, min_n=3, verbose=T) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  if(length(dataset$factors)<2) {stop('Only 1 factor in this dataset. No other factors to stratify by')}

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  stratifiers <- stratifiers[stratifiers %in% names(dataset$factors) & !(stratifiers %in% factor)]
  if(is.null(stratifiers)) {
    stratifiers <- names(dataset$factors)[!(names(dataset$factors) %in% factor)]
    while(length(stratifiers)>2) {
      message('\nPlease pick up to 2 stratifiers')
      stratifiers <- select.list(stratifiers,multiple = T,title = 'Pick up to 2 stratifiers',graphics = T)
    }
  }

  melted <- mvmelt(dataset,rank=rank)

  for(stratum in stratifiers) {
    melted <- melted %>% group_by(get(stratum),.add=TRUE)
    melted <- melted %>% select(-stratum)
    colnames(melted)[ncol(melted)] <- stratum
  }

  return(melted)
}
