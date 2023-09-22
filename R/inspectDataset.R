#' View the active dataset in RStudio Viewer
#'
#' @return NULL
#' @export
#'
viewactive <- function() {
  active_dataset <- get('active_dataset',envir = mvEnv)
  print(active_dataset)
  View(as.list(active_dataset))
}

#' Get Number of Included Samples in a Dataset
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#' @param factors Factors along which samples will be split into groups.
#'     Default is the active factor of the dataset.
#' @param stratifiers (Optional) Factors along which groups will be stratified.
#'     Default is none.
#' @param getSizes Whether or not to return a nested list with the sizes of each
#'     groups. Default is FALSE.
#' @param min_n Passed along to countSamples.base (non-exported function)
#' @param verbose (Optional) If TRUE, will print out sizes of each group in each
#'     factor.
#'
#' @return If "getSizes" flag is TRUE, will report group sizes. If not, returns
#'     nothing.
#' @export
#'
countSamples <- function(dataset=NULL, factors=NULL, stratifiers=NULL, getSizes=F, min_n=3, verbose=T) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset', envir=mvEnv)
    dataset_name <- 'the active dataset'
  } else {
    dataset_name <- paste0('"',dataset$name,'"')
  }

  factors <- factors[factors %in% names(dataset$factors)]
  if(!length(factors)) factors <- names(dataset$factors)

  stratifiers <- stratifiers[stratifiers %in% names(dataset$factors) & !(stratifiers %in% factors)]

  metadata <- mvmelt(dataset)[1:ncol(dataset$metadata)]

  size_report <- countSamples.base(metadata=metadata,
                                   factors=factors,
                                   stratifiers=stratifiers,
                                   min_n=min_n,
                                   dataset=dataset,
                                   verbose=verbose)

  # Lastly, list any ignored samples
  if(verbose) {
    if(length(dataset$data$proc$ignored_samples)) {
      cat('\n  Ignored samples:\n')
      cat(paste(dataset$data$proc$ignored_samples,collapse = '\t'))
    }
    cat('\n')
  }

  if(getSizes) return(size_report)
}

#' Core function for calculating sample sizes
#'
#' @param metadata Metadata of the dataset
#' @param factors Factors along which samples will be counted in groups
#' @param stratifiers Factors along which groups will be stratified
#' @param min_n Threshold at which to warn user of a small sample size in a group
#'     or a stratified group.
#' @param dataset (Optional) Dataset is passed to relay information about factors that is
#'     otherwise lost.
#' @param verbose Whether or not to print out the size report for the groups/
#'     stratified groups.
#'
#' @return Returns a nested list of the sample size for each group/stratified group
#'
countSamples.base <- function(metadata,
                              factors,
                              stratifiers=NULL,
                              min_n=3,
                              dataset=NULL,
                              verbose=T) {
  metadata$Size <- rep(1,nrow(metadata))
  sample_names <- metadata$sample

  if(verbose) cat(paste0('\n>>> ',length(sample_names),' total samples:'))

  size_report <- list()

  # Only print out unstratified group sizes if no stratifiers were specified
  for(f in factors) {
    if(verbose & !length(stratifiers)) cat(paste0('\n  by ',f,':'))

    grp_sizes <- aggregate(as.formula(paste0('Size ~ ',f)), metadata[c(f,'Size')],sum)
    if(verbose & !length(stratifiers)) for(grp in grp_sizes[[f]]) {
      cat(paste0('\n    ',grp,' contains ',grp_sizes$Size[grp_sizes[[f]]==grp],' samples'))
    }

    if(!is.null(dataset)) {
      zero_grps <- dataset$factors[[f]]$grps[!(dataset$factors[[f]]$grps %in% grp_sizes[[f]])]
      if(verbose & length(zero_grps) & !length(stratifiers)) {
        cat(paste0('\n\n    No samples remaining in ',
                   paste0(zero_grps,collapse = ', ')))
      }
    }

    size_report[[f]] <- grp_sizes
    if(verbose & !length(stratifiers)) cat('\n')
  }

  # Count groups in each stratum (up to 2 levels of stratification), if stratifying
  if(length(stratifiers)) {
    for(f in factors) {
      strata_sizes <- aggregate(as.formula(paste0('Size ~ ',
                                                  paste0(c(f,stratifiers),collapse = '+'))),
                                metadata[c(f,stratifiers,'Size')],sum)
      if(verbose) {
        cat(paste0('\n  by ',f,' and stratified by ',
                   paste0(stratifiers,collapse = ' and '),':'))

        for(grp in unique(strata_sizes[[f]])) {
          cat(paste0('\n    ',grp,' contains ',sum(strata_sizes$Size[strata_sizes[[f]]==grp]),
                     ' samples:'))
          cat(paste0('\n      ',strata_sizes$Size[strata_sizes[[f]]==grp],' ',
                     paste0(interaction(strata_sizes[stratifiers],
                                        sep = ' and ')[strata_sizes[[f]]==grp])))
        }
      }
      # Identify any strata that are zero in at least one group
      possible_strata <- levels(interaction(strata_sizes[c(f,stratifiers)],sep = ' and '))
      present_strata <- as.character(interaction(strata_sizes[c(f,stratifiers)],sep = ' and '))
      zero_strata <- possible_strata[!(possible_strata %in% present_strata)]
      if(verbose & length(zero_strata)) cat(paste0('\n\n    No samples remaining in:\n      ',
                                                   paste0(zero_strata,collapse = '\n      ')))

      # Add strata with zero samples to the summary table
      for(stratum in zero_strata) {
        temp <- data.frame(t(strsplit(zero_strata,' and ')[[1]]),0)
        colnames(temp) <- c(f,stratifiers,'Size')
        strata_sizes <- rbind(strata_sizes,temp)
      }

      size_report$stratified[[f]] <- strata_sizes
      if(verbose) cat('\n')
    }
  }

  return(size_report)
}

#' Get Number of Included Samples in a Dataset
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#'
#' @return Number of features in a dataset at its active rank.
#' @export
#'
countFeatures <- function(dataset=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset', envir=mvEnv)
    dataset_name <- 'the active dataset'
  } else {
    dataset_name <- paste0('"',dataset$name,'"')
  }

  rank <- dataset$data$proc$active_rank
  num_features <- ncol(dataset$data$proc[[rank]])
  return(num_features)
}

#' Get read depth summary statistics for each group or for all samples
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#' @param factor (Optional) Factor along which samples will be grouped
#' @param individual (Optional) Whether to return read depths for each individual
#'     sample. Defaults to FALSE
#' @param all (Optional) Whether to calculate read depth statistics for all the
#'     samples as one group
#'
#' @return Dataframe of read depth statistics including min, max, median, and mean
#' @export
#'
readdepth <- function(dataset=NULL, factor=NULL, individual=F, all=F) {
  if(is.null(dataset)) dataset <- get('active_dataset', envir=mvEnv)

  if(!is.null(dataset$name)) dataset_name <- paste0('"',dataset$name,'"')
  else dataset_name <- 'active_dataset'

  taxaRanks <- get('taxaRanks', envir = mvEnv)
  highest_rank <- taxaRanks[taxaRanks %in% getRanks(dataset)][[1]]
  abun.unproc <- mvmelt(clearProcessing(dataset, temp = T, silent = T),
                        rank = highest_rank)

  mdcols <- colnames(dataset$metadata)
  fts <- colnames(abun.unproc)[!(colnames(abun.unproc) %in% mdcols)]

  samplereads <- data.frame(abun.unproc[mdcols],
                            'Total_Reads'=rowSums(abun.unproc[fts]))

  if(individual) {
    return(samplereads)
  } else if(all) {
    readsummary <- data.frame(Min=min(samplereads$Total_Reads),
                              Max=max(samplereads$Total_Reads),
                              Median=median(samplereads$Total_Reads),
                              Mean=mean(samplereads$Total_Reads),
                              StDev=sd(samplereads$Total_Reads))
  } else {
    factor <- factor[factor %in% names(dataset$factors)]
    if(is.null(factor)) factor <- dataset$active_factor

    readsummary <- samplereads %>%
      dplyr::group_by(get(factor)) %>%
      summarise(Min=min(Total_Reads),
                Max=max(Total_Reads),
                Median=median(Total_Reads),
                Mean=mean(Total_Reads),
                StDev=sd(Total_Reads))
  }

  return(readsummary)
}


#' Get Data with/without Attached Metadata
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#' @param rank (Optional) Rank to get abundance values of. Default is active rank.
#' @param group This isn't used and can be removed
#' @param metadata (Optional) Whether or not to attach metadata.
#'
#' @return Abundance data at specified rank with/without metadata.
#' @export
#'
getdata <- function(dataset=NULL, rank=NULL, group=NULL, metadata=T) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset', envir=mvEnv)
    dataset_name <- 'the active dataset'
  } else {
    dataset_name <- paste0('"',dataset$name,'"')
  }
  rank <- rank[tolower(rank) %in% tolower(c(get('taxaRanks',envir = mvEnv),'single_rank'))]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  data <- dataset$data$proc[[rank]]
  if(metadata) {
    data$sample <- rownames(data)
    data <- merge(dataset$metadata,data,by='sample')
  }

  return(data)
}

#' View Abundance Table
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#'
#' @return Does not return anything.
#' @export
#'
viewtable <- function(dataset=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset', envir=mvEnv)
    dataset_name <- 'the active dataset'
  } else {
    dataset_name <- paste0('"',dataset$name,'"')
  }
  rank <- dataset$data$proc$active_rank

  cat('\nBringing up abundance table for',dataset_name,'at',rank,'level\n\n')

  View(dataset$data$proc[[rank]])
}

#' View Statistics of All Features
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#' @param type Type of statistical test to view results for. Defaults to the first
#'     available of univariate, deseq, and lefse in that order
#' @param factor (Optional) Factor to count samples by groups. Default is the
#'     active factor of the dataset.
#' @param rank (Optional) Rank to get abundance values of. Default is active rank.
#'
#' @return View of the statistics of the dataset
#' @export
#'
viewstats <- function(dataset=NULL, type=NULL, factor=NULL, rank=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else {
    dataset_name <- dataset$name
    dataset <- get(dataset$name, pos = 1)
  }

  factor <- names(dataset$factors)[tolower(names(dataset$factors)) %in% tolower(factor)]
  if(!length(factor)) factor <- dataset$active_factor

  rank <- rank[rank %in% getRanks(dataset)]
  if(!length(rank)) rank <- dataset$data$proc$active_rank

  stat_types <- names(dataset$stats[[factor]])
  if(is.null(stat_types)) {
    stop('No analyses have been performed for ',dataset_name)
  }

  if(is.null(type)) type <- stat_types[1]
  else if(is.null(type[type %in% stat_types])) {
    stop('Please choose one of the following for "type":\n ',
         paste(stat_types, collapse=', '))
  }

  cat('\nPulling up',type,'statistical results for',factor,'at',rank,'rank\n\n')

  stats <- dataset$stats[[factor]][[type]][[rank]]
  if(type=='univar') View(stats$stats)
  else View(stats)
}

#' Check Colors of Groups in a Dataset
#'
#' @param dataset (Optional) MicroVis dataset (mvdata object). If not specified,
#'     defaults to active dataset.
#' @param assigned (Optional) Whether or not to show what groups the colors
#'     are assigned to. Default is FALSE.
#'
#' @return Plot of colors with/without the associated groups.
#' @export
#'
checkColors <- function(dataset=NULL, assigned=FALSE) {
  defCols <- get('defCols',envir = mvDefaults)
  clrs <- dataset$colors

  if(assigned) {
    dummydf <- data.frame(x=names(clrs),
                          y=rep(1,length(clrs))
    )
    ggplot(data=dummydf, aes(x=x,y=y,fill=x))+
      geom_bar(stat='identity')+
      scale_fill_manual(values=clrs)
  } else {
    dummydf <- data.frame(x=letters[1:length(defCols)],
                          y=rep(1,length(defCols))
    )
    ggplot(data=dummydf, aes(x=x,y=y,fill=x))+
      geom_bar(stat='identity')+
      scale_fill_manual(values=defCols)
  }
}
