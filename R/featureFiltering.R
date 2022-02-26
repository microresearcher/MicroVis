#' @title Feature Filter
#'
#' @description This function is one of the processing functions called by
#'     processDataset() and identifies and removes features based on parameters
#'     set by the filter___ feature filtering functions
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
runFeatureFilter <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(!is.null(dataset$data$proc$rarefied)) return(dataset)

  if(is.null(dataset$data$proc$unranked)) dataset <- runNormalization(dataset, silent = T)

  if(is.null(dataset$data$proc$filter_rank)) dataset$data$proc$filter_rank <- getLowestRank(dataset)

  filtering <- dataset$data$proc$filtering

  # If there is something in the dataset's filtering history besides
  #   "filterlist" and "ftstats", proceed with filtering
  if(length(names(filtering)[!(names(filtering) %in% c('filterlist','ftstats'))])) {
    if(!silent) cat('\n\n|~~~~~~~~~~~~~  FILTERING FEATURES  ~~~~~~~~~~~~~|\n')

    filter_rank <- dataset$data$proc$filter_rank

    ft_data <- getFtStats(dataset,justStats = F)
    abd_temp <- ft_data$proc$unranked
    if(dataset$features=='taxa') abd_temp <- agglomTaxa(ft_data, abd_temp,
                                                        from_rank='asv',
                                                        to_rank=filter_rank)

    nfts <- nrow(ft_data$proc$filtering$ftstats)

    filterlist <- list()
    ftstats <- ft_data$proc$filtering$ftstats

    ### Identify Low Prevalence ###
    #-----------------------------#
    if(!is.null(filtering$top_prevalence)) {
      low_prevalence <- (ftstats %>% slice_min(Prevalence,
                                               n=(nfts-filtering$top_prevalence)))[['Feature']]

      if(!silent) cat(paste0('\n  Identified top ',
                             filtering$top_prevalence,' features by prevalence'))

      filterlist$low_prevalence <- low_prevalence

    } else if(!is.null(filtering$min_prevalence)) {
      low_prevalence <- ftstats[ftstats$Prevalence<filtering$min_prevalence,]$Feature

      if(!silent) cat(paste0('\n  Identified ',
                             length(low_prevalence),' features with < ',
                             filtering$min_prevalence,' prevalence'))

      filterlist$low_prevalence <- low_prevalence
    }
    ### Identify Low Relative Abundance ###
    #-------------------------------------#
    if(!is.null(filtering$top_relabun)) {
      low_relabun <- (ftstats %>% slice_min(Mean_Relative_Abundance,
                                            n=(nfts-filtering$top_relabun)))[['Feature']]

      if(!silent) cat(paste0('\n  Identified top ',
                             filtering$top_relabun,' features by relative abundance'))

      filterlist$low_relabun <- low_relabun

    } else if(!is.null(filtering$min_relabun)) {
      low_relabun <- ftstats[ftstats$Mean_Relative_Abundance<filtering$min_relabun,]$Feature

      if(!silent) cat(paste0('\n  Identified ',
                             length(low_relabun),' features with < ',
                             filtering$min_relabun,' relative abundance'))

      filterlist$low_relabun <- low_relabun
    }
    ### Identify Low Total Abundance ###
    #----------------------------------#
    if(!is.null(filtering$top_totabun)) {
      low_totabun <- (ftstats %>% slice_min(Total_Abundance,
                                            n=(nfts-filtering$top_totabun)))[['Feature']]

      if(!silent) cat(paste0('\n  Identified top ',
                             filtering$top_totabun,' features by total abundance'))

      filterlist$low_totabun <- low_totabun

    } else if(!is.null(filtering$min_totabun)) {
      low_totabun <- ftstats[ftstats$Total_Abundance<filtering$min_totabun,]$Feature

      if(!silent) cat(paste0('\n  Identified ',
                             length(low_totabun),' features with < ',
                             filtering$min_totabun,' total abundance'))

      filterlist$low_totabun <- low_totabun
    }
    ### Identify Low Variance ###
    #---------------------------#
    if(!is.null(filtering$top_var)) {
      low_var <- (ftstats %>% slice_min(StDev, n=(nfts-filtering$top_var)))[['Feature']]

      if(!silent) cat(paste0('\n  Identified the top ',
                             filtering$top_var,' features by standard deviation'))

      filterlist$low_var <- low_var

    } else if(!is.null(filtering$low_var_percentile)) {
      low_sd <- slice_min(data.frame(sd=unique(ftstats$StDev)),
                          sd, n=floor(filtering$low_var_percentile*nrow(ftstats)/100))$sd

      low_var <- ftstats[ftstats$StDev %in% low_sd,]$Feature

      if(!silent) cat(paste0('\n  Identified ',
                             length(low_var),' features with variance in the bottom ',
                             filtering$low_var_percentile,' percentile by standard deviation'))

      filterlist$low_variance <- low_var
    }
    ### Identify Low Abundance ###
    #----------------------------#
    if(!is.null(filtering$low_abun$min_abun) & !is.null(filtering$low_abun$min_prop)) {
      ftstats$LowAbunProp=apply(abd_temp, 2,
                                function(x) sum(x<filtering$low_abun$min_abun)/nrow(abd_temp))

      low_abun <- ftstats[ftstats$LowAbunProp>(1-filtering$low_abun$min_prop/100),]$Feature

      if(get('keepSigFisher',envir = mvEnv)) {
        low_abun <- low_abun[!(low_abun %in% findSigFisher(dataset,low_abun,silent=silent))]
      }

      if(!silent) cat(paste0('\n  Identified ',
                             length(low_abun),' features with < ',
                             filtering$low_abun$min_abun,' abundance in >',
                             filtering$low_abun$min_prop,'% of samples'))

      filterlist$low_abun <- low_abun
    }

    taxaranks <- c('species','genus','family','order','class','phylum','single_rank','pathways')
    lowest_rank <- getLowestRank(dataset)

    ### Identify Taxa with NAs ###
    #----------------------------#
    # If keepNA exists and is false, identify NAs at the active rank
    taxa_names_tab <- ft_data$taxa_names
    active_rank <- ft_data$proc$active_rank
    if(!is.null(filtering$NAfilter$ranks)) {
      if(!is.null(taxa_names_tab)) if(ncol(taxa_names_tab)>1) for(rank in filtering$NAfilter$ranks) {
          filterlist$NAs <- unique(taxa_names_tab[grep(paste0('_',rank),
                                                       taxa_names_tab[[rank]]),][[lowest_rank]])
      }

      if(!silent) cat(paste0('\n  Identified ',
                             length(filterlist$NAs),' features without assigned ',
                             paste0(filtering$NAfilter$ranks,collapse = ', ')))
    } else filterlist$NAs <- NULL

    # Reload the normalized ft_data and abundance table
    ft_data <- dataset$data
    abd_temp <- ft_data$proc$unranked

    if(length(unlist(filterlist))) {
      ### Filter all Identified Features ###
      #------------------------------------#
      # First, temporarily separate the NA list from the nested filterlist as
      #   NAs will be dealt with separately
      unklist <- filterlist$NAs
      filterlist$NAs <- NULL

      # Replace the taxa names with all ASV numbers corresponding to them
      if(!is.null(taxa_names_tab)) {
        filterlist.ids <- TaxatoASV(ft_data,
                                    unique(unlist(filterlist)),
                                    filter_rank)
        unk.ids <- TaxatoASV(ft_data,
                             unique(unklist),
                             filter_rank)
      } else {
        filterlist.ids <- unique(unlist(filterlist))
        unk.ids <- unique(unklist)
      }

      # Put the NA list back in filterlist since we won't be using filterlist
      #   for the rest of this
      filterlist$NAs <- unklist

      # Remove any ASVs that are not in the current unranked abundance table either
      #   because they were excluded prior to processing or during normalization.
      #   They are unnecessary to have in filterlist.ids or unk.ids and they will
      #   cause issues
      filterlist.ids <- filterlist.ids[filterlist.ids %in% colnames(abd_temp)]
      unk.ids <- unk.ids[unk.ids %in% colnames(abd_temp)]

      # Make a new abundance table of just the features that made it through
      #   the filter
      abd_filtered <- abd_temp[!(colnames(abd_temp) %in% c(filterlist.ids,unk.ids))]

      # Now, pool all the filter list features into 'Other' and add it to the table
      abd_filtered$Other <- rowSums(abd_temp[filterlist.ids])
      abd_filtered$Unknown <- rowSums(abd_temp[unk.ids])
    } else abd_filtered <- abd_temp

    ft_data$proc$unranked <- abd_filtered

    # Record the filter list in the dataset's history
    ft_data$proc$filtering$filterlist <- filterlist

    num_removed <- ncol(abd_temp) - ncol(abd_filtered) + any(colnames(abd_filtered) %in% c('Other','Unknown'))

    if(!silent) cat('\n\n>>> Removed',num_removed,'features based on filtering parameters <<<\n')
  } else if(length(dataset$data$proc$selected)) {
    ft_data <- dataset$data
    abd_temp <- ft_data$proc$unranked
    selected <- ft_data$proc$selected

    if(dataset$features=='taxa') {
      if(is.null(names(selected))) selected.ids <- TaxatoASV(dataset$data, selected, 'single_rank')
      else selected.ids <- unique(unlist(lapply(getRanks(dataset),
                                                function(rank) {
                                                  if(length(selected[[rank]])) {
                                                    TaxatoASV(dataset$data,
                                                              selected[[rank]],
                                                              rank)
                                                  }
                                                })))
    } else selected.ids <- unique(unlist(selected))

    abd_selected <- abd_temp[selected.ids]
    # Now, pool all the filter list features into 'Other' and add it to the table
    abd_selected$Other <- rowSums(abd_temp[!(colnames(abd_temp) %in% selected.ids)])

    ft_data$proc$unranked <- abd_selected

    num_selected <- length(selected.ids)

    if(!silent) cat('\n\n>>> Selected',num_selected,'features in total <<<\n')
  } else {
    ft_data <- dataset$data
    if(!silent) cat('\n~~~ No feature filtering performed ~~~\n')
  }

  ft_data <- makeRankTabs(ft_data)
  # No need for the unranked table anymore, as all downstream data manipulation
  #   (i.e. normalization methods) will be performed separately on each ranked
  #   abundance table
  ft_data$proc$unranked <- NULL
  dataset$data <- ft_data

  return(dataset)
}

#' Clear Feature Filtering
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp If set to TRUE, it tells processDataset() to NOT update the active
#'     dataset.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with all feature filtering parameters
#'     cleared. This function then calls processDataset() which will call the
#'     3 "run___" functions and therefore add all the features back to the dataset
#'     via runFeatureFilter()
#'
#' @export
#'
clearFeatureFilt <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  dataset$data$proc$filtering <- NULL

  dataset <- processDataset(dataset, temp=temp, silent=silent)

  return(dataset)
}

#' Filter Low Prevalence Features
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param top (Optional) Top n number of features to be kept based on prevalence.
#'     Specifying a value for this will override the min_prevalence value.
#' @param min_prevalence (Optional) (Default) Prevalence threshold. Features with
#'     prevalence lower than this threshold will be filtered out.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with updated low prevalence filtering
#'     parameters.
#'
#' @export
#'
filterLowPrev <- function(dataset=NULL, top=NULL, min_prevalence=NULL, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(length(dataset$data$proc$selected)) {
    message('\nRemoving selected feature list')
    dataset$data$proc$selected <- NULL
  }

  filtering.defaults <- get('filtering.defaults',envir = mvDefaults)

  dataset$data$proc$rarefied <- NULL

  nfts <- ncol(dataset$data$orig)
  ft_data <- dataset$data

  if(is.null(top) & is.null(min_prevalence)){
    if(!is.null(ft_data$proc$filtering$top_prevalence)) {
      top <- ft_data$proc$filtering$top_prevalence
    } else if(!is.null(ft_data$proc$filtering$min_prevalence)) {
      min_prevalence <- ft_data$proc$filtering$min_prevalence
    } else {
      min_prevalence <- filtering.defaults$min_prevalence
    }
  } else {
    if(!is.null(top)) {
      if(top > 0 & top < nfts) min_prevalence <- NULL
      else top <- NULL
    }
    if(!is.null(min_prevalence)) {
      if(min_prevalence <= 0) min_prevalence <- NULL
    }
  }

  # Load these parameters to the dataset
  ft_data$proc$filtering$top_prevalence <- top
  ft_data$proc$filtering$min_prevalence <- min_prevalence

  dataset$data <- ft_data

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, silent=silent)

  return(dataset)
}

#' Filter Low Relative Abundance Features
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param top (Optional) Top n number of features to be kept based relative
#'     abundance. Specifying a value for this will override the min_relabun
#'     value.
#' @param min_relabun (Optional) (Default) Relative abundance threshold. Features
#'     with prevalence lower than this threshold will be filtered out.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with updated low relative abundance
#'     filtering parameters.
#'
#' @export
#'
filterLowRelAbun <- function(dataset=NULL, top=NULL, min_relabun=NULL, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(length(dataset$data$proc$selected)) {
    message('\nRemoving selected feature list')
    dataset$data$proc$selected <- NULL
  }

  filtering.defaults <- get('filtering.defaults',envir = mvDefaults)

  dataset$data$proc$rarefied <- NULL

  nfts <- ncol(dataset$data$orig)
  ft_data <- dataset$data

  if(is.null(top) & is.null(min_relabun)){
    if(!is.null(ft_data$proc$filtering$top_relabun)) {
      top <- ft_data$proc$filtering$top_relabun
    } else if(!is.null(ft_data$proc$filtering$min_relabun)) {
      min_relabun <- ft_data$proc$filtering$min_relabun
    } else {
      min_relabun <- filtering.defaults$min_relabun
    }
  } else {
    if(!is.null(top)) {
      if(top > 0 & top < nfts) min_relabun <- NULL
      else top <- NULL
    }
    if(!is.null(min_relabun)) {
      if(min_relabun <= 0) min_relabun <- NULL
    }
  }

  # Whatever threshold is being set for this run will be stored in the dataset's filtering "history"
  ft_data$proc$filtering$top_relabun <- top
  ft_data$proc$filtering$min_relabun <- min_relabun

  dataset$data <- ft_data

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, silent=silent)

  return(dataset)
}

#' Filter Low Total Abundance Features
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param top (Optional) Top n number of features to be kept based on total abundance.
#'     Specifying a value for this will override the min_totabun value.
#' @param min_totabun (Optional) (Default) Total abundance threshold. Features
#'     with total abundance lower than this threshold will be filtered out.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with updated low total abundance
#'     filtering parameters.
#'
#' @export
#'
filterLowTotAbun <- function(dataset=NULL, top=NULL, min_totabun=NULL, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(length(dataset$data$proc$selected)) {
    message('\nRemoving selected feature list')
    dataset$data$proc$selected <- NULL
  }

  filtering.defaults <- get('filtering.defaults',envir = mvDefaults)

  dataset$data$proc$rarefied <- NULL

  nfts <- ncol(dataset$data$orig)
  ft_data <- dataset$data

  if(is.null(top) & is.null(min_totabun)){
    if(!is.null(ft_data$proc$filtering$top_totabun)) {
      top <- ft_data$proc$filtering$top_totabun
    } else if(!is.null(ft_data$proc$filtering$min_totabun)) {
      min_totabun <- ft_data$proc$filtering$min_totabun
    } else {
      min_totabun <- filtering.defaults$min_totabun
    }
  } else {
    if(!is.null(top)) {
      if(top > 0 & top < nfts) min_totabun <- NULL
      else top <- NULL
    }
    if(!is.null(min_totabun)) {
      if(min_totabun <= 0) min_totabun <- NULL
    }
  }

  # Whatever threshold is being set for this run will be stored in the dataset's filtering "history"
  ft_data$proc$filtering$top_totabun <- top
  ft_data$proc$filtering$min_totabun <- min_totabun

  dataset$data <- ft_data

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, silent=silent)

  return(dataset)
}

#' Filter Low Variance Features
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param top (Optional) Top n number of features to be kept based on relative
#'     abundance. Specifying a value for this will override the min_relabun value.
#' @param var_type Method for assessing variance of samples. Current choices are
#'     standard deviation ('sd') and interquartile range ('iqr')
#' @param low_var_percentile (Optional) (Default) Variance threshold. Features
#'     with variance (standard deviation) lower than this threshold will be
#'     filtered out.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with updated low variance filtering
#'     parameters.
#'
#' @export
#'
filterLowVar <- function(dataset=NULL, top=NULL, var_type='sd', low_var_percentile=NULL, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(length(dataset$data$proc$selected)) {
    message('\nRemoving selected feature list')
    dataset$data$proc$selected <- NULL
  }

  filtering.defaults <- get('filtering.defaults',envir = mvDefaults)

  dataset$data$proc$rarefied <- NULL

  nfts <- ncol(dataset$data$orig)
  ft_data <- dataset$data

  if(is.null(top) & is.null(low_var_percentile)){
    if(!is.null(ft_data$proc$filtering$top_var)) {
      top <- ft_data$proc$filtering$top_var
    } else if(!is.null(ft_data$proc$filtering$low_var_percentile)) {
      low_var_percentile <- ft_data$proc$filtering$low_var_percentile
    } else {
      low_var_percentile <- filtering.defaults$low_var_percentile
    }
  } else {
    if(!is.null(top)) {
      if(top > 0 & top < nfts) low_var_percentile <- NULL
      else top <- NULL
    }
    if(!is.null(low_var_percentile)) {
      if(low_var_percentile <= 0) low_var_percentile <- NULL
    }
  }

  # Whatever threshold is being set for this run will be stored in the
  #   dataset's filtering "history"
  ft_data$proc$filtering$top_var <- top
  ft_data$proc$filtering$low_var_percentile <- low_var_percentile

  dataset$data <- ft_data

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, silent=silent)

  return(dataset)
}

#' Filter Low Abundance Features
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param min_abun (Optional) (Default = 1) Abundance threshold. Only features
#'     with at least this much abundance in min_proportion percent of samples
#'     will be KEPT.
#' @param min_proportion (Optional) (Default = 20) Low abundance proportion
#'     threshold. Only features with at least min_abun abundance in this percentage
#'     of samples will be KEPT.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with updated low abundance filtering
#'     parameters.
#'
#' @export
#'
filterLowAbun <- function(dataset=NULL, min_abun=NULL, min_proportion=NULL, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(length(dataset$data$proc$selected)) {
    message('\nRemoving selected feature list')
    dataset$data$proc$selected <- NULL
  }

  filtering.defaults <- get('filtering.defaults',envir = mvDefaults)

  dataset$data$proc$rarefied <- NULL

  ft_data <- dataset$data

  # Is a new threshold being specified?
  if(is.null(min_abun)) {
    # If not, has this filtering been done previously?
    if(is.null(ft_data$proc$filtering$low_abun$min_abun)) {
      # If not, set the threshold to the default value
      min_abun <- filtering.defaults$low_abun$min_abun
    } else {
      # If this filtering has been done previously, set the threshold to the previously used one
      #   since no threshold was specified during this call of the function
      min_abun <- ft_data$proc$filtering$low_abun$min_abun
    }
  } else if(min_abun==0) {
    # If the threshold is being set to 0, set it to NULL so it disappears
    #   in the dataset's filtering history
    min_abun <- NULL
    min_proportion <- 0
  }

  # Is a new threshold being specified?
  if(is.null(min_proportion)) {
    # If not, has this filtering been done previously?
    if(is.null(ft_data$proc$filtering$low_abun_prop)) {
      # If not, set the threshold to the default value
      min_proportion <- filtering.defaults$low_abun$min_prop
    } else {
      # If this filtering has been done previously, set the threshold to the previously used one
      #   since no threshold was specified during this call of the function
      min_proportion <- ft_data$proc$filtering$low_abun$min_prop
    }
  } else if(min_proportion==0) {
    # If the threshold is being set to 0, set it to NULL so it disappears
    #   in the dataset's filtering history
    min_proportion <- NULL
    min_abun <- NULL
  }

  # Whatever threshold is being set for this run will be stored in the dataset's filtering "history"
  if(is.null(c(min_abun,min_proportion))) ft_data$proc$filtering$low_abun <- NULL
  else ft_data$proc$filtering$low_abun <- list(min_abun=min_abun,
                                               min_prop=min_proportion)

  dataset$data <- ft_data

  # Now process the dataset
  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, silent=silent)

  return(dataset)
}

#' @title Filter Unassigned Taxa
#'
#' @description Filter out taxa that are unassigned at certain ranks
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param keepNAs Logical telling whether or not to keep taxa that aren't assigned
#'     at the ranks specified
#' @param ranks Any taxa that have an unknown assignment at this/these rank(s)
#'     will be filtered out.
#' @param temp If set to TRUE, it tells processDataset() to NOT update the active
#'     dataset.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with an updated list of taxa to be
#'     removed because they were unassigned at specified ranks
#' @export
#'
filterNAs <- function(dataset=NULL, keepNAs=F, ranks=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(length(dataset$data$proc$selected)) {
    message('\nRemoving selected feature list')
    dataset$data$proc$selected <- NULL
  }

  if(dataset$features!='taxa') return(dataset)

  if(keepNAs) dataset$data$proc$filtering$NAfilter <- NULL
  else {
    # If not keeping NAs, then by default remove any taxa in just the highest
    #   rank of the dataset
    dataset$data$proc$filtering$NAfilter$ranks <- colnames(dataset$data$taxa_names)[1]

    ranks <- tolower(ranks)[tolower(ranks) %in% getRanks(dataset)]

    if(!length(ranks)) ranks <- dataset$data$proc$filtering$NAfilter$ranks
    else dataset$data$proc$filtering$NAfilter$ranks <- ranks
  }

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, temp=temp, silent=silent)

  return(dataset)
}

#' Select Specific Features
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param features Vector of features (at any rank) to select
#' @param ranks (Optional) Specify a rank from which to select features. Searches
#'     through all ranks by default
#' @param temp This parameter has no use in this function and can be removed
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return Dataset with list of selected features that is passed to runFeatureFilter
#' @export
#'
selectFeatures <- function(dataset=NULL, features, ranks=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(length(names(dataset$data$proc$filtering)[!(names(dataset$data$proc$filtering)
                                                 %in% c('filterlist','ftstats'))])) {
    if(temp) removefilts <- T
    else removefilts <- ifelse(select.list(c('Yes','No'), title='\nSelecting specific features will undo any filtering. Continue?')=='Yes', yes = T, no = T)

    if(removefilts) dataset$data$proc$filtering <- NULL
    else stop('Cannot select specific features from a dataset that has been filtered')
  }

  if(dataset$features=='taxa') alltaxa <- unique(unname(unlist(dataset$data$taxa_names)))
  else alltaxa <- colnames(dataset$data$orig)

  features <- unique(alltaxa[tolower(alltaxa) %in% tolower(features)])
  if(!length(features)) stop('None of the features were found in "',dataset$name,'"')

  if(dataset$features=='taxa') {
    ft_names <- dataset$data$taxa_names
    if(ncol(ft_names)>1) selected <- sapply(getRanks(dataset),
                                            function(rank) features[features %in% ft_names[[rank]]])
    else selected <- list('single_rank'=features)
  } else selected <- list('functional'=features)

  ranks <- ranks[ranks %in% getRanks(dataset)]
  if(!is.null(ranks)) selected[!(names(selected) %in% ranks)] <- NULL

  dataset$data$proc$selected <- selected

  dataset <- processDataset(dataset, temp=temp, silent=silent)

  return(dataset)
}

#' Undo Feature Selection
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param temp This parameter has no use in this function and can be removed
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset with the selected feature list cleared
#' @export
#'
unselectFeatures <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  dataset$data$proc$selected <- NULL

  dataset <- processDataset(dataset, temp=temp, silent=silent)

  return(dataset)
}

#' @title Run Fisher Test on Features
#'
#' @description Run Fisher test on features that are labeled for filtering. This
#'     test is run on presence-absence data and not abundance data.
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param fts List of features to test.
#' @param lowabun_thresh Abundance value above which a feature is considered to
#'     be present in a given sample.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return List of features that had a statistically significant Fisher test
#'     when performed across the subsetted groups of the active factor.
#'
#' @export
#'
findSigFisher <- function(dataset, fts, lowabun_thresh=0, silent=F) {
  if(is.null(dataset$data$proc$unranked)) dataset <- runSampleFilter(dataset,temp = T,silent = T)

  abun <- dataset$data$proc$unranked

  if(dataset$features=='taxa') {
    abun <- agglomTaxa(dataset$data, abun, from_rank = 'asv',
                       to_rank=dataset$data$proc$filter_rank)
  }
  active_factor <- dataset$active_factor

  if(!silent) cat(paste0('\n Running Fisher test on low abundance features...'))

  sig_fisher <- c()
  for(ft in fts) {
    abd <- data.frame(sample=rownames(abun),
                      present=as.numeric(abun[[ft]]>lowabun_thresh),
                      absent=as.numeric(abun[[ft]]<=lowabun_thresh))
    df <- merge(dataset$metadata, abd)

    form <- formula(paste('. ~',active_factor))
    xtab <- aggregate(form, df[c(active_factor,'present','absent')], function(x) sum(x))

    if(fisher_test(xtab[2:3])$p<0.05) sig_fisher <- c(sig_fisher,ft)
  }
  return(sig_fisher)
}

#' Get Feature Statistics for Filtering
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param rank Rank at which to get feature summary stats
#' @param raw Use raw counts to calculate feature statistics. Defaults to TRUE
#' @param justStats Whether to just output the a table of feature statistics
#'     instead of all data processing information
#' @param doCLR Calculate average CLR-transformed abundances for each feature
#'
#' @return List containing original abundance table, feature names table, and
#'     a dataframe of feature statistics used for filtering
#'
#' @export
#'
getFtStats <- function(dataset=NULL, rank=NULL, raw=T, justStats=T, doCLR=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  # If a pre-filtered, unranked abundance table is not available, run sample
  #   filtering or normalization to get one
  if(is.null(dataset$data$proc$unranked)) {
    if(raw) dataset <- runSampleFilter(dataset,temp = T,silent = T)
    else dataset <- runNormalization(dataset,temp = T,silent = T)
  }

  if(is.null(rank)) rank <- dataset$data$proc$filter_rank

  ft_data <- dataset$data
  abd <- ft_data$proc$unranked

  if(!is.null(ft_data$taxa_names)) colnames(abd) <- ASVtoTaxa(ft_data, colnames(abd), rank)

  # Combine any features with the same name at this rank (used for agglomerating taxa at higher ranks)
  abd <- t(rowsum(t(abd),group = colnames(abd)))

  # Compute prevalence, total abundance, proportional prevalence, and relative abundance of each feature
  ftstats <- data.frame(Prevalence=apply(abd, 2, function(x) sum(x>0,na.rm = T)),
                        Total_Abundance=apply(abd, 2, function(x) sum(x,na.rm = T)),
                        Mean_Abundance=apply(abd, 2, function(x) mean(x,na.rm = T)),
                        Relative_Abundance=apply(apply(abd,1,function(x) x/sum(x,na.rm = T)),
                                                 1,function(x) max(x,na.rm = T)),
                        Mean_Relative_Abundance=apply(apply(abd,1,function(x) x/sum(x,na.rm = T)),
                                                 1,function(x) mean(x,na.rm = T)),
                        StDev=apply(abd, 2, function(x) sd(x,na.rm = T)),
                        IQR=apply(abd, 2, function(x) IQR(x,na.rm = T)))

  # What percentage of all samples is each feature present in
  ftstats$Prevalence_Proportion <- ftstats$Prevalence/nrow(abd)
  ftstats$Feature <- colnames(abd)

  if(doCLR) {
    # Mean CLR-transformed abundance for each feature
    abd.nozeros <- zeroReplace(abd)

    ftstats$Mean_CLR=apply(apply(abd.nozeros, 1, function(x) log(x/exp(mean(log(x))), base=10)),
                           1, function(x) mean(x,na.rm = T))
  }

  ft_data$proc$filtering$ftstats <- ftstats

  if(justStats) return(ftstats)
  else return(ft_data)
}
