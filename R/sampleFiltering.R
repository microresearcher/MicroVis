#' @title Sample Filter
#'
#' @description This function is one of the processing functions called by
#'     processDataset() removes samples from low quality samples identified by
#'     removeLowQuality(), groups excluded by chooseGrps()/addGrps()/removeGrps()
#'     , and samples specifically excluded by removeSamples()
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp This parameter has no use in this function and can be removed
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runFeatureFilter(), and runNormalizer(). If TRUE, it tells these functions
#'     to NOT print out any of their processing messages.
#'
#' @return MicroVis dataset (mvdata object) with the following samples removed:
#'     1) Low quality samples identified by removeLowQuality()
#'     2) Samples in groups excluded by chooseGrps()/addGrps()/removeGrps()
#'     3) Specific samples removed by the user using removeSamples()
#'
runSampleFilter <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(!silent) cat(paste0('\n\n|~~~~~~~~~~~~~  FILTERING SAMPLES  ~~~~~~~~~~~~~~|\n'))
  # Reference variables
  metadata <- dataset$metadata
  active_factor <- dataset$active_factor
  # End of reference variables

  ft_data <- dataset$data

  excluded <- ft_data$proc$low_quality$low_quality
  if(!silent & length(excluded)) {
    cat(paste0('\nRemoving ',length(excluded),' low quality samples:'))
    for(s in excluded) {
      cat(paste0('\n  Removing sample ',s))
      if(!is.null(active_factor)) {
        cat(paste0(' (from ',as.character(metadata[metadata$sample==s,][[active_factor]]),')'))
      }
    }
    cat('\n')
  }

  factors <- dataset$factors
  for(f in factors) {
    excluded_samples <- list()
    excluded_grps <- f$groups[!(f$groups %in% f$subset)]
    if(!silent & length(excluded_grps)) cat(paste0('\nSubsetting data by ',f$name,':'))
    for(grp in excluded_grps) {
      excluded_samples[[grp]] <- metadata[metadata[[f$name]]==grp,]$sample
      if(!silent & length(excluded_samples[[grp]])) {
        cat(paste0('\n  Removing ',length(excluded_samples[[grp]]),' samples from ',grp))
      }
    }
    excluded_samples <- unname(unlist(excluded_samples))
    excluded <- union(excluded,excluded_samples)
  }
  if(!silent) cat('\n')

  ignored_samples <- ft_data$proc$ignored_samples
  if(!silent & length(ignored_samples)) {
    cat(paste0('\nRemoving ',length(ignored_samples),' samples that were selected for exclusion:'))
    for(s in ignored_samples) {
      cat(paste0('\n  Removing sample ',s))
      if(!silent & !is.null(active_factor)) {
        cat(paste0(' (from ',as.character(metadata[metadata$sample==s,][[active_factor]]),')'))
      }
    }
    cat('\n')
  }
  excluded <- union(excluded,ignored_samples)

  # The sample filter is the first step, and thus starts with the original
  #   abundance table
  abd <- ft_data$orig
  abd_subset <- abd[!(rownames(abd) %in% excluded),]
  ft_data$proc$unranked <- abd_subset

  dataset$data <- makeRankTabs(ft_data)

  # Remove any groups from each factor's subset list if their sample sizes are 0 after pruning
  if(!is.null(factors)) grp_sizes <- countSamples(dataset, getSizes = T, verbose = !silent)
  for(f in factors) {
    grp_stats <- grp_sizes[[f$name]]
    grps <- as.character(grp_stats[[f$name]])
    factors[[f$name]]$subset <- f$groups[f$groups %in% grps]
  }
  dataset$factors <- factors

  # If, after filtering out samples, only one group of the primary factor
  #   has any samples remaining, provide a warning
  if(!is.null(active_factor)) {
    num_grps_for_stats <- sum(grp_sizes[[active_factor]]$Size >= 3)
    if(num_grps_for_stats < 2) {
      message('\nWARNING: Only ',num_grps_for_stats,' group(s) left in the primary factor with at least 3 samples.\nTo run comparative analyses, please rerun filtSamples() with a different reads threshold\n')
      assign('warning_list',c(get('warning_list',envir = mvEnv), paste0('WARNING: Only ',num_grps_for_stats,' group(s) left in the primary factor with at least 3 samples.\nTo run comparative analyses, please rerun filtSamples() with a different reads threshold')),envir = mvEnv)
    }
  }

  return(dataset)
}

#' Remove Low Quality Samples
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runFeatureFilter(), and runNormalizer(). If TRUE, it tells these functions
#'     to NOT print out any of their processing messages.
#'
#' @return MicroVis dataset (mvdata object) with an updated list of low-quality
#'     features for runSampleFilter() to filter out
#'
#' @export
#'
removeLowQuality <- function(dataset=NULL, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(!is.null(dataset$data$proc$low_quality)) {
    rthresh <- dataset$data$proc$low_quality$reads_thresh
  } else rthresh <- get('rthresh',envir = mvDefaults)

  # Process input dataset
  ft_data <- dataset$data
  abd <- ft_data$orig
  sample_names <- rownames(abd)
  metadata <- dataset$metadata
  cmpgrp <- setFVar(dataset)

  samplestats <- data.frame(sample=as.character(sample_names),
                            data.frame(Tot.Reads=apply(abd, 1, function(x) sum(x))))
  samplestats <- samplestats %>% arrange(-Tot.Reads)
  samplestats$sample <- factor(samplestats$sample, levels=samplestats$sample)

  cat(paste0('\n\nDataset has ',length(sample_names),' samples\n'))
  p_samplereads <- ggplot(samplestats, aes(sample, Tot.Reads, fill=Tot.Reads))+
    geom_bar(stat='identity',)+
    geom_hline(yintercept = rthresh, color='coral', size=1)+
    coord_flip()
  print(p_samplereads)

  # Ask user if they would like to change the reads threshold 'rthresh'
  cat(paste0('\n\nCurrent read count threshold: ',rthresh,'\n'))

  changethresh <- ifelse(select.list(choices = c('Keep','Change'),
                                     title = 'Change threshold? Samples with read counts lower than this will be removed from analysis.')=='Change',TRUE,FALSE)

  if(changethresh) {
    rthresh <- as.numeric(readline(prompt = 'Enter a cutoff value for sample read depth: '))
    while(is.na(rthresh)) {
      rthresh <- as.numeric(readline(prompt = 'Please enter a number for the read depth cutoff: '))
    }
    p_samplereads <- ggplot(samplestats, aes(sample, Tot.Reads, fill=Tot.Reads))+
      geom_bar(stat='identity',)+
      geom_hline(yintercept = rthresh, color='coral', size=1)+
      coord_flip()
    print(p_samplereads)
  }

  # Save the read depths image to a png file
  if(get('.loading',envir = mvEnv)) {ggsave(filename = file.path(get('project_dir',envir = mvEnv),
                                                                 'read_depth.png'),height = nrow(samplestats)*.15,
                       plot = p_samplereads,device = 'png')}

  lowqual <- as.character(samplestats[samplestats$Tot.Reads<rthresh,]$sample)

  # Store the reads threshold and low quality sample names in dataset
  ft_data$proc$low_quality <- list(reads_threshold=rthresh,
                                   low_quality=lowqual)

  dataset$data <- ft_data

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, silent = silent)

  return(dataset)
}

#' Choose Groups to be Analyzed in Each Factor
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param factor_names (Optional) Names of the factors for which to specify groups
#'
#' @return MicroVis dataset (mvdata object) with updated "subset" lists for each
#'     of its factors, and an updated list of samples for runSampleFilter() to
#'     remove based on the groups chosen for analysis
#' @export
#'
chooseGrps <- function(dataset=NULL, factor_names=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  # This function alters these variables, so they will need to be passed back to
  #   "dataset" at the end before "dataset" is returned. We use the original
  #   abundance table because other functions permanently remove groups from the
  #   "$proc" abundance tables
  factors <- dataset$factors
  # If dataset contains no factors, then just return the dataset
  if(!length(factors)) return(dataset)

  # If no valid factor names were provided in the function call then default to asking about all factors
  if(!length(factor_names %in% names(factors))) factor_names <- names(factors)

  # These variables are just for reference so they do not need to be passed back
  #   to "dataset"
  metadata <- dataset$metadata

  for(f in factors[names(factors) %in% factor_names]) {
    if(any(grepl('\u2264',f$groups))) {
      isrange <- T
      symbol <- paste(expression('\u2264'))
    } else if(any(grepl('\u2265',f$groups))) {
      isrange <- T
      symbol <- paste(expression('\u2265'))
    } else isrange <- F
    message(paste0('Currently selected groups for "', f$name,'":'), paste0('\n  ', f$subset))
    grps <- select.list(f$groups,
                        multiple = TRUE,
                        title = 'Select groups to analyze in this factor. Press Enter to keep current groups.',
                        graphics = TRUE)

    if(isrange) grps <- unlist(lapply(grps, function(x) gsub('=',symbol,x)))
    # If no groups were selected for a given factor(meaning no samples will be selected)
    #   then default to all the groups in the factor
    if(!length(grps)) {
      grps <- f$subset
      cat(paste0('Keeping the same groups for "',f$name,'"\n\n'))
      # message('WARNING: No groups were selected for ', f$name, ' defaulting to all groups in this factor\n')
    } else if(!all(grps %in% f$groups)) {
      # This is just a safety in case some group names were changed internally
      #   In this case, the subsetted groups will not change
      grps <- f$groups
    }
    factors[[f$name]]$subset <- grps
  }
  dataset$factors <- factors

  dataset <- processDataset(dataset)

  return(dataset)
}


#' Remove Groups
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param factor_name (Optional) Names of the factors for which to specify groups
#' @param grps (Optional) Names of the groups to remove
#'
#' @return MicroVis dataset (mvdata object) with updated "subset" lists for each
#'     of its factors, and an updated list of samples for runSampleFilter() to
#'     remove based on the groups chosen for analysis
#'
#' @export
#'
removeGrps <- function(dataset=NULL, factor_name=NULL, grps=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  # This function alters these variables, so they will need to be passed back to
  #   "dataset" at the end before "dataset" is returned. We use the original
  #   abundance table because other functions permanently remove groups from the
  #   "$proc" tables
  factors <- dataset$factors

  # These variables are just for reference so they do not need to be passed back
  #   to "dataset"
  metadata <- dataset$metadata

  if(length(factors)) {
    for(f in factors) {
      if(any(grepl('\u2264',f$groups))) {
        isrange <- T
        symbol <- paste(expression('\u2264'))
      } else if(any(grepl('\u2265',f$groups))) {
        isrange <- T
        symbol <- paste(expression('\u2265'))
      } else isrange <- F
      exclude <- select.list(f$subset,
                             multiple = TRUE,
                             title = 'Select group(s) to remove',
                             graphics = TRUE)

      if(isrange) exclude <- unlist(lapply(exclude, function(x) gsub('=',symbol,x)))
      grps <- f$subset[!(f$subset %in% exclude)]
      # If no groups were selected for a given factor(meaning no samples will be selected)
      #   then default to all the groups in the factor
      if(!length(grps)) {
        grps <- f$subset
        message('\nWARNING: All groups were removed for ', f$name, ' will not remove any groups')
      } else if(!all(grps %in% f$subset)) {
        # This is just a safety in case some group names were changed internally
        #   In this case, the subsetted groups will not change
        grps <- f$subset
      }
      factors[[f$name]]$subset <- grps
    }
    dataset$factors <- factors
  }

  dataset <- processDataset(dataset)

  return(dataset)
}

#' Add Groups
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param factor_name (Optional) Names of the factors for which to specify groups
#' @param grps (Optional) Names of the groups to remove
#'
#' @return MicroVis dataset (mvdata object) with updated "subset" lists for each
#'     of its factors, and an updated list of samples for runSampleFilter() to
#'     remove based on the groups chosen for analysis
#'
#' @export
#'
addGrps <- function(dataset=NULL, factor_name=NULL, grps=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  # This function alters these variables, so they will need to be passed back to
  #   "dataset" at the end before "dataset" is returned. We use the original
  #   abundance table because other functions permanently remove groups from the
  #   "$proc" tables
  factors <- dataset$factors

  # These variables are just for reference so they do not need to be passed back
  #   to "dataset"
  metadata <- dataset$metadata

  if(length(factors)) {
    for(f in factors) {
      if(any(grepl('\u2264',f$groups))) {
        isrange <- T
        symbol <- paste(expression('\u2264'))
      } else if(any(grepl('\u2265',f$groups))) {
        isrange <- T
        symbol <- paste(expression('\u2265'))
      } else isrange <- F

      # Get the list of groups that aren't currently included
      missing <- f$groups[!(f$groups %in% f$subset)]
      if(!length(missing)) {
        cat('\nNo missing groups in ',f$name)
        next
      }

      include <- select.list(missing,
                             multiple = TRUE,
                             title = 'Select group(s) to remove',
                             graphics = TRUE)

      if(isrange) include <- unlist(lapply(include, function(x) gsub('=',symbol,x)))
      grps <- union(f$subset,include)
      # If no groups were selected for a given factor(meaning no samples will be selected)
      #   then default to all the groups in the factor
      if(!length(grps)) {
        grps <- f$subset
        message('\nWARNING: No groups were added for ', f$name, ' so the current subset will be kept')
      } else if(!all(grps %in% f$subset)) {
        # This is just a safety in case some group names were changed internally
        #   In this case, the subsetted groups will not change
        grps <- f$subset
      }
      factors[[f$name]]$subset <- grps
    }
    dataset$factors <- factors
  }

  dataset <- processDataset(dataset)

  return(dataset)
}

#' Select specific samples to analyze
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param samples List of samples to analyze
#' @param exclude List of samples to exclude. Samples in this list will be excluded even if they are in the "samples" list.
#' @param includeLowQual Include samples even if they are low quality? Defaults
#'     to FALSE
#'
#' @return Dataset with updated list of samples to ignore
#' @export
#'
chooseSamples <- function(dataset=NULL, samples, exclude=c(), includeLowQual=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  samples <- samples[samples %in% dataset$metadata$sample]
  if(!length(samples)) stop('None of the specified samples are in the dataset')

  current_samples <- mvmelt(dataset)$sample

  low_quality <- dataset$data$proc$low_quality$low_quality
  ignored_samples <- dataset$data$proc$ignored_samples

  if(any(samples %in% low_quality)) {
    if(includeLowQual) low_quality <- low_quality[!(low_quality %in% samples)]
    else message('\nNote: The following samples are low quality at a reads threshold of ',
                 dataset$data$proc$low_quality$reads_threshold,' and will not be included:\n',
                 paste0(samples[samples %in% low_quality], collapse = '\t'),
                 '\n\nTo include these samples, re-run this function with `includeLowQual=TRUE`\n')
  }

  ignored_samples <- union(setdiff(current_samples, samples), exclude)

  dataset$data$proc$low_quality$low_quality <- low_quality
  dataset$data$proc$ignored_samples <- ignored_samples

  dataset <- processDataset(dataset)

  return(dataset)
}

#' Remove Specific Samples
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param samples Vector of names of the samples to be removed
#'
#' @return MicroVis dataset (mvdata object) with an updated "ignored" lists of
#'     samples for runSampleFilter() to remove
#'
#' @export
#'
removeSamples <- function(dataset=NULL, samples) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(samples)) {
    message('\nERROR: Please enter a sample name (in quotes) or vector of sample names to exclude')
    return(dataset)
  }

  ft_data <- dataset$data
  rank <- ft_data$proc$active_rank
  remaining_samples <- rownames(ft_data$proc[[rank]])
  ignored_samples <- ft_data$proc$ignored_samples

  excludelist <- c()
  for(s in samples) {
    if(s %in% remaining_samples) {
      ignored_samples <- union(ignored_samples,s)
    } else message(s,' is not currently in the dataset. It may be misspelled or already filtered out')
  }

  ft_data$proc$ignored_samples <- ignored_samples

  dataset$data <- ft_data

  dataset <- processDataset(dataset)

  return(dataset)
}

#' Add Samples
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param samples Vector of names of the samples to be added back in. These must
#'     be samples previously removed by removeSamples(). Low-quality samples will
#'     not be added back in using this function.
#'
#' @return MicroVis dataset (mvdata object) with an updated "ignored" lists of
#'     samples for runSampleFilter() to remove
#'
#' @export
#'
addSamples <- function(dataset=NULL, samples) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(samples)) {
    message('\nERROR: Please enter a sample name (in quotes) or vector of sample names to add back')
    return(dataset)
  }

  ft_data <- dataset$data
  low_quality <- ft_data$proc$low_quality$low_quality
  ignored_samples <- ft_data$proc$ignored_samples

  for(s in samples) {
    if(s %in% low_quality) {
      low_quality <- low_quality[!(low_quality %in% s)]
    }
    if(s %in% ignored_samples) {
      ignored_samples <- ignored_samples[!(ignored_samples %in% s)]
    }
  }
  ft_data$proc$low_quality$low_quality <- low_quality
  ft_data$proc$ignored_samples <- ignored_samples

  dataset$data <- ft_data

  dataset <- processDataset(dataset)

  return(dataset)
}
