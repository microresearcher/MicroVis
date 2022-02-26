#' Print function for mvdata (MicroVis dataset) objects
#'
#' @param x MicroVis dataset (mvdata object)
#' @param ... Arguments passed to print() function
#'
#' @return Print-friendly printout of dataset details
#' @export
#'
print.mvdata <- function(x, ...) {
  # Get relevant dataset info
  if(is.null(x$name)) saved <- crayon::red$bold$italic('(Not saved)\n')
  else if(is.na(x$name)) saved <- crayon::italic('(Saved)\n')
  else saved <- crayon::italic('(Saved to',paste0('"',x$name,'")'),'\n')

  rank <- x$data$proc$active_rank

  if(x$features=='taxa') feature_type <- rank
  else feature_type <- x$features

  n_samples <- nrow(x$data$proc[[rank]])
  n_ignored_samples <- length(x$data$proc$ignored_samples)

  n_features <- ncol(x$data$proc[[rank]])
  if('Other' %in% colnames(x$data$proc[[rank]])) n_features <- n_features - 1
  if('Unknown' %in% colnames(x$data$proc[[rank]])) n_features <- n_features - 1
  if(x$features=='taxa') totn_features <- length(unique(ASVtoTaxa(x$data,
                                                                  colnames(x$data$orig),
                                                                  taxa_rank=rank)))
  else totn_features <- ncol(x$data$orig)

  factors <- x$factors
  active_factor <- x$active_factor

  excluded <- x$data$proc$excluded_features

  rarefied <- x$data$proc$rarefied

  filtering <- x$data$proc$filtering[!(names(x$data$proc$filtering) %in% c('filterlist','ftstats'))]
  filter_rank <- x$data$proc$filter_rank

  normalization <- x$data$proc$normalization

  # Now print out the info
  cat('\n~~~\n')
  cat(paste0(' | ',crayon::bold$underline('MicroVis Dataset:'),' ',saved))

  # Print number of included samples (and how many are specifically ignored, if any)
  cat(paste0(' |   > ',crayon::green$bold(n_samples),' samples '))
  if(length(n_ignored_samples)) cat(crayon::italic(paste0('(',crayon::yellow$bold(n_ignored_samples),
                                                          ' ignored samples)\n')))

  # Print number of included features and what type (taxa or pathways), as well as total number of features
  cat(paste0(' |   > ',crayon::green$bold(n_features), ' ',crayon::green(feature_type),
             ' out of ',crayon::cyan$bold(totn_features),'\n'))

  # Print all the independent variables (factors)
  cat(paste0(' |   > ',length(factors),' independent variable(s):\n'))
  for(f in factors) {
    if(f$name==active_factor) bullet <- '-> '
    else bullet <- ' * '
    cat(paste0(' |      ',bullet,f$name,' with ',crayon::green$bold(length(f$subset)),
               ' out of ',crayon::cyan$bold(length(f$groups)),' groups\n'))
  }

  # Print excluded features and number of ASVs in each
  if(length(excluded)) {
    cat(paste0(' |   > Excluded Features:\n'))
    for(ft in names(excluded)) cat(paste0(' |       * ',crayon::yellow(ft),
                                          ' (',length(excluded[[ft]]),' ASVs)\n'))
  }

  # Print rarefaction minimum richness
  if(length(rarefied)) {
    cat(paste0(' |   > Rarefied to ',crayon::yellow$bold(rarefied),' reads\n'))
  }

  # Print the normalization parameters
  if(length(normalization)) {
    cat(paste0(' |   > Normalization:\n'))
    for(norm in names(normalization)) {
      if(norm=='sample_scale') {
        cat(paste0(' |       * Samples ',crayon::yellow('scaled by',
                                                        normalization[[norm]][[1]]),'\n'))
      }
      if(norm=='feature_scale') {
        cat(paste0(' |       * Features ',crayon::yellow('scaled by',
                                                         normalization[[norm]][[1]]),'\n'))
      }
      if(norm=='transformation') {
        cat(paste0(' |       * Data ',crayon::yellow('transformed by',
                                                     normalization[[norm]][[1]]),'\n'))
      }
    }
  }

  # Print the filtering parameters
  if(length(filtering[!(names(filtering) %in% 'filterlist')])) {
    cat(paste0(' |   > Filtering (by ',filter_rank,'): \n'))
    for(filt in names(filtering)) {
      if(filt=='prevalence_proportion') {
        cat(paste0(' |       * Prevalence percent threshold: ',
                   crayon::yellow(expression("\u2265"),paste0(filtering[[filt]][[1]]*100,'%')),
                   ' of samples\n'))
      }
      if(filt=='min_prevalence') {
        cat(paste0(' |       * Low prevalence threshold: ',
                   crayon::yellow(expression("\u2265"),filtering[[filt]][[1]]),'\n'))
      }
      if(filt=='top_prevalence') {
        cat(paste0(' |       * Restricted to ',
                   crayon::yellow('top',filtering[[filt]]),' features by prevalence\n'))
      }
      if(filt=='min_relabun') {
        cat(paste0(' |       * Low relative abundance threshold: ',
                   crayon::yellow(expression("\u2265"),filtering[[filt]][[1]]),'\n'))
      }
      if(filt=='top_relabun') {
        cat(paste0(' |       * Restricted to ',
                   crayon::yellow('top',filtering[[filt]]),' features by relative abundance\n'))
      }
      if(filt=='min_totabun') {
        cat(paste0(' |       * Low total abundance threshold: ',
                   crayon::yellow(expression("\u2265"),filtering[[filt]][[1]]),'\n'))
      }
      if(filt=='top_totabun') {
        cat(paste0(' |       * Restricted to ',
                   crayon::yellow('top',filtering[[filt]]),' features by total abundance\n'))
      }
      if(filt=='low_var_percentile') {
        cat(paste0(' |       * Low variance percentile threshold: ',
                   crayon::yellow(expression("\u2265"),filtering[[filt]][[1]]),'\n'))
      }
      if(filt=='top_var') {
        cat(paste0(' |       * Restricted to ',
                   crayon::yellow('top',filtering[[filt]]),' features by standard deviation\n'))
      }
      if(filt=='low_abun') {
        cat(paste0(' |       * Abundance of at least ',
                   crayon::yellow(filtering[[filt]][[1]]),' in at least ',
                   crayon::yellow(filtering[[filt]][[2]]),' percent of samples\n'))
      }
      if(filt=='NAfilter') {
        cat(paste0(' |       * Removed taxa with unassigned ',
                   crayon::yellow(paste0(filtering[[filt]][[1]],collapse = ', ')),'\n'))
      }
    }
  }

  cat('~~~\n\n')
}

#' Print function for mvmerged (merged MicroVis dataset) objects
#'
#' @param x Merged MicroVis dataset (mvmerged object)
#' @param ... Arguments passed to print() function
#'
#' @return Print-friendly printout of dataset details
#' @export
#'
print.mvmerged <- function(x, ...) {
  # Get relevant dataset info
  if(is.null(x$name)) saved <- crayon::red$bold$italic('(Not saved)\n')
  else if(is.na(x$name)) saved <- crayon::italic('(Saved)\n')
  else saved <- crayon::italic('(Saved to',paste0('"',x$name,'")'),'\n')

  rank <- x$data$proc$active_rank

  if(x$features=='taxa') feature_type <- rank
  else feature_type <- x$features

  n_samples <- nrow(x$data$proc[[rank]])
  n_ignored_samples <- length(x$data$proc$ignored_samples)

  tot_nfeatures <- ncol(x$data$proc[[rank]])
  if('Other' %in% colnames(x$data$proc[[rank]])) tot_nfeatures <- tot_nfeatures - 1
  if('Unknown' %in% colnames(x$data$proc[[rank]])) tot_nfeatures <- tot_nfeatures - 1

  ds1_nfts <- length(x$data$features[[1]])
  ds2_nfts <- length(x$data$features[[2]])

  factors <- x$factors
  active_factor <- x$active_factor

  # Now print out the info
  cat('\n~~~\n')
  cat(paste0(' | ',crayon::bold$underline('Merged MicroVis Dataset:'),' ',saved))

  # Print number of included samples (and how many are specifically ignored, if any)
  cat(paste0(' |   > ',crayon::green$bold(n_samples),' samples '))
  if(length(n_ignored_samples)) cat(crayon::italic(paste0('(',crayon::yellow$bold(n_ignored_samples),
                                                          ' ignored samples)\n')))

  # Print number of total included features as well as the breakdown from each dataset
  cat(paste0(' |   > ',crayon::green$bold(tot_nfeatures), ' features from ',
             paste(crayon::italic(names(x$data$features)),collapse = ' and '),'\n'))
  cat(paste0(' |       * ',crayon::green$bold(ds1_nfts),' from ',
             crayon::italic(names(x$data$features)[[1]]),'\n'))
  cat(paste0(' |       * ',crayon::green$bold(ds2_nfts),' from ',
             crayon::italic(names(x$data$features)[[2]]),'\n'))

  # Print all the independent variables (factors)
  cat(paste0(' |   > ',length(factors),' independent variable(s):\n'))
  for(f in factors) {
    if(f$name==active_factor) bullet <- '-> '
    else bullet <- ' * '
    cat(paste0(' |      ',bullet,f$name,' with ',crayon::green$bold(length(f$subset)),
               ' out of ',crayon::cyan$bold(length(f$groups)),' groups\n'))
  }
  cat('~~~\n\n')
}
