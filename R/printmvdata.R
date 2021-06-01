#' Print function for mvdata (MicroVis dataset) objects
#'
#' @param x MicroVis x (mvdata object)
#' @param ... Arguments passed to print() function
#'
#' @return Print-friendly printout of x details
#' @export
#'
print.mvdata <- function(x, ...) {
  # Get relevant dataset info
  rank <- x$data$proc$active_rank
  if(x$features=='taxa') feature_type <- rank
  else feature_type <- x$features
  n_samples <- nrow(x$data$proc[[rank]])
  n_ignored_samples <- length(x$data$proc$ignored_samples)

  n_features <- ncol(x$data$proc[[rank]])
  if('Other' %in% colnames(x$data$proc[[rank]])) n_features <- n_features - 1
  if(x$features=='taxa') totn_features <- length(unique(ASVtoTaxa(x$data,
                                                                        colnames(x$data$orig),
                                                                        taxa_rank=rank)))
  else totn_features <- ncol(x$data$orig)

  factors <- x$factors
  active_factor <- x$active_factor

  rarefied <- x$data$proc$rarefied

  filtering <- x$data$proc$filtering[!(names(x$data$proc$filtering) %in% c('filterlist','ftstats'))]
  filter_rank <- x$data$proc$filtering$filter_rank

  normalization <- x$data$proc$normalization

  # Now print out the info
  cat('\n~~~\n')
  cat(paste0(' | MicroVis Dataset:\n'))

  # Print number of included samples (and how many are specifically ignored, if any)
  cat(paste0(' |   > ',n_samples,' samples '))
  if(length(n_ignored_samples)) cat(paste0('(',n_ignored_samples,' ignored samples)\n'))

  # Print number of included features and what type (taxa or pathways), as well as total number of features
  cat(paste0(' |   > ',n_features, ' ', feature_type,' out of ',totn_features,'\n'))

  # Print all the independent variables (factors)
  cat(paste0(' |   > ',length(factors),' independent variable(s):\n'))
  for(f in factors) {
    if(f$name==active_factor) bullet <- '-> '
    else bullet <- ' * '
    cat(paste0(' |      ',bullet,f$name,' with ', length(f$subset),' out of ',length(f$groups),' groups\n'))
  }

  # Print rarefaction minimum richness
  if(length(rarefied)) {
    cat(paste0(' |   > Rarefied to ',rarefied,' reads\n'))
  }

  # Print the normalization parameters
  if(length(normalization)) {
    cat(paste0(' |   > Normalization:\n'))
    for(norm in names(normalization)) {
      if(norm=='sample_scale') {
        cat(paste0(' |       * Samples scaled by ',normalization[[norm]][[1]],'\n'))
      }
      if(norm=='feature_scale') {
        cat(paste0(' |       * Features scaled by ',normalization[[norm]][[1]],'\n'))
      }
      if(norm=='transformation') {
        cat(paste0(' |       * Data transformed by ',normalization[[norm]][[1]],'\n'))
      }
    }
  }

  # Print the filtering parameters
  if(length(filtering)>1) {
    cat(paste0(' |   > Filtering (by ',filter_rank,'): \n'))
    for(filt in names(filtering)) {
      if(filt=='min_prevalence') {
        cat(paste0(' |       * Low prevalence threshold: ',expression("\u2265"),' ',filtering[[filt]][[1]],'\n'))
      }
      if(filt=='top_prevalence') {
        cat(paste0(' |       * Restricted to top ',filtering[[filt]],' features by prevalence\n'))
      }
      if(filt=='min_relabun') {
        cat(paste0(' |       * Low relative abundance threshold: ',expression("\u2265"),' ',filtering[[filt]][[1]],'\n'))
      }
      if(filt=='top_relabun') {
        cat(paste0(' |       * Restricted to top ',filtering[[filt]],' features by relative abundance\n'))
      }
      if(filt=='min_totabun') {
        cat(paste0(' |       * Low total abundance threshold: ',expression("\u2265"),' ',filtering[[filt]][[1]],'\n'))
      }
      if(filt=='top_totabun') {
        cat(paste0(' |       * Restricted to top ',filtering[[filt]],' features by total abundance\n'))
      }
      if(filt=='low_var_percentile') {
        cat(paste0(' |       * Low variance percentile threshold: ',expression("\u2265"),' ',filtering[[filt]][[1]],'\n'))
      }
      if(filt=='top_var') {
        cat(paste0(' |       * Restricted to top ',filtering[[filt]],' features by standard deviation\n'))
      }
      if(filt=='low_abun') {
        cat(paste0(' |       * Abundance of at least ', filtering[[filt]][[1]], ' in at least ',filtering[[filt]][[2]],' percent of samples\n'))
      }
      if(filt=='NAfilter') {
        cat(paste0(' |       * Removed taxa with unassigned ', paste0(filtering[[filt]][[1]],collapse = ', '),'\n'))
      }
    }
  }

  cat('~~~\n\n')
}
