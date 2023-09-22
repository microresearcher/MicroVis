#' Check that each group has at least a certain number of samples
#'
#' @param dataset MicroVis dataset. Default is the active dataset
#' @param factor Factor along which to group the samples
#' @param stratifiers Stratifiers along which to stratify the samples
#' @param min_n Minumum number of samples expected in each group
#' @param verbose If set to TRUE, prints out processing text
#'
#' @return TRUE or FALSE indicating whether any groups in the specified factor
#'    have fewer than the specified number of samples
#' @export
#'
checkGroups <- function(dataset=NULL,factor=NULL,stratifiers=NULL,min_n=3,verbose=T) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  stratifiers <- stratifiers[stratifiers %in% names(dataset$factors) & !(stratifiers %in% factor)]

  # Count the number of samples in each group
  grp_sizes <- countSamples(dataset,factor,stratifiers,getSizes=T,min_n=min_n,verbose=verbose)
  anysmallgrps <- F
  if(is.null(grp_sizes$stratified)) {
    for(i in 1:nrow(grp_sizes[[factor]])) if(grp_sizes[[factor]][i,'Size'] < min_n) {
      message(' ',grp_sizes[[factor]][i,factor],' has fewer than ',min_n,' samples')
      anysmallgrps <- T
    }
  } else {
    for(i in 1:nrow(grp_sizes$stratified[[factor]])) {
      if(grp_sizes$stratified[[factor]][i,'Size'] < min_n) {
        message(' ',as.character(interaction(grp_sizes$stratified[[factor]][i,c(factor,stratifiers)],
                                         sep=' and ')),
                ' has fewer than ',min_n,' samples')
        anysmallgrps <- T
      }
    }
  }
  if(verbose) cat('\n')

  return(anysmallgrps)
}
