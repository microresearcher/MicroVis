#' Switch to a Different Factor
#'
#' @param dataset Microvis dataset (mvdata)
#' @param factor_name (Optional) Name of factor to switch to
#' @param temp If set to TRUE, this will not update the active dataset.
#'
#' @return Dataset with new primary factor
#' @export
#'
changeFactor <- function(dataset=NULL, factor_name=NULL, temp=F) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- dataset$name
  }

  factors <- names(dataset$factors)
  active_factor <- dataset$active_factor

  if(!length(factors)) {return(message('\nERROR: No factors to choose from'))}

  if(is.null(factor_name)) {
    factor_name <- select.list(factors,'Which factor would you like to analyze by?',graphics=TRUE)
  } else if(!(tolower(factor_name) %in% tolower(factors))) {
    factor_name <- select.list(factors,'Which factor would you like to analyze by?',graphics=TRUE)
  }

  new_dataset <- dataset
  new_dataset$active_factor <- factor_name

  if(!validFactor(new_dataset)) {
    return(dataset)
  }

  if(!is.null(new_dataset$name)) assign(new_dataset$name,new_dataset,1)
  return(activate(new_dataset))
}

#' @title Factor Validation
#'
#' @description Function for checking factors to make sure they are valid for
#'     analysis. If the currently active factor only has 1 group, a new active
#'     factor will be chosen.
#'
#' @param dataset MicroVis dataset (mvdata object)
#'
#' @return Dataset with a new active factor.
#' @export
#'
validFactor <- function(dataset) {
  active_factor <- dataset$active_factor
  factors <- dataset$factors

  cat('\nCurrent active factor is',active_factor,'\n')

  valid <- T
  if(length(factors[[active_factor]]$subset) < 2) {
    message('\nWARNING: Fewer than 2 groups remaining in ', active_factor, ' comparative analyses will not work')
    valid <- F
  }

  grp_sizes <- countSamples(dataset,getSizes = T)
  if(sum(grp_sizes[[active_factor]]$Size>2)<2) {
    message('\nERROR: This factor has fewer than 2 groups remaining with at least 3 samples each\n')
    valid <- F
  }

  # If, after trying to reset the active factor, the number of subsetted groups is still < 1
  #   then warn the user
  # if(length(factors[[active_factor]]$subset) < 2) {
  #   message('\nWARNING: Fewer than 2 groups remaining in ', active_factor, ' comparative analyses will not work')
  # }

  return(valid)
}
