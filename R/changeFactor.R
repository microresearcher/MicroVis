#' Switch to a Different Factor
#'
#' @param dataset Microvis dataset (mvdata).
#' @param factor_name (Optional) Name of factor to switch to.
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

  if(!length(factors)) return(message('\nERROR: No factors to choose from'))

  if(is.null(factor_name)) {
    factor_name <- select.list(factors,'Which factor would you like to analyze by?',graphics=TRUE)
  } else if(!(tolower(factor_name) %in% tolower(factors))) {
    factor_name <- select.list(factors,'Which factor would you like to analyze by?',graphics=TRUE)
  }

  new_dataset <- dataset

  if(validFactor(dataset, factor_name)) new_dataset$active_factor <- factor_name
  else stop('\nFactor has not been changed. See warning message below.\n')

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
#' @param factor Candidate factor to switch to.
#'
#' @return Dataset with a new active factor.
#' @export
#'
validFactor <- function(dataset, factor_name=NULL) {
  factors <- dataset$factors
  factor_name <- factor_name[factor_name %in% names(factors)]

  if(is.null(factor_name)) chosen_factor <- dataset$active_factor
  else chosen_factor <- factor_name

  cat('\nCurrent active factor is',chosen_factor,'\n')

  valid <- T
  if(length(factors[[chosen_factor]]$subset) < 2) {
    warning('\n Fewer than 2 groups remaining in "', chosen_factor, '". Comparative analyses will not work')
    valid <- F
  }

  grp_sizes <- countSamples(dataset,getSizes = T)
  if(sum(grp_sizes[[chosen_factor]]$Size>2)<2) {
    warning('\n This factor has fewer than 2 groups remaining with at least 3 samples each\n')
    valid <- F
  }

  # If, after trying to reset the active factor, the number of subsetted groups is still < 1
  #   then warn the user
  # if(length(factors[[dataset$active_factor]]$subset) < 2) {
  #   message('\nWARNING: Fewer than 2 groups remaining in ', active_factor, ' comparative analyses will not work')
  # }

  return(valid)
}
