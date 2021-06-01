#' @title Set Factor Variable
#'
#' @description This function just assigns a factor in the dataset to a variable
#'     to be used within a function. If no factor name is specified in the function
#'     call, it will choose the dataset's active factor.
#'
#' @param dataset Microvis dataset (mvdata object)
#' @param factor_name (Optional) Name of the desired factor
#'
#' @return Factor variable (a list object) with its name, groups, subsetted groups,
#'     and print-friendly name (spaces instead of underscores)
#' @export
#'
setFVar <- function(dataset, factor_name=NULL) {
  varname <- list()
  factors <- dataset$factors
  active_factor <- dataset$active_factor

  # If there are no factors in flist, the user-friendly variable is set to point to the "sample" column of metadata (which is renamed internally in processMDFile())
  if(!length(dataset$factors)) {
    varname$name <- 'sample'
  } else {
    if(is.null(factor_name)) {
      factor_name <- active_factor
    }
    if(tolower(factor_name) %in% tolower(names(factors))) {
      varname <- factors[tolower(names(factors)) %in% tolower(factor_name)][[factor_name]]
    } else {
      cat(paste('\nCould not find',factor_name,'among the independent variables being analyzed\n'))
    }
  }

  return(varname)
}
