#' @title Switch Saving Mode
#'
#' @description Cycles between three options for saving results of plotting:
#'     - Offer to save after each plotting function is run
#'     - Automatically save without asking
#'     - Do not save anything and don't ask the user either
#'
#' @return Does not return anything
#' @export
#'
savetoggle <- function() {
  if(get('autosave',envir = mvEnv)) {
    assign('autosave',F,envir = mvEnv)
    assign('offerSave',T,envir = mvEnv)
    cat('\nYou will be prompted before saving each figure\n\n')
  } else if(get('offerSave',envir = mvEnv)) {
    assign('offerSave',F,envir = mvEnv)
    cat('\nNo figures will be saved\n\n')
  } else {
    assign('autosave',T,envir = mvEnv)
    cat('\nFigures will automatically be saved\n\n')
  }
}

#' @title Automatically Name Results
#'
#' @description Cycles between naming results automatically or prompting user
#'     for a name
#'
#' @return Does not return anything
#' @export
#'
autonametoggle <- function() {
  assign('autoNameResults',!get('autoNameResults',envir = mvEnv),envir = mvEnv)
  cat('\nAuto-naming of figures is set to',get('autoNameResults',envir = mvEnv),'\n\n')
}

#' Switch between Tukey/Games-Howell or T-test/Wilcox for post-hoc analysis
#'
#' @return Does not return anything
#' @export
#'
posthoctoggle <- function() {
  assign('tukey_games',!get('tukey_games',envir = mvEnv),envir = mvEnv)
  if(get('tukey_games',envir = mvEnv)) {
    cat('\nPost-hoc analysis will be performed using either Tukey (parametric) or Games-Howell (nonparametric)\n\n')
  } else {
    cat('\nPost-hoc analysis will be performed using either T-test (parametric) or Wilcox (nonparametric) and then corrected for multiple comparisons\n\n')
  }
}

#' Function for setting whether or not to run Fisher tests to prune filtering
#'
#' @param keep Keep features with significant Fisher test if they would have been
#'     filtered out otherwise? Default is TRUE
#'
#' @return NULL
#' @export
#'
keepSigFisher <- function(keep=T) {
  assign('keepSigFisher',keep,envir = mvEnv)
  if(keep) cat('\nSignificant features by fisher test between the groups will be kept')
  else cat('\nFisher test will not be run to retain features with significantly different proportions among groups')
}
