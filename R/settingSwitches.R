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

#' @title Switch Method for Replacing Zeros
#'
#' @description Cycles between two options for dealing with zero counts:
#'     - replace: Replacing zero counts with a value between 1/100 and 1/10 of the minimum count value
#'     - impute: Use zCompositions package to impute a distribution of values to replace zeros
#'
#' @param method (Optional) The method of zero-replacement to switch to
#' @param div (Optional) The divisor by which to divide the lowest count value by
#' @param rdist (Optional) How to replace zeros, either using a uniform distribution
#'     or a point value
#'     to determine the maximum value to replace zeros with for the "replace" method
#' @return Does not return anything
#' @export
#'
zerostoggle <- function(method=c('replace','impute'),
                        div=NULL,
                        rdist=c('runif','point')) {
  if(length(method)>1) {
    method <- get('zeroReplaceMethod',envir = mvEnv)$method
    if(is.null(div) & length(rdist)>1) {
      if(method=='replace') method <- 'impute'
      else if(method=='impute') method <- 'replace'
    }
  } else method <- match.arg(method)

  if(!is.null(div)) {
    if(!is.numeric(div) | div<=1) stop('"div" must be a number greater than 1')
  } else div <- get('zeroReplaceMethod',envir = mvEnv)$div

  if(length(rdist)==1) rdist <- match.arg(rdist)
  else rdist <- get('zeroReplaceMethod',envir = mvEnv)$rdist

  assign('zeroReplaceMethod',list(method=method,
                                  div=div,
                                  rdist=rdist),envir = mvEnv)

  cat('Zero replacement method during normalization steps:',method,'\n')
}

#' @title Automatically Name Results
#'
#' @description Cycles between naming results automatically or prompting user
#'     for a name
#'
#' @param setting (Optional) Whether to set naming of files to auto-naming by
#'     groups in analysis, auto-naming by the name of the dataset, or manual
#'     naming by the user. By default, this function cycles through these 3 settings
#'
#' @return Does not return anything
#' @export
#'
autonametoggle <- function(setting=c('groups','dataset','manual')) {
  if(length(setting)>1) setting <- get('autoNameResults',envir = mvEnv)
  else setting <- match.arg(setting)

  if(setting=='groups') setting <- 'dataset'
  else if(setting=='dataset') setting <- 'manual'
  else if(setting=='manual') setting <- 'groups'

  assign('autoNameResults',setting,envir = mvEnv)

  cat('\nAuto-naming of figures is set to',setting,'\n\n')
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
