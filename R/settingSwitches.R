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
