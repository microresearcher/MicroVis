#' Capitalize a Single-Word String
#'
#' @param string Single-word string.
#'
#' @return String with first letter capitalized and all other letters in lowercase.
#' @export
#'
capitalize <- function(string) {
  capitalized <- paste0(toupper(substr(string,1,1)),tolower(substr(string,2,nchar(string))))
  return(capitalized)
}
