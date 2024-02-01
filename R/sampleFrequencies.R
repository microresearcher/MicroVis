#' Frequency table of samples
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param characteristics A vector of 2 column names in the metadata with which to create the frequency table.
#'    First variable will be rows and second will be columsn.
#' @param unique (Optional) Column to use to select unique samples
#' @param subset (Optional) A named list where names are additional column names, and their values are values to subset the samples by.
#' @param raw Whether to create frequency table from unprocessed data so that any samples excluded due to low quality or other sample filtering are included.
#'    Defaults to False.
#'
#' @return Returns a frequency table of two @characteristics
#' @export
#'
sampleFreq <- function(dataset=NULL, characteristics, unique=NULL, subset=NULL, raw=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(!raw) data <- mvmelt(dataset)[1:ncol(dataset$metadata)]
  else data <- dataset$metadata

  characteristics <- characteristics[characteristics %in% colnames(data)]
  if(length(characteristics) != 2) stop('Must select 2 column names from the metadata')

  if(length(unique)) {
    unique <- unique[unique %in% colnames(data)]
    unique <- unique[!(unique %in% characteristics)]
  }

  if(length(subset)) {
    subset <- subset[subset %in% colnames(data)]
    subset <- subset[!(subset %in% characteristics)]
  }

  tab.long <- data[c(characteristics, unique)]
  tab.long <- tab.long[!duplicated(tab.long),][characteristics]
  tab.long <- cbind(tab.long, 'count' = rep(1, nrow(tab.long)))

  freq_tab <- table(tab.long)

  return(freq_tab)
}
