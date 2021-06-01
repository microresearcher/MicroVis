#' @title Process Dataset for Analysis
#'
#' @description Called by a number of functions to process a dataset based on
#'     set parameters.
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp If TRUE, will not set the processed dataset as the active dataset.
#' @param silent If TRUE, none of the processing steps will print out.
#'
#' @return Processed MicroVis dataset (mvdata object)
#' @export
#'
processDataset <- function(dataset, temp=F, silent=F) {
  if(!get('detailed_taxa_names',envir = mvEnv) & dataset$features=='taxa') {
    dataset <- cleanUnkTaxa(dataset)
  }

  dataset <- runSampleFilter(dataset, silent=silent)

  if(!silent) p_beforeproc <- plotSampleDensity(dataset)

  dataset <- runRarefaction(dataset, silent=silent)
  dataset <- runNormalization(dataset, silent=silent)
  dataset <- runFeatureFilter(dataset, silent=silent)

  if(!silent) {
    p_afterproc <- plotSampleDensity(dataset)
    show(ggarrange(p_beforeproc,p_afterproc,ncol=1))
  }

  dataset$stats <- list()
  dataset$data$proc$filtering$ftstats <- NULL
  class(dataset) <- 'mvdata'

  if(!temp) {assign('active_dataset',dataset,envir = mvEnv)}
  return(dataset)
}

#' Clear Normalization and Feature Filtering
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param clearSampleFilter Not currently functional
#' @param temp If TRUE, will not set the unprocessed dataset as the active dataset.
#' @param silent If TRUE, none of the processing steps will print out.
#'
#' @return Processed MicroVis dataset (mvdata)
#' @export
#'
clearProcessing <- function(dataset, clearSampleFilter=F, temp=F, silent=F) {
  dataset <- clearNormalization(dataset, silent=silent, temp = T)
  dataset <- clearFeatureFilt(dataset, silent=silent, temp = T)

  dataset$stats <- list()
  dataset$data$proc$filtering$ftstats <- NULL
  class(dataset) <- 'mvdata'

  if(!temp) {assign('active_dataset',dataset,envir = mvEnv)}
  return(dataset)
}


#' Activate Dataset
#'
#' @param dataset MicroVis dataset (mvdata object)
#'
#' @return MicroVis dataset (mvdata object)
#' @export
#'
activate <- function(dataset) {
  dataset_name <- deparse(substitute(dataset))

  cat(paste0('\n\n  <|> Active Dataset: "',dataset_name,'" <|>\n\n'))
  print(dataset)
  assign('active_dataset',dataset,envir = mvEnv)
}
