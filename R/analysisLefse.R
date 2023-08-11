#' LEfSe Analysis
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param dataset_name (Not Recommended) Name of the dataset to save results
#'     to. This should not need to be used by users since the function can
#'     determine the name of the dataset directly passed to it, but not when
#'     it is called within another function.
#'
#' @return MicroVis dataset containing results of LEfSe analysis
#' @export
#'
mvlefse <- function(dataset=NULL, dataset_name=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  factor <- setFVar(dataset)

  phylodata <- makePS(dataset)
  mm <- microbiomeMarker::run_lefse(phylodata,
                                    group=factor$name,
                                    norm='CPM',
                                    bootstrap_n=1000)

  mm@marker_table$BHadjust <- p.adjust(mm@marker_table$pvalue,method = 'BH')

  dataset$stats[[factor$name]]$lefse <- mm

  assign('active_dataset', dataset,envir = mvEnv)
  if(dataset_name!='active_dataset') assign(dataset_name,dataset,1)

  return(dataset)
}
