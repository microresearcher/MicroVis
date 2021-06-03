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
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else if(is.null(dataset_name)) {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- setFVar(dataset)

  phylodata <- makePS(dataset)
  mm <- lefse(phylodata,
              normalization=1e6,
              class=factor$name,
              bootstrap_n=1000,
              correct='2')

  mm@marker_table$p_orig <- mm@marker_table$p_value
  mm@marker_table$p_value <- p.adjust(mm@marker_table$p_value,method = 'BH')

  dataset$stats[[factor$name]]$lefse <- mm

  if(dataset_name=='active_dataset') assign(dataset_name, dataset,envir = mvEnv)
  else assign(dataset_name,dataset,1)

  return(dataset)
}
