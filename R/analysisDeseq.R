#' DESeq2 Analysis
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor along which to perform analysis. Defaults to active factor
#' @param rank Rank at which to select features to analyze. Defaults to active rank
#' @param dataset_name (Not recommended) Name of the dataset to save results
#'     to. This should not need to be used by users since the function can
#'     determine the name of the dataset directly passed to it, but not when
#'     it is called within another function.
#'
#' @return Deseq analysis results
#' @export
#'
mvdeseq <- function(dataset=NULL,
                    factor=NULL,
                    rank=NULL,
                    dataset_name=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else if(is.null(dataset_name)) {
    dataset_name <- deparse(substitute(dataset))
  }

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- setFVar(dataset)

  dds <- DESeq(makeDeseq(dataset))
  result_names <- resultsNames(dds)
  result_names <- result_names[2:length(result_names)]

  for(contrast in result_names) {
    res <- results(dds,name=contrast)
    comparison <- strsplit(gsub(paste0('.*',factor$name,' '),'',res@elementMetadata$description[2]),
                           split = ' vs ')[[1]]
    base_grp <- comparison[2]
    contrast <- comparison[1]

    temp <- cbind(Reference=rep(base_grp,nrow(res)),
                  Contrast=rep(contrast,nrow(res)),
                  data.frame(Feature=rownames(res)),
                  res@listData)
    if(!exists('deseq_res',inherits = F)) deseq_res <- temp
    else deseq_res <- rbind(deseq_res,temp)
  }

  dataset$stats[[factor$name]][[rank]]$deseq <- deseq_res
  if(dataset_name=='active_dataset') assign(dataset_name,dataset,envir = mvEnv)
  else assign(dataset_name,dataset,1)

  return(deseq_res)
}
