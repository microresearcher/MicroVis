#' DESeq2 Analysis
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor along which to perform analysis. Defaults to active factor
#' @param rank Rank at which to select features to analyze. Defaults to active rank
#' @param compareAll Perform pairwise comparisons between all groups or just
#'     between a reference group and all other groups? Defaults to FALSE
#'     (compares a reference group to all other groups)
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
                    compareAll=F,
                    dataset_name=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  # if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  # else dataset_name <- dataset$name

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- orig_rank <- dataset$data$proc$active_rank

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor
  factor <- setFVar(dataset,factor_name = factor)

  if(!compareAll) dds <- DESeq2::DESeq(makeDeseq(dataset,rank = rank))
  else dds <- DESeq2::DESeq(makeDeseq(dataset,factor$subset[1]))

  result_names <- DESeq2::resultsNames(dds)
  result_names <- result_names[2:length(result_names)]

  # First, extract the comparisons between contrast groups and the reference group
  for(contrast in result_names) {
    res <- DESeq2::results(dds,name=contrast)
    comparison <- strsplit(gsub(paste0('.*',factor$name_text,' '),'',res@elementMetadata$description[2]),
                           split = ' vs ')[[1]]

    base_grp <- factor$subset[gsub('_',' ',
                                   gsub('[- ;,/]','\\.',factor$subset)) %in% comparison[2]]
    comp_grp <- factor$subset[gsub('_',' ',
                                   gsub('[- ;,/]','\\.',factor$subset)) %in% comparison[1]]

    cat(paste0('\nAnalyzing "',comp_grp,'" vs "',base_grp,'"'))

    temp <- cbind(Reference=rep(base_grp,nrow(res)),
                  Contrast=rep(comp_grp,nrow(res)),
                  data.frame(Feature=rownames(res)),
                  res@listData)

    if(!exists('deseq_res',inherits = F)) deseq_res <- temp
    else deseq_res <- rbind(deseq_res,temp)
  }

  if(compareAll) {
    # If making all comparisons, extract the comparisons between the contrast groups
    contrast_levels <- gsub('_vs_.*','',gsub(paste0(factor$name,'_'),'',result_names))

    for(i in 1:(length(contrast_levels)-1)) {
      for(j in (i+1):length(contrast_levels)) {
        base_grp <- contrast_levels[i]
        comp_grp <- contrast_levels[j]

        cat('\nAnalyzing ',comp_grp,' vs ',base_grp)

        res <- DESeq2::results(dds,contrast = c(factor$name,comp_grp,base_grp))

        temp <- cbind(Reference=rep(base_grp,nrow(res)),
                      Contrast=rep(comp_grp,nrow(res)),
                      data.frame(Feature=rownames(res)),
                      res@listData)

        deseq_res <- rbind(deseq_res,temp)
      }
    }
  }

  cat('\n')

  dataset$stats[[factor$name]]$deseq[[rank]] <- deseq_res

  assign('active_dataset',dataset,envir = mvEnv)
  if(!is.null(dataset$name)) assign(dataset$name,dataset,1)

  return(dataset)
}
