#' Correlate abundance of features between paired samples
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param ids Column(s) that uniquely identify each pair of samples
#' @param compare Column with the 2 groups to correlate between. Each pair of
#'     samples has a sample in each of the two groups in this column. If there
#'     are more than 2 groups in the column, user will be asked to select 2 to
#'     compare
#' @param fts Vector of features to correlate between the two sets
#' @param rank Taxonomic rank at which to perform the paired correlation. Defaults
#'     to the second highest rank in the dataset
#' @param alpha Significance threshold. Defaults to 0.05
#' @param padj Adjust the p-values for multiple testing? Defaults to TRUE
#' @param rthresh R value cutoff for coloring the figures of each figure.
#'     Defaults to 0.7
#' @param data.paired Provide a dataframe with data already paired according to
#'     "ids" and "compare"
#' @param groups If providing a data.paired dataframe, must provide groups as well
#' @return Scatter plot with trendline
#' @export
#'
pairedCor <- function(dataset=NULL, ids, compare, fts=NULL, rank=NULL,
                      rthresh=0.7, alpha=0.05, padj=T,
                      data.paired=NULL, groups=NULL) {
  if(is.null(data.paired)) {
    if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

    mdcols <- colnames(dataset$metadata)

    # If no valid rank was specified, default to the second highest rank
    rank <- rank[rank %in% getRanks(dataset)]
    if(!length(rank)) rank <- getRanks(dataset)[2]

    data <- mvmelt(dataset, rank = rank)
    data$Other <- NULL
    data$Unknown <- NULL
    paired_samples <- getSamples(dataset, id_cols = ids, complete = compare)$sample

    data.paired <- data[data$sample %in% paired_samples,]
    cat('\n',nrow(data.paired),'out of',nrow(data),'samples form complete pairs\n\n')

    compare <- compare[compare %in% colnames(data.paired)]
    groups <- levels(data.paired[[compare]])
    while(length(groups)!=2) {
      groups <- select.list(groups, multiple=T,
                            title='Please select two groups to correlate')
    }

    data.paired <- data.paired[data.paired[[compare]] %in% groups,]
    ids <- ids[ids %in% colnames(data.paired)]
  } else if(length(groups)!=2) {
    stop('Must specify 2 groups within ',compare,' to correlate')
  }
  else if(any(is.na(data.paired %>% pivot_wider(id_cols = ids,
                                                  names_from = compare,
                                                  values_from = 'sample')))) {
    stop('Data for each ',paste0(ids,collapse = '-'),' must have exactly 2 values for ',compare)
  }

  fts <- fts[fts %in% colnames(data.paired)]
  if(!length(fts)) fts <- colnames(data.paired)[!(colnames(data.paired) %in% c(ids, compare, mdcols))
                                                & unlist(lapply(data.paired, is.numeric))]

  cor_rvals <- c()
  cor_pvals <- c()
  for(ft in fts) {
    pivoted <- data.paired[c(ids,compare,ft)] %>% pivot_wider(id_cols=ids,
                                                              names_from=compare,
                                                              values_from=ft)
    corvals <- rcorr(as.matrix(pivoted[groups]))
    cor_rvals[[ft]] <- corvals$r[1,2]
    cor_pvals[[ft]] <- corvals$P[1,2]
  }
  if(padj) cor_pvals <- p.adjust(cor_pvals, method = 'BH')

  cor_rvals[is.na(cor_rvals)] <- 0
  cor_pvals[is.na(cor_pvals)] <- 1.1

  stats <- merge(data.frame(unlist(cor_rvals)),data.frame(unlist(cor_pvals)),by = 0)
  colnames(stats) <- c(rank, 'R', 'p(adj)')

  if(!is.null(dataset)) {
    dataset$stats$paired_cor[[paste0(groups,collapse = '_')]][[rank]] <- stats

    assign('active_dataset',dataset,envir = mvEnv)
    if(!is.null(dataset$name)) assign(dataset$name,dataset,1)
  }

  return(stats)
}
