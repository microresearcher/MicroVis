#' Get mean and SE of features
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param factor Factor along which to conduct beta-diversity analysis
#' @param ranks Ranks at which to summarize features
#' @param allRanks If set to TRUE, will summarize features at all ranks. Defaults
#'     to FALSE
#' @param meansOnly Only return the means and not SE
#' @param pivoted If set to TRUE, the returned table will have each group of the
#'     factor in its own column
#'
#' @return Table with the means +/- SE for features in each group of the dataset
#' @export
#'
summariseTaxa <- function(dataset=NULL, factor=NULL,
                          ranks=NULL, allRanks=F,
                          meansOnly=F, pivoted=T) {

  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  if(allRanks) ranks <- getRanks(dataset)
  else {
    ranks <- ranks[ranks %in% getRanks(dataset)]
    if(is.null(ranks)) ranks <- dataset$data$proc$active_rank
  }

  summary_tab <- data.frame(Rank=0,Feature=0,Group=0,Mean=0,SE=0)[0,]
  for(rank in ranks) {
    fts <- getFeatures(dataset, ranks=rank)

    abun <- mvmelt(dataset, rank=rank, features=fts)

    abun_bygroup <- split(abun,abun[[factor]])
    for(grp in abun_bygroup) {
      grp_name <- as.character(grp[[factor]][1])
      grptab <- data.frame(lapply(grp[fts],function(x) as.numeric(as.character(x))))
      colnames(grptab) <- fts
      means <- colMeans(grptab[fts])
      SE <- unlist(lapply(grptab,function(x) sd(x)/sqrt(length(x))))
      grpsumtab <- data.frame(rep.int(rank,length(fts)),fts,rep.int(grp_name,length(fts)),means,SE)
      row.names(grpsumtab) <- c()
      colnames(grpsumtab) <- c('Rank','Feature','Group','Mean','SE')
      summary_tab <- rbind(summary_tab, grpsumtab)
    }
  }

  if(meansOnly) {
    summary_tab$SE <- NULL
    summary_tab <- pivot_wider(summary_tab,
                               names_from = 'Group',
                               values_from = 'Mean')
  } else if(pivoted) {
    summary_tab <- pivot_wider(summary_tab,
                               names_from = 'Group',
                               values_from = c('Mean','SE'),
                               names_glue='{Group} {.value}')
  }

  return(summary_tab)
}
