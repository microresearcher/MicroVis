#' Correlate abundance of features between paired samples
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param ids Column(s) that uniquely identify each pair of samples
#' @param compare Column with the 2 groups to correlate between. Each pair of
#'     samples has a sample in each of the two groups in this column. If there
#'     are more than 2 groups in the column, user will be asked to select 2 to
#'     compare
#' @param fts Vector of features to correlate between the two sets
#' @param invert Switch the x and y axes?
#'
#' @return Scatter plot with trendline
#' @export
#'
plotPairedCor <- function(dataset=NULL, ids, compare, fts=NULL, invert=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  data <- mvmelt(dataset)
  data$Other <- NULL
  paired_samples <- getSamples(dataset, id_cols = ids, complete = compare)$sample

  data.paired <- data[data$sample %in% paired_samples,]
  cat('\n',nrow(data.paired),'out of',nrow(data),'samples form complete pairs\n')

  fts <- fts[fts %in% colnames(data.paired)]
  if(!length(fts)) fts <- colnames(data.paired)[!(colnames(data.paired) %in% c(ids, compare))
                                         & unlist(lapply(data.paired, is.numeric))]

  compare <- compare[compare %in% colnames(data.paired)]
  groups <- levels(data.paired[[compare]])
  while(length(groups)!=2) {
    groups <- select.list(groups, multiple=T,
                          title='Please select two groups to correlate')
  }
  if(invert) groups <- rev(groups)
  data.paired <- data.paired[data.paired[[compare]] %in% groups,]
  ids <- ids[ids %in% colnames(data.paired)]

  cor_pvals <- c()
  for(ft in fts) {
    pivoted <- data.paired[c(ids,compare,ft)] %>% pivot_wider(id_cols=ids,
                                                              names_from=compare,
                                                              values_from=ft)
    cor_pvals[[ft]] <- rcorr(as.matrix(pivoted[groups]))$P[1,2]
  }
  cor_pvals <- p.adjust(cor_pvals, method = 'BH')
  cor_pvals[is.na(cor_pvals)] <- 1

  p <- list()
  for(ft in fts[cor_pvals<=0.05]) {
    pivoted <- data.paired[c(ids,compare,ft)] %>% pivot_wider(id_cols=ids,
                                                              names_from=compare,
                                                              values_from=ft)
    p[[ft]] <- ggscatter(pivoted, groups[1], groups[2], color='blue', size=4,
                         add='reg.line', add.params=list(size=1.2), conf.int=T,
                         title=ft)+
      stat_cor(size=7,show.legend=F)+
      theme(plot.title = element_text(size=25, hjust=0.5),
            axis.title = element_text(size=25),
            axis.text = element_text(size=20),
            legend.title = element_blank(),
            legend.key.size = unit(1,'cm'),
            legend.text = element_text(size=20))

    show(p[[ft]])

    if(!exists('save_one_all',inherits = F)) save_one_all <- NULL
    save_one_all <- multisave(save_one_all)

    if(save_one_all %in% c('Yes','Yes to all figures')) saveFig <- T
    else saveFig <- F

    if(saveFig) {
      save_directory <- saveResults(save_path,
                                    foldername='Paired Correlations',
                                    filename=ft,
                                    forcesave=T,
                                    verbose=F,
                                    width=8,height=6)
    }
  }
}
