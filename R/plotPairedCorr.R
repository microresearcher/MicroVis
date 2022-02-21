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
#' @param rank Taxonomic rank at which to perform the paired correlation. Defaults
#'     to the second highest rank in the dataset
#' @param alpha Significance threshold. Defaults to 0.05
#' @param padj Adjust the p-values for multiple testing? Defaults to TRUE
#' @param rthresh R value cutoff for coloring the figures of each figure.
#'     Defaults to 0.7
#' @param showstats Whether or not to show R and p-value. Defaults to TRUE
#'
#' @return Scatter plot with trendline
#' @export
#'
plotPairedCor <- function(dataset=NULL, ids, compare, fts=NULL, rank=NULL,
                          invert=F,
                          rthresh=0.7, alpha=0.05, padj=T,
                          showstats=T) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  # If no valid rank was specified, default to the second highest rank
  rank <- rank[rank %in% getRanks(dataset)]
  if(!length(rank)) rank <- getRanks(dataset)[2]

  data <- mvmelt(dataset, rank = rank)
  data$Other <- NULL
  data$Unknown <- NULL
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

  p <- list()
  for(ft in fts[cor_pvals<=1]) {
    pivoted <- data.paired[c(ids,compare,ft)] %>% pivot_wider(id_cols=ids,
                                                              names_from=compare,
                                                              values_from=ft)
    if(cor_rvals[match(ft,fts)] >= rthresh & cor_pvals[match(ft,fts)] <= alpha) {
      ptemp <- ggscatter(pivoted, groups[1], groups[2], color='#00c700', size=4,
                         add='reg.line', add.params=list(size=1.2), conf.int=T,
                         title=ft)
    } else {
      ptemp <- ggscatter(pivoted, groups[1], groups[2], color='red', size=4, shape=18,
                         add='reg.line', add.params=list(size=1.2), conf.int=T,
                         title=ft)
    }

    if(showstats) ptemp <- ptemp + stat_cor(size=7,show.legend=F)

    p[[ft]] <- ptemp+
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
      save_directory <- saveResults(dataset$results_path,
                                    foldername='Paired Correlations',
                                    filename=ft,
                                    forcesave=T,
                                    verbose=F,
                                    width=8,height=6)
    }
  }

  if(exists('save_directory')) {
    dir.create(file.path(save_directory,'Statistics'),showWarnings = FALSE)

    write.csv(cor_pvals,
              file=file.path(save_directory,'Statistics','Statistics.csv'))

    cat('Figures and any associated statistics saved to:\n ',save_directory)
  }
}
