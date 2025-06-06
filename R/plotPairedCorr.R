#' Correlate abundance of features between paired samples
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param data.paired (Alternative to dataset) Provide a dataframe with data already paired according to
#'     "ids" and "compare"
#' @param ids Column(s) that uniquely identify each pair of samples
#' @param compare Column with the 2 groups to correlate between. Each pair of
#'     samples has a sample in each of the two groups in this column. If there
#'     are more than 2 groups in the column, user will be asked to select 2 to
#'     compare
#' @param groups (Needed if specifying data.paired) Groups being compared.
#' @param fts (Optional) Specify vector of features to correlate between the two sets
#' @param invert (Optional) Switch the x and y axes?
#' @param corr_type (Optional) Specify method of calculating correlation. Either Pearson or Spearman. Defaults to Pearson.
#' @param rank (Optional) Specify taxonomic rank at which to perform the paired correlation. Defaults
#'     to the second highest rank in the dataset
#' @param alpha (Optional) Specify significance threshold. Defaults to 0.05
#' @param padj (Optional) Adjust the p-values for multiple testing? Defaults to TRUE
#' @param rthresh (Optional) Specify R value cutoff for coloring the figures of each figure.
#'     Defaults to 0.7
#' @param showstats (Optional) Whether or not to show R and p-value. Defaults to TRUE
#'
#' @return Scatter plot with trendline
#' @export
#'
plotPairedCor <- function(dataset=NULL,
                          data.paired=NULL, groups=NULL,
                          ids, compare,
                          fts=NULL, rank=NULL,
                          invert=F,
                          corr_type='pearson',
                          rthresh=0.7, alpha=0.05, padj=T,
                          showstats=T) {
  # if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)
  #
  # mdcols <- colnames(dataset$metadata)
  #
  # # If no valid rank was specified, default to the second highest rank
  # rank <- rank[rank %in% getRanks(dataset)]
  # if(!length(rank)) rank <- getRanks(dataset)[2]
  #
  # data <- mvmelt(dataset, rank = rank)
  # data$Other <- NULL
  # data$Unknown <- NULL
  # paired_samples <- getSamples(dataset, id_cols = ids, complete = compare)$sample
  #
  # data.paired <- data[data$sample %in% paired_samples,]
  # cat('\n',nrow(data.paired),'out of',nrow(data),'samples form complete pairs\n')
  #
  # compare <- compare[compare %in% colnames(data.paired)]
  # groups <- levels(data.paired[[compare]])
  # while(length(groups)!=2) {
  #   groups <- select.list(groups, multiple=T,
  #                         title='Please select two groups to correlate')
  # }
  #
  # data.paired <- data.paired[data.paired[[compare]] %in% groups,]
  # ids <- ids[ids %in% colnames(data.paired)]
  #
  # fts <- fts[fts %in% colnames(data.paired)]
  # if(!length(fts)) fts <- colnames(data.paired)[!(colnames(data.paired) %in% c(ids, compare, mdcols))
  #                                               & unlist(lapply(data.paired, is.numeric))]
  #
  # corr_type <- match.arg(corr_type, c('pearson','spearman'))

  results <- pairedCor(dataset=dataset,
                       data.paired=data.paired, groups=groups,
                       ids=ids, compare=compare,
                       fts=fts, rank=rank,
                       corr_type=corr_type,
                       rthresh=rthresh, alpha=alpha, padj=padj,
                       forPlot=T)

  stats <- results$stats
  data.paired <- results$data.paired
  groups <- results$groups
  ids <- results$ids
  compare <- results$compare
  fts <- results$fts
  rank <- results$rank
  corr_type <- results$corr_type

  if(invert) groups <- rev(groups)

  p <- list()
  for(ft in fts[stats$`p(adj)`<=1]) {
    pivoted <- data.paired[c(ids,compare,ft)] %>% tidyr::pivot_wider(id_cols=ids,
                                                                     names_from=compare,
                                                                     values_from=ft)
    if(stats$R[stats[[rank]]==ft] >= rthresh &
       stats$`p(adj)`[stats[[rank]]==ft] <= alpha) {
      ptemp <- ggpubr::ggscatter(pivoted, groups[1], groups[2], color='#00c700', size=4,
                                 add='reg.line', add.params=list(size=1.2), conf.int=T,
                                 title=ft)
    } else {
      ptemp <- ggpubr::ggscatter(pivoted, groups[1], groups[2], color='red', size=4, shape=18,
                                 add='reg.line', add.params=list(size=1.2), conf.int=T,
                                 title=ft)
    }

    if(showstats) ptemp <- ptemp + ggpubr::stat_cor(size=7, show.legend=F,
                                                    cor.coef.name = ifelse(corr_type=='pearson','R','rho'))

    p[[ft]] <- ptemp+
        theme(plot.title = element_text(size=30, hjust=0.5),
              axis.title = element_text(size=25),
              axis.text = element_text(size=20),
              legend.title = element_blank(),
              legend.key.size = unit(1,'cm'),
              legend.text = element_text(size=20))

    # Show the figure only if asking user if they want to save it, otherwise all figures will be shown at the end anyway
    if(get('autosave', envir = mvEnv) | get('offerSave', envir = mvEnv)) show(p[[ft]])

    if(!exists('save_one_all',inherits = F)) save_one_all <- NULL
    save_one_all <- multisave(save_one_all)

    if(save_one_all %in% c('Yes','Yes to all figures')) saveFig <- T
    else saveFig <- F

    if(saveFig) {
      save_directory <- saveResults(dataset,
                                    foldername='Paired Correlations',
                                    filename=ft,
                                    forcesave=T,
                                    verbose=F,
                                    width=8,height=6)
    }
  }

  assign('active_dataset',dataset,envir = mvEnv)
  if(!is.null(dataset$name)) assign(dataset$name,dataset,1)

  if(exists('save_directory')) {
    dir.create(file.path(save_directory,'Statistics'),showWarnings = FALSE)

    write.csv(stats,
              file=file.path(save_directory,'Statistics','Statistics.csv'),
              row.names = F)

    cat('Figures and any associated statistics saved to:\n ',save_directory)
  }

  cat('\n')

  return(p)
}
