#' Plot Heat Tree
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param test Type of statistical test to use. One of either "univar" or "deseq".
#'     Defaults to univariate
#' @param ftlist (Optional) List of specific features to show
#' @param showAll Show all features? Defaults to FALSE
#' @param color_limit Any log2FC values with absolute value above this will be
#'     treated as if they are at this value for coloring
#' @param alpha Significance threshold. Defaults to 0.05
#' @param layout Type of tree layout
#'
#' @return MicroVis dataset
#' @export
#'
plotHeatTree <- function(dataset=NULL,
                         test='univar',ftlist=NULL,showAll=F,
                         color_limit=2,
                         alpha=0.05,
                         initial_layout='fruchterman-reingold',
                         layout='reingold-tilford') {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  rank <- getLowestRank(dataset)
  ftlist <- ftlist[ftlist %in% getFeatures(dataset, ranks=rank)]

  if(test=='univar') {
    sigfts <- listsigs(dataset, allRanks=T,
                       dataset_name=dataset_name,
                       silent=T, alpha=alpha)
  } else if(is.null(ftlist)) ftlist <- getFeatures(dataset, ranks=rank)

  fc_tab <- foldChange(dataset, allRanks=T)
  fc_tab$FC[is.na(fc_tab$FC) | is.infinite(fc_tab$FC)] <- 0
  # Change non-real values to 0
  fc_tab$Log2FC[is.na(fc_tab$Log2FC)] <- 0
  fc_tab$Log2FC[is.infinite(fc_tab$Log2FC)] <- 0
  # Change insignificant values to 0
  if(exists('sigfts',inherits=F)) fc_tab$Log2FC[!(fc_tab$Feature %in% sigfts)] <- 0

  mvtaxmap <- makeTaxMap(dataset, ftlist=ftlist)

  color_limit <- abs(color_limit)
  #TODO: Unfortunately, it seems node_color must be a number... so we can only
  #       compare 2 groups at a time
  set.seed(999)
  for(grp in unique(fc_tab$Comparison)) {
    ht <- heat_tree(mvtaxmap,
                    node_label = mvtaxmap$taxon_names(),
                    node_size = mvtaxmap$n_obs(),
                    node_size_axis_label = 'Prevalence',
                    node_color = fc_tab$Log2FC[fc_tab$Comparison==grp],
                    node_color_interval = c(-1*color_limit,color_limit),
                    node_color_axis_label = paste0('Log2FC ',grp,'/',
                                                   fc_tab$Reference[[1]]),
                    node_color_range = c('cyan', 'gray', 'red'),
                    layout = layout)
    print(ht)

    if(!exists('save_one_all',inherits = F)) save_one_all <- NULL
    save_one_all <- multisave(save_one_all)

    if(save_one_all %in% c('Yes','Yes to all figures')) saveFig <- T
    else saveFig <- F

    if(saveFig) {
      save_directory <- saveResults(dataset$results_path,
                                    foldername = paste0('Heattree_',test),
                                    filename = paste0('log2FC_',grp,'v',fc_tab$Reference[[1]]),
                                    width = 8, height = 6,
                                    forcesave = T,
                                    verbose = F)
    }
  }

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))
  return(dataset)
}
