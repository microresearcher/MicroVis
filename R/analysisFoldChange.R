#' Calculate Fold Change
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param factor Factor along which to conduct beta-diversity analysis
#' @param ref_group Reference group for fold change computation. This will be the
#'     denominator. User will be asked if none is specified
#' @param comparison_groups Comparison group for fold change computation. This
#'     will be the numerator. User will be asked if none is specified
#' @param ranks Ranks at which to calculate fold change of features
#' @param allRanks If set to TRUE, will compute feature fold change at all ranks.
#'     Defaults to FALSE
#' @param pairSamples Should fold change be calculated on paired samples instead
#'     of between the means of groups? Defaults to FALSE. This option is useful
#'     when you have 2 or more time points from the same patient or something
#'     similar
#'
#' @return Table with the fold change and log2(fold change) for the features
#' @export
#'
foldChange <- function(dataset=NULL, factor=NULL,
                       ref_group=NULL, comparison_groups=NULL,
                       ranks=NULL, allRanks=F,
                       pairSamples=F) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  group_options <- dataset$factors[[factor]]$subset

  ref_group <- ref_group[ref_group %in% group_options]
  if(length(ref_group)!=1) ref_group <- select.list(group_options, graphics=T,
                                                    title='Select a reference group')
  group_options <- group_options[!(group_options %in% ref_group)]

  comparison_groups <- comparison_groups[comparison_groups %in% group_options]
  if(!length(comparison_groups)) comparison_groups <- group_options

  if(allRanks) ranks <- getRanks(dataset)
  else {
    ranks <- ranks[ranks %in% getRanks(dataset)]
    if(is.null(ranks)) ranks <- dataset$data$proc$active_rank
  }

  if(pairSamples) {
    #TODO
    stop('Sorry: Fold-change calculation with paired samples is still under development')
    pairing_factor <- select.list(colnames(metadata),graphics = T,
                                  title = 'Select the factor to pair samples by')
  } else {
    summary_tab <- summarizeTaxa(dataset, factor=factor, ranks=ranks)
    fc_tab <- data.frame(summary_tab[1:2],Reference=rep(ref_group,nrow(summary_tab)))
    for(grp in comparison_groups) {
      fc_tab[[grp]] <- summary_tab[[paste0(grp,' Mean')]]/summary_tab[[paste0(ref_group,' Mean')]]
    }
    fc_tab <- pivot_longer(fc_tab, cols=comparison_groups,
                           names_to='Comparison', values_to='FC')
    fc_tab$Log2FC <- log(fc_tab$FC, base=2)
  }

  return(fc_tab)
}
