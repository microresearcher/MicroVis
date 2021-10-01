#' Plot distribution curves of feature abundances in each sample
#'
#' @param dataset MicroVis dataset. Default is the active dataset
#' @param bySample Split curves by sample?
#'
#' @return Density plots of feature abundances of each sample
#' @export
#'
plotSampleDensity <- function(dataset=NULL,bySample=T) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  if(bySample) {
    colorby <- 'Sample'
    opaqueness <- 0.02
  } else {
    colorby <- 'turquoise'
    opaqueness <- 0.5
  }

  ft_data <- dataset$data
  lowest_rank <- getLowestRank(dataset)

  histocounts <- ft_data$proc[[lowest_rank]]
  histocounts$Other <- NULL
  ftnames <- colnames(histocounts)
  histocounts$Sample <- rownames(histocounts)
  histocounts <- histocounts %>% pivot_longer(cols = ftnames,
                                              names_to = 'Feature',
                                              values_to = 'Normalized Abundance Values')

  p_histocounts <- ggdensity(histocounts,x='Normalized Abundance Values',y='..count..',
                             color=colorby,size=1,
                             fill = 'lightblue',alpha=opaqueness)+
    labs(title='Sample Abundance Distributions',y='Percent of Features')+
    theme(plot.title = element_text(hjust=0.5,size=25),legend.position = 'none',
          axis.title = element_text(size=22), axis.text = element_text(size=18),
          axis.ticks.length = unit(0.25,'cm'))+
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

  return(p_histocounts)
}

#' Plot distribution curves of abundances of each feature across all samples
#'
#' @param dataset MicroVis dataset. Default is the active dataset
#' @param byFeature Split curves by feature?
#' @param rank Rank at which to select features to plot
#' @param ftlist Features to plot
#'
#' @return Density plots of each feature across all samples
#' @export
#'
plotFeatureDensity <- function(dataset=NULL,byFeature=T,
                               rank=NULL,ftlist=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  if(byFeature) {
    colorby <- 'Feature'
    opaqueness <- 0.15
  } else {
    colorby <- 'lightblue'
    opaqueness <- 0.7
  }

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- getLowestRank(dataset)

  ftlist <- ftlist[ftlist %in% getFeatures(dataset,ranks=rank)]
  if(is.null(ftlist)) ftlist <- getFeatures(dataset,ranks=rank)

  histocounts <- dataset$data$proc[[rank]][ftlist]
  histocounts$Other <- NULL

  samplenames <- rownames(histocounts)
  histocounts <- data.frame(t(histocounts))
  colnames(histocounts) <- samplenames

  histocounts$Feature <- rownames(histocounts)
  histocounts <- histocounts %>% pivot_longer(cols = samplenames,
                                              names_to = 'Sample',
                                              values_to = 'Normalized Abundance Values')

  p_histocounts <- ggdensity(histocounts,x='Normalized Abundance Values',y='..count..',
                             color='#03fcbe33',size=1,
                             fill = colorby,alpha=opaqueness)+
    labs(title='Feature Abundance Distributions',y='Percent of Samples')+
    theme(plot.title = element_text(hjust=0.5,size=25),legend.position = 'none',
          axis.title = element_text(size=22), axis.text = element_text(size=18),
          axis.ticks.length = unit(0.25,'cm'))+
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))

  return(p_histocounts)
}
