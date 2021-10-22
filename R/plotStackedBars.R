#' Plot stacked abundace barplots
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param proportional Whether to plot relative/proportional abundance. Defaults
#'     to TRUE
#' @param factor Factor to group samples by
#' @param stratify Whether to stratify the groups. Defaults to FALSE
#' @param bySample Do not aggregate samples by groups and instead show a bar for
#'     each sample. Defaults to FALSE
#' @param unfiltered Whether to use data before or after filtering out features.
#'     Defaults to TRUE (uses data before filtering)
#' @param rank Rank at which to use features
#' @param top Number of features to show based on greatest relative abundance.
#'     Defaults to top 15 features by relative abundance
#' @param top_by Select the top features either by the max, min, or mean relative
#'     abundance in the groups, or the summed relative abundance across all groups.
#'     Defaults to max (in other words, for each feature the maximum of its
#'     relative abundance values across the groups will be used to determine
#'     whether it is a top feature).
#' @param ftlist (Optional) List of features to consider. Defaults to all features
#'     at the rank
#' @param plotSigs Plot only the significant features? Defaults to FALSE
#' @param alpha Significance threshold for selecting significant features. Defaults
#'     to 0.05
#' @param separateLegend Whether to separate the legend from the plot. Defaults
#'     to FALSE
#'
#' @return Stacked bar plot of feature abundances of selected features with the
#'     other features lumped into an "Other" category
#' @export
#'
plotStackedBars <- function(dataset=NULL, proportional=T,
                            factor=NULL, stratify=F,
                            bySample=F,
                            unfiltered=T,rank=NA,
                            top=15, top_by='max',
                            ftlist=c(),plotSigs=F,alpha=0.05,
                            separateLegend=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  taxaranks <- get('taxaRanks',envir = mvEnv)

  if(unfiltered) tempds <- clearProcessing(dataset, temp = T, silent = T)
  else tempds <- clearNormalization(dataset, temp = T, silent = T)

  if(!(rank %in% intersect(taxaranks,names(tempds$data$proc)))) {
    rank <- tempds$data$proc$active_rank
  } else tempds$data$proc$active_rank <- rank

  data <- mvmelt(tempds)
  factor <- dataset$active_factor
  mdcolnum <- ncol(tempds$metadata)
  fts <- colnames(data[(mdcolnum+1):ncol(data)])
  fts <- fts[fts!='Other']

  if(length(ftlist)) suffix <- paste0('_selected_',rank)
  else suffix <- paste0('_',rank)

  if(plotSigs) {
    sigfts <- listsigs(dataset=dataset,
                       factor=factor,
                       ranks=rank,
                       alpha=alpha,
                       silent=T,
                       dataset_name=dataset_name)
    if(!length(sigfts)) message('\nNo significant features were found\n')
    else suffix <- paste0('_sig-alpha_',alpha,suffix)
    ftlist <- c(ftlist,sigfts)
  }

  ftlist <- ftlist[ftlist %in% fts]
  if(length(ftlist)) {
    if(length(ftlist)>1) {
      fts <- ftlist
    } else {
      data$Other <- rowSums(data[fts[!(fts %in% ftlist)]])
      data[fts[!(fts %in% ftlist)]] <- NULL
      fts <- c(ftlist,'Other')
    }
  }

  if(bySample) {
    data$sample <- as.character(data$sample)
    suffix <- paste0(suffix,'_bySample')
  } else {
    if(stratify) {
      facet <- select.list(names(dataset$factors)[!(names(dataset$factors) %in% factor)],graphics = T)
      if(is.null(facet)) stratify <- F
      else {
        # compareby <- paste(factor,facet,sep = '+')
        data <- data[c(factor,facet,fts)]
        suffix <- paste0(suffix,'_stratified')
      }
    } else {
      data <- data[c(factor,fts)]
    }
    aggformula <- formula(paste('. ~',factor))
    # data <- aggregate(aggformula, data, function(x) sum(x))
    agglist <- list()
    for(ft in fts) agglist[[ft]] <- aggregate(aggformula,data[c(factor,ft)],
                                              function(x) mean(x,na.rm=T))
    data <- purrr::reduce(agglist,merge)
  }

  if(proportional) {
    data[fts] <- data.frame(t(apply(data[fts], 1, function(x) x/sum(x))))
    abundance_type <- 'Proportional Abundance'
  } else {
    abundance_type <- 'Absolute Abundance'
  }

  if(top<length(fts)) {
    if(!(top_by %in% c('max','min','mean','sum'))) top_by <- 'max'
    ftstats <- data.frame(max=apply(data[fts],2,function(x) max(x)),
                          min=apply(data[fts],2,function(x) min(x)),
                          mean=apply(data[fts],2,function(x) mean(x)),
                          sum=apply(data[fts],2,function(x) sum(x)))
    ftstats <- ftstats[!(rownames(ftstats) %in% 'Other'),]

    if(top_by=='max') low_abun <- rownames(slice_min(ftstats, order_by=max, n=(length(fts)-top)))
    if(top_by=='min') low_abun <- rownames(slice_min(ftstats, order_by=min, n=(length(fts)-top)))
    if(top_by=='mean') low_abun <- rownames(slice_min(ftstats, order_by=mean, n=(length(fts)-top)))
    if(top_by=='sum') low_abun <- rownames(slice_min(ftstats, order_by=sum, n=(length(fts)-top)))

    if(is.null(data$Other)) data$Other <- rep(0,nrow(data))
    data$Other <- data$Other + rowSums(data[low_abun])
    data[low_abun] <- NULL
    fts <- c(fts[!(fts %in% low_abun)],'Other')
    suffix <- paste0('_top',top,suffix)
  }

  data_pivoted <- data %>% pivot_longer(fts,
                                        names_to=rank,
                                        values_to=abundance_type)

  namedfts <- unique(data_pivoted[[rank]][!(data_pivoted[[rank]] %in% c('Other'))])
  data_pivoted[[rank]] <- factor(data_pivoted[[rank]],
                                 levels=c('Other',namedfts))

  if(bySample) {
    p <- ggplot(data_pivoted, aes(x=as.numeric(factor(.data$sample)), y=.data[[abundance_type]]))+
      geom_area(aes(fill=.data[[rank]]))+
      theme_pubr()+
      labs(fill=capitalize(rank))+
      facet_wrap(facets = factor, scales = 'free_x',ncol = 1)+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=22),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'right',
            legend.title = element_text(size = 22),
            legend.text = element_text(size=15),
            legend.key.size = unit(1,'cm'),
            strip.text = element_text(size = 18))
  } else {
    p <- ggbarplot(data_pivoted,x=factor,y=abundance_type,
                   fill=tempds$data$proc$active_rank,color='white')+
      labs(fill=capitalize(rank),
           x=dataset$factors[[factor]]$name_text)+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=22),
            legend.position = 'right',
            legend.title = element_text(size = 22),
            legend.text = element_text(size=15),
            legend.key.size = unit(1,'cm'))
  }

  if(stratify) {
    p <- facet(p,facet.by=facet)+
      theme(strip.text = element_text(size = 18))
  }

  if(separateLegend) {
    if(!exists('p_legend',inherits = F)) p_legend <- as_ggplot(get_legend(p))
    p <- p+theme(legend.position = 'none')
    suffix <- paste0(suffix,'_nolegend')
    legend_output_location <- paste0(dataset$results_path,'/Results_',Sys.Date(),'/Stacked Bar Graphs/')

    show(p_legend)
  }

  show(p)

  save_directory <- saveResults(dataset$results_path,foldername = 'Stacked Bar Graphs',
                                factors = dataset$factors,
                                active_factor = dataset$active_factor,
                                width = 12, height = 8,
                                suffix = paste0(suffix,'_',abundance_type),
                                verbose = F)

  if(!is.null(save_directory)) {
    if(exists('p_legend')) {
      ggsave(save_directory,
             filename=paste0('Legend_top',top,'.png'),
             plot=p_legend,
             device = 'png',
             width = 20,
             height = 8,
             dpi=600)
    }
    cat('\nFigure(s) saved to:\n ',save_directory,'\n')
  }

  activate(dataset)

  total_props <- aggregate(data_pivoted$`Proportional Abundance`[data_pivoted[[rank]] %in% namedfts],by=list(data_pivoted[[factor]][data_pivoted[[rank]] %in% namedfts]), sum)

  cat('Top',top,'features make up:\n')
  cat('',paste0(signif(100*mean(total_props$x),3),'%'),
      'across all samples')
  for(grp in total_props$Group.1) cat('\n',paste0(signif(100*mean(total_props$x[total_props$Group.1==grp]),3),'%'),'in',grp)
  cat('\n\n')

  return(p)
}
