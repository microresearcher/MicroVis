#' Plot stacked abundace barplots
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param relative Whether to plot relative abundance. Defaults to TRUE
#' @param factor Factor to group samples by
#' @param stratify Whether to stratify the groups. Defaults to FALSE
#' @param bySample Whether to show a bar for each sample instead of aggregating
#'     samples by groups. Defaults to FALSE
#' @param allSamples Whether to plot a single graph for all the samples in the
#'     dataset. Defaults to FALSE
#' @param unfiltered Whether to use data before or after filtering out features.
#'     Defaults to TRUE (uses data before filtering)
#' @param label_unknown Whether to group unknown taxa separate. Defaults to TRUE
#' @param rank Rank at which to use features
#' @param top Number of features to show based on greatest relative abundance.
#'     Defaults to top 15 features by relative abundance
#' @param top_by Select the top features either by the max, min, or mean relative
#'     abundance in the groups, or the summed relative abundance across all groups.
#'     Defaults to mean (in other words, for each feature its mean
#'     relative abundance values across the samples will be used to determine
#'     whether it is a top feature).
#' @param ftlist (Optional) List of features to consider. Defaults to all features
#'     at the rank
#' @param plotSigs Plot only the significant features? Defaults to FALSE
#' @param alpha Significance threshold for selecting significant features. Defaults
#'     to 0.05
#' @param separateLegend Whether to separate the legend from the plot. Defaults
#'     to FALSE
#' @param paired_ids Plot samples paired by "factor", with each pair having
#'     a unique combination of IDs determined by this variable (as a vector)
#'
#' @return Stacked bar plot of feature abundances of selected features with the
#'     other features lumped into an "Other" category
#' @export
#'
plotStackedBars <- function(dataset=NULL, relative=T,
                            factor=NULL, stratify=F,
                            bySample=F, allSamples=F, paired_ids = c(),
                            unfiltered=T,label_unknown=T,rank=NA,
                            top=15, top_by='mean',
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

  if(label_unknown) tempds <- filterNAs(tempds,
                                        ranks = c(rank),
                                        temp=T, silent = T)

  data <- mvmelt(tempds)
  factor <- dataset$active_factor
  mdcolnum <- ncol(tempds$metadata)
  fts <- colnames(data[(mdcolnum+1):ncol(data)])
  fts <- fts[!(fts %in% c('Other','Unknown'))]
  paired_ids <- paired_ids[paired_ids %in% colnames(tempds$metadata)]

  if(length(ftlist)) fig_name <- paste0('selected_',rank)
  else fig_name <- rank

  if(plotSigs) {
    sigfts <- listsigs(dataset=dataset,
                       factor=factor,
                       ranks=rank,
                       alpha=alpha,
                       silent=T,
                       dataset_name=dataset_name)
    if(!length(sigfts)) message('\nNo significant features were found\n')
    else fig_name <- paste0('sig-alpha_',alpha,fig_name)
    ftlist <- c(ftlist,sigfts)
  }

  # If pairing samples, this step should occur after getting a list of
  #   significant features, but before manipulating the data for representation
  #   (ie aggregating less common features into the "Other" group and
  #   determining relative abundance).
  if(length(paired_ids)) {
    pairedSamples <- T
    paired_samples <- getSamples(dataset, id_cols = paired_ids, complete = factor)$sample

    data.paired <- data[data$sample %in% paired_samples,]
    cat('\n',nrow(data.paired),'out of',nrow(data),'samples form complete pairs\n')
    data <- data.paired

    id <- paste0(paired_ids,collapse = '_')
    data[[id]] <- apply(data[paired_ids], 1, paste, collapse = '_')
  } else pairedSamples <- F

  ftlist <- ftlist[ftlist %in% fts]
  if(length(ftlist)) {
    if(length(ftlist)>1) {
      fts <- ftlist
    } else {
      data$Other <- rowSums(data[fts[!(fts %in% ftlist)]])
      data[fts[!(fts %in% ftlist)]] <- NULL
      fts <- c(ftlist,'Other','Unknown')
    }
  } else if(label_unknown) fts <- c(fts,'Unknown')

  if(relative) {
    data[fts] <- data.frame(t(apply(data[fts], 1, function(x) x/sum(x))))
    abundance_type <- 'Relative Abundance'
  } else {
    abundance_type <- 'Absolute Abundance'
  }

  if(bySample) {
    data$sample <- as.character(data$sample)
    if(!pairedSamples) id <- 'sample'
    fig_name <- paste0(fig_name,'_bySample')
  } else if(allSamples) {
    data <- data.frame(t(colMeans(data[fts])))
    colnames(data) <- fts
    abundance_type <- paste0('Mean ',abundance_type)
  } else {
    if(stratify) {
      facet <- select.list(names(dataset$factors)[!(names(dataset$factors) %in% factor)],
                           graphics = T)
      if(is.null(facet)) stratify <- F
      else {
        # compareby <- paste(factor,facet,sep = '+')
        data <- data[c(factor,facet,fts)]
        fig_name <- paste0(fig_name,'_stratified')
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
    abundance_type <- paste0('Mean ',abundance_type)
  }

  if(top<length(fts)) {
    if(allSamples | !(top_by %in% c('max','min','mean','sum'))) top_by <- 'mean'
    ftstats <- data.frame(max=apply(data[fts],2,function(x) max(x)),
                          min=apply(data[fts],2,function(x) min(x)),
                          mean=apply(data[fts],2,function(x) mean(x)),
                          sum=apply(data[fts],2,function(x) sum(x)))
    ftstats <- ftstats[!(rownames(ftstats) %in% c('Other','Unknown')),]

    # If unknowns are to be included separately, make sure they don't count
    #   towards the number of named features shown
    fts <- fts[!(fts %in% 'Unknown')]

    if(top_by=='max') low_abun <- rownames(dplyr::slice_min(ftstats, order_by=max, n=(length(fts)-top)))
    if(top_by=='min') low_abun <- rownames(dplyr::slice_min(ftstats, order_by=min, n=(length(fts)-top)))
    if(top_by=='mean') low_abun <- rownames(dplyr::slice_min(ftstats, order_by=mean, n=(length(fts)-top)))
    if(top_by=='sum') low_abun <- rownames(dplyr::slice_min(ftstats, order_by=sum, n=(length(fts)-top)))

    # Add the unknown feature back to the feature list (fts) if it is to be labeled
    #   separately
    if(label_unknown) fts <- c(fts,'Unknown')

    if(is.null(data$Other)) data$Other <- rep(0,nrow(data))
    data$Other <- data$Other + rowSums(data[low_abun])
    data[low_abun] <- NULL
    fts <- c(fts[!(fts %in% low_abun)],'Other','Unknown')
    fig_name <- paste0('top',top,fig_name)
  }

  data_pivoted <- data %>% tidyr::pivot_longer(all_of(fts),
                                               names_to=rank,
                                               values_to=abundance_type)

  namedfts <- unique(data_pivoted[[rank]][!(data_pivoted[[rank]] %in% c('Other','Unknown'))])
  data_pivoted[[rank]] <- factor(data_pivoted[[rank]],
                                 levels=c('Unknown','Other',namedfts))

  if(bySample) {
    p <- ggplot(data_pivoted, aes(x=as.numeric(factor(.data[[id]])), y=.data[[abundance_type]]))+
      geom_area(aes(fill=.data[[rank]]))+
      scale_x_continuous(id, labels = unique(data_pivoted[[id]]),
                         breaks = unique(as.numeric(factor(data_pivoted[[id]]))))+
      ggpubr::theme_pubr()+
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
    if(pairedSamples) {
      p <- p+theme(axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 60, vjust = 0.5, size=15),
                   axis.ticks.x = element_line())
      separateLegend <- T
    }
  } else if(allSamples) {
    p <- ggpubr::ggbarplot(data_pivoted,y=abundance_type,
                           fill=tempds$data$proc$active_rank,color='white')+
      labs(fill=capitalize(rank),
           x=dataset$factors[[factor]]$name_text)+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=22),
            legend.position = 'right',
            legend.title = element_text(size = 22),
            legend.text = element_text(size=15),
            legend.key.size = unit(1,'cm'))
  } else {
    p <- ggpubr::ggbarplot(data_pivoted,x=factor,y=abundance_type,
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
    fig_name <- paste0(fig_name,'_nolegend')
    legend_output_location <- paste0(dataset$results_path,'/Results_',Sys.Date(),'/Stacked Bar Graphs/')

    show(p_legend)
  }

  show(p)

  save_directory <- saveResults(dataset,foldername = 'Stacked Bar Graphs',
                                filename = fig_name,
                                factors = dataset$factors,
                                active_factor = dataset$active_factor,
                                width = 12, height = 8,
                                stat_results = list(stats=data_pivoted),
                                suffix = paste0('_',abundance_type),
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

  if(allSamples) {
    total_props <- data_pivoted[data_pivoted[[rank]] %in% namedfts,]
    colnames(total_props) <- c('Group.1','x')
  } else if(bySample) {
    total_props <- aggregate(data_pivoted[[abundance_type]][data_pivoted[[rank]] %in% namedfts],
                             by=list(data_pivoted[[rank]][data_pivoted[[rank]] %in% namedfts]),
                             mean)
  } else {
    total_props <- aggregate(data_pivoted[[abundance_type]][data_pivoted[[rank]] %in% namedfts],
                             by=list(data_pivoted[[factor]][data_pivoted[[rank]] %in% namedfts]),
                             sum)
  }

  cat('\nTop',top,'features make up:\n')
  if(allSamples | bySample) cat('',paste0(signif(100*sum(total_props$x),3),'%'),
                                'across all samples\n')

  if(!bySample & !allSamples) for(grp in total_props$Group.1) {
    cat('',paste0(signif(100*mean(total_props$x[total_props$Group.1==grp]),
                           3),'%'),'in',grp,'\n')
  }

  if('Unknown' %in% data_pivoted[[rank]]) {
    cat('\n',
        paste0(signif((100*mean(data_pivoted[[abundance_type]][data_pivoted[[rank]]=='Unknown'])),
                      3),'%'),'of',rank,'are unknown')
  }

  cat('\n\n')

  return(p)
}
