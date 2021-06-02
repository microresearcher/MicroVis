#' Plot box-plots of univariate analysis
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param jitter Show "jitter" points for each sample? Defaults to TRUE
#' @param dotplot Create dotplot? Defaults to FALSE
#' @param violin Create violinplot? Defaults to FALSE
#' @param showStats Show bars w/asterisks for significance labels? Defaults to
#'     TRUE
#' @param raw Show raw abundance values? Defaults to FALSE. If TRUE, showStats
#'     is forced to FALSE
#' @param stratify Stratify the data along a second factor? Defaults to FALSE
#' @param facet.x (Optional) Name of a factor to stratify horizontally
#' @param facet.y (Optional) Name of a factor to stratify vertically
#' @param flattenFactors (Not ready) Combine two factors into one
#' @param rank Rank at which to select features
#' @param ftlist List of specific features to plot. Defaults to all features at
#'     the given rank
#' @param unique_groups (Optional) Only plot features that are uniquely over/under-
#'     expressed in these groups
#' @param plotAll Plot all features regardless of significance? Defaults to FALSE
#' @param alpha Significance threshold. Defaults to 0.05
#' @param param Perform parametrized analysis? Defaults to FALSE (non-parametrized
#'     analysis by default)
#' @param scalePlot Scale plot by log2? Defaults to FALSE
#' @param add_xaxis Add group labels to the x-axis? Defaults to FALSE
#' @param separateLegend Make one separate figure of just the legend instead of
#'     including it in all the boxplots? Defaults to FALSE
#'
#' @return MicroVis dataset
#' @export
#'
plotUnivar <- function(dataset=NULL,
                       jitter=T,
                       dotplot=F,
                       violin=F,
                       showStats=T,
                       raw=F,
                       stratify=F, facet.x=NULL, facet.y=NULL,
                       flattenFactors=F,
                       rank=NULL,
                       ftlist=NULL, unique_groups=NULL, plotAll=F,
                       alpha=0.05, param=F,
                       scalePlot=F,
                       add_xaxis=F, separateLegend=F) {
  #TODO: Need to deal with when flattenFactors is true

  ### Load Dataset ###
  #------------------#
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor
  factor <- setFVar(dataset, factor_name=factor)

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  clrs <- dataset$colors

  if(raw) {
    dataset.raw <- clearNormalization(dataset, temp=T, silent=T)
    melted <- mvmelt(dataset.raw, rank=rank)
    showStats <- FALSE
    israw <- ' (Raw)'
    isnrml <- ''
  } else {
    melted <- mvmelt(dataset,rank=rank)
    israw <- ''
    isnrml <- 'Normalized '
  }

  melted$Other <- NULL

  allfts <- getFeatures(dataset,rank=rank)
  # Make sure ftlist only has features that actually exist (case-insensitive)
  ftlist <- allfts[tolower(allfts) %in% tolower(ftlist)]

  # Determine if there are any factors to facet the plot by
  facets <- parseStratifiers(factor$name, dataset$factors, stratify, facet.x, facet.y)

  ### Get statistics for the chosen feature type ###
  #------------------------------------------------#
  if(!length(facets)) stats <- dataset$stats[[factor$name]][[rank]]
  else stats <- dataset$stats[[factor$name]][[facets$txt]][[rank]]

  if(is.null(stats)) {
    dataset <- univar(dataset=dataset,
                      stratifiers = c(facets$x,facets$y),
                      rank = rank,
                      param = param,
                      dataset_name = dataset_name)
    if(!length(facets)) stats <- dataset$stats[[factor$name]][[rank]]
    else stats <- dataset$stats[[factor$name]][[facets$txt]][[rank]]
  }
  sigfts <- unique(stats$stats$.y.[stats$stats$p.adj<=alpha])

  ### Determine Which Features to Plot ###
  #--------------------------------------#
  if(plotAll) fts <- allfts
  else if(length(ftlist)) fts <- ftlist
  else if(length(sigfts)) {
    unique_groups <- unique_groups[unique_groups %in% factor$subset]
    if(length(unique_groups)) {
      uniques <- listUniques(dataset,dataset_name=dataset_name)
      sigfts <- unlist(uniques[unique_groups])
    }
    fts <- sigfts
  }
  else {
    print(dataset)
    message('\nNo significant features were found for this dataset at a significance threshold of ',alpha,'.\n')
    return(NULL)
  }

  ### Generate Figures for the Selected Features ###
  #------------------------------------------------#
  if(plotAll | length(ftlist)) {
    cat('\n\nGenerating box-plots for',length(fts),'features:')
  } else {
    cat(paste('\n\nGenerating box-plots for',length(fts),'significant features:\n'))
  }

  # Loop through and plot the desired features
  for(ft in fts) {
    # Create datatable for the feature
    ftTab <- melted[c(ft, names(dataset$factors))]

    colnames(ftTab)[1] <- 'Abundance'

    addlist <- c()
    addlistparam <- list()
    if(dotplot) {
      addlist <- c(addlist, 'dotplot')
      addlistparam[['size']] <- 1
      addlistparam[['alpha']] <- 0.5
    } else if(jitter) {
      addlist <- c(addlist, 'jitter')
      addlistparam[['size']] <- 5
      addlistparam[['alpha']] <- 0.5
    }

    # Graph the abundance data for each significant feature
    if(violin) {
      p <- ggviolin(ftTab,x=factor$name,y='Abundance',color=factor$name,size=1,
                    add = c('mean_sd'))+
        scale_color_manual(values=clrs)+
        labs(y=paste0(isnrml,'Abundance'),title=paste0(gsub('\\.',' ',ft),israw),colour=factor$name_text)+
        theme(plot.title = element_text(hjust = 0.5,size = 24),
              axis.title.y = element_text(size=20),
              axis.text.y = element_text(size=20),
              legend.position = 'none',
              axis.title.x = element_blank(),
              axis.text.x = element_text(size=20))+
        expand_limits(y=0)
    } else {
      p <- ggboxplot(ftTab,x=factor$name,y='Abundance',color=factor$name,size=1,
                     add = addlist, add.params = addlistparam)+
        scale_color_manual(values=clrs)+
        labs(y=paste0(isnrml,'Abundance'),title=paste0(gsub('\\.',' ',ft),israw),colour=factor$name_text)+
        theme(plot.title = element_text(hjust = 0.5,size = 30),
              axis.title.y = element_text(size=25),
              axis.text.y = element_text(size=22),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'top',
              legend.title = element_blank(),
              legend.text = element_text(size=25,margin=margin(r=12)),
              legend.key.size = unit(3,'line'))+
        expand_limits(y=0)
    }

    if(add_xaxis) {
      p <- p+theme(axis.title.x = element_text(size=25),
                   axis.text.x = element_text(size=22))
    }

    if(separateLegend) {
      if(!exists('p_legend',inherits = F)) p_legend <- as_ggplot(get_legend(p))
      p <- p+theme(legend.position = 'none')
      suffix <- '_nolegend'
    } else suffix <- ''

    ### Facet the Plot if Applicable ###
    #----------------------------------#
    if(!is.null(facets$x)) p <- facet(p,facet.by=facets$x,nrow=1)+
      theme(strip.text.x = element_text(size=18))
    if(!is.null(facets$y)) p <- facet(p,facet.by=facets$y,ncol=1)+
      theme(strip.text.y = element_text(size=18))

    ### Add Significance Markers ###
    #------------------------------#
    # If statistical data is to be shown, AND this feature wasn't skipped in
    #   statistical analysis, then add the statistical data to the plot
    if(showStats & !(ft %in% stats$skipped)) {
      if(!is.null(stats$pw_stats)) {
        pw_stats <- stats$pw_stats[stats$pw_stats$.y.==ft,]
        # pw_stats$y.position <- 1.1*pw_stats$y.position
        p <- p+stat_pvalue_manual(pw_stats,
                                  label='p.adj.signif',
                                  label.size = 9,
                                  bracket.size = 1,
                                  tip.length = 0,
                                  step.increase = 0.1,
                                  hide.ns = T)
      } else {
        tot_stats <- stats$stats[stats$stats$.y.==ft,]
        # tot_stats$y.position <- 1.1*tot_stats$y.position
        p <- p+stat_pvalue_manual(tot_stats,
                                  label='p.adj.signif',
                                  label.size = 9,
                                  bracket.size = 1,
                                  tip.length = 0,
                                  step.increase = 0.1,
                                  hide.ns = T)
      }
    }
    show(p)

    ### Save the Results ###
    #----------------------#
    if(!exists('save_one_all',inherits = F)) save_one_all <- NULL
    save_one_all <- multisave(save_one_all)

    if(save_one_all %in% c('Yes','Yes to all figures')) saveFig <- T
    else saveFig <- F

    if(saveFig) {
      save_directory <- saveResults(dataset$results_path,
                                    foldername = paste0('Boxplots_',rank),
                                    filename = ft,
                                    factors = dataset$factors,
                                    active_factor = factor$name,
                                    facets = facets,
                                    suffix = paste0(israw,suffix),
                                    width = 8, height = 6,
                                    forcesave = T,
                                    verbose = F)
    }
  }

  if(exists('save_directory')) {
    if(exists('p_legend')) ggsave(filename=file.path(save_directory,'Legend.png'),
                                  plot=p_legend,
                                  device = 'png')

    dir.create(file.path(save_directory,'Statistics'),showWarnings = FALSE)
    if(!is.null(stats$stats)) {
      write.csv(apply(stats$stats,2,function(x) as.character(x)),
                file=file.path(save_directory,'Statistics','Overall_Statistics.csv'),
                row.names = F)
    }
    if(!is.null(stats$pw_stats)) {
      write.csv(apply(stats$pw_stats,2,function(x) as.character(x)),
                file=file.path(save_directory,'Statistics','Pairwise_Statistics.csv'),
                row.names = F)
    }

    cat('Figures and any associated statistics saved to:\n ',save_directory)
  }

  cat(paste('\n\nSuccessfully plotted:\n',paste(fts,collapse = '\n '),'\n\n'))

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))
  return(dataset)
}
