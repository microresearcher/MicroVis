#' Plot grouped bar plot of significant taxa at each rank of all groups
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param abun_thresh Abundance cutoff for what features to show. Defaults to 0
#'     (no threshold)
#' @param plotSigs Plot significant features only. Defaults to TRUE
#' @param ftlist (Optional) List of specific features to plot
#' @param addfts (Optional) Additional features to plot in addition to significant
#'     ones
#' @param factor Factor to perform statistical analysis by. Defaults to the
#'     active factor
#' @param stratify Stratify by another factor? Defaults to FALSE
#' @param flattenFactors Whether to combine factors. Defaults to FALSE
#' @param facet.x (Optional) Factor to stratify horizontally
#' @param facet.y (Optional) Factor to stratify vertically
#' @param width Width of the saved figure (in inches). Defaults to 12 inches
#' @param alpha Significance threshold. Defaults to 0.05
#' @param param Whether to perform parametrized or nonparametrized univariate
#'     analysis. Defaults to FALSE (nonparametrized)
#' @param byrank Facet the plot by rank? Defaults to FALSE
#'
#' @return Grouped bar plot of significant features
#' @export
#'
plotGroupedBars <- function(dataset=NULL,
                            abun_thresh=0,
                            plotSigs=T, ftlist=NULL, addfts=NULL,
                            factor=NULL,
                            stratify=FALSE, flattenFactors=FALSE,
                            facet.x=FALSE, facet.y=FALSE,
                            width=12,
                            alpha=0.05, param=F,
                            byrank=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  factor <- dataset$active_factor
  colors <- dataset$colors
  colors <- colors[names(colors) %in% dataset$factors[[factor]]$subset]

  metadata <- dataset$metadata

  dataset.rel <- scaleSamples(clearNormalization(dataset, temp=T, silent=T),
                              scaling='relative', temp=T, silent=T)

  plottab <- data.frame()
  skipped_list <- list()

  ranks <- getRanks(dataset)
  stats_allranks <- list()
  for(rank in ranks) {
    if(rank=='domain') next

    rankfts <- getFeatures(dataset, ranks=rank)

    if(plotSigs) {
      if(is.null(dataset$stats[[factor]]$univar[[rank]])) dataset <- univar(data=dataset,
                                                                            rank=rank,
                                                                            param=param,
                                                                            dataset_name = dataset_name)

      stats <- dataset$stats[[factor]]$univar[[rank]]

      sigfts <- c(stats$stats$.y.[stats$stats$p.adj <= alpha])
      skipped_list[[rank]] <- stats$skipped

      fts <- unique(c(sigfts,addfts[addfts %in% rankfts]))

      if(length(stats$stats)) {
        temp <- cbind(data.frame(Rank=rep(rank,nrow(stats$stats))),
                      stats$stats)

        if(!length(stats_allranks$stats)) stats_allranks$stats <- temp
        else stats_allranks$stats <- rbind(stats_allranks$stats, temp)
      }

      if(length(stats$pw_stats)) {
        temp <- cbind(data.frame(Rank=rep(rank,nrow(stats$pw_stats))),
                      stats$pw_stats)

        if(!length(stats_allranks$pw_stats)) stats_allranks$pw_stats <- temp
        else stats_allranks$pw_stats <- rbind(stats_allranks$pw_stats, temp)
      }

    } else {
      stats <- NULL
      fts <- ftlist[ftlist %in% rankfts]
      if(!length(fts)) fts <- rankfts
    }

    if(length(fts)) {
      # If any significant fts were found at this taxonomy rank, add them to the table
      abd.rel <- mvmelt(dataset.rel, rank=rank, features=fts)

      abd_bygroup <- split(abd.rel,abd.rel[[factor]])
      for(grp in abd_bygroup) {
        grp_name <- as.character(grp[[factor]][1])
        grptab <- data.frame(lapply(grp[fts],function(x) as.numeric(as.character(x))))
        colnames(grptab) <- fts
        means <- colMeans(grptab[fts])
        SE <- unlist(lapply(grptab,function(x) sd(x)/sqrt(length(x))))
        grpsumtab <- data.frame(rep.int(rank,length(fts)),fts,rep.int(grp_name,length(fts)),means,SE)
        row.names(grpsumtab) <- c()
        colnames(grpsumtab) <- c('Rank','Feature','Group','Mean','SE')
        plottab <- rbind(plottab,grpsumtab)
      }
    }
  }
  if(!nrow(plottab)) {
    message('\nNo significant features found for the following dataset:')
    return(dataset)
  }

  pivoted_tab <- plottab
  pivoted_tab$Mean <- sapply(pivoted_tab$Mean,function(x) formatC(x,digits=2,format='e'))
  pivoted_tab$SE <- sapply(pivoted_tab$SE,function(x) formatC(x,digits=2,format='e'))
  pivoted_tab <- pivot_wider(pivoted_tab,names_from = 'Group',values_from = c('Mean','SE'))
  # Identify fts with very low abundance counts for any group and remove them
  lowab <- unique(as.character(plottab[plottab$Mean<abun_thresh,]$Feature))

  plottab <- plottab[!(plottab$Feature %in% lowab),]
  plottab$Feature <- factor(plottab$Feature,levels=unique(plottab$Feature))
  plottab$Group <- factor(plottab$Group,levels=dataset$factors[[factor]]$subset)
  if(byrank) plottab$Rank <- factor(capitalize(as.character(plottab$Rank)),
                                    levels = capitalize(levels(plottab$Rank)))

  p<-ggplot(plottab,aes(x=Feature,y=Mean,fill=Group))+
    geom_bar(position='dodge',stat='identity',width = 0.7)+
    geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE),
                  colour='black',
                  position=position_dodge(.68),
                  size=0.4, width=0.4)+
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust=1))+
    labs(y='Mean Relative Abundance')+
    scale_y_continuous(trans='sqrt',breaks=scales::pretty_breaks(n=6))+
    scale_fill_manual(values=colors)

  p <- p+theme(axis.title = element_text(size=25),
               axis.text = element_text(size=20),
               axis.text.x = element_text(size=22,angle = 60),
               legend.text = element_text(size=22,margin = margin(r=20,unit='pt')),
               legend.title = element_blank(),
               legend.position = 'top',
               legend.key.size = unit(1.5,'line'))

  suffix <- paste0('_alpha_',alpha,'_abun_',abun_thresh)
  if(byrank) {
    p <- p+facet_grid(cols = vars(Rank),scales = 'free_x',space = 'free_x')+
      theme(strip.text.x = element_text(size=18))
    suffix <- paste0(suffix,'_byrank')
  }

  show(p)

  stats_allranks$stats <- stats_allranks$stats[2:nrow(stats_allranks$stats),]
  if(length(stats_allranks$pw_stats)) {
    stats_allranks$pw_stats <- stats_allranks$pw_stats[2:nrow(stats_allranks$pw_stats),]
  }

  save_directory <- saveResults(dataset$results_path,foldername = 'Grouped Bar Graphs',
                                factors = dataset$factors,
                                active_factor = factor,
                                stat_results = stats_allranks,
                                other_results = list(Mean_SE=pivoted_tab),
                                width = width, height = 10,
                                suffix = suffix)

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))
  print(dataset)

  return(p)
}
