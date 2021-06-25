#' Plot feature correlation matrix
#'
#' @param dataset1 MicroVis dataset. Defaults to the active dataset
#' @param dataset2 (Optional) MicroVis dataset with features to correlate with
#'     the dataset1
#' @param method Correlation method. One of either "pearson", "spearman", or
#'     "kendall". Defaults to "spearman"
#' @param factor Factor to group samples by
#' @param rank Rank of features to use for calculation
#' @param sigsOnly Select only significant correlations? Defaults to TRUE
#' @param alpha Significance threshold. Defaults to 0.05
#' @param adjustp Use the adjusted p-values? Defaults to TRUE
#' @param r_cutoff Only select correlations with an absolute R-value above a certain
#'     threshold? This value must be between -1 and 1. Defaults to 0
#' @param showlowrs (If r_cutoff is non-zero) Still shade correlations that are
#'     not above the correlation cutoff (r_cutoff)? Defaults to F
#' @param features1 (Optional) Specific features to consider from first dataset.
#'     If only one dataset is being used, then these features will be correlated
#'     with features of the same dataset
#' @param features2 (Optional) Specific features to consider from second dataset.
#' @param matchFts Show the same features in the correlation matrices for all groups?
#'     This can be helpful when comparing correlations between groups. Defaults
#'     to FALSE
#' @param circles Use circles instead of squares to make the matrix? Defaults to
#'     FALSE
#' @param theme Coloring theme. One of "classic" or "modern". Defaults to classic
#'
#' @return List of the one or two datasets used for this plot
#' @export
#'
plotFtCormat <- function(dataset1=NULL, dataset2=NULL,
                         method=c('spearman','pearson','kendall'),
                         factor=NULL,
                         rank=NULL,
                         sigsOnly=T, alpha=0.05, adjustp=T,
                         r_cutoff=0, showlowrs=F,
                         features1=NULL, features2=NULL,
                         matchFts=F,
                         theme=c('classic','modern'),
                         circles=F) {
  if(is.null(dataset)) dataset1 <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset1$name)) dataset1_name <- 'active_dataset'
  else dataset1_name <- dataset1$name

  if(!is.null(dataset2)) {
    dataset2_name <- dataset2$name
    filename <- paste(dataset1_name,dataset2_name,collapse = '_and_')
  } else filename <- dataset1$data$proc$active_rank

  if(length(method) > 1) method <- 'spearman'
  # Function will accept a positive or negative number
  #   as long as its absolute value is between 0 and 1
  r_cutoff <- abs(r_cutoff)
  if(!(r_cutoff>0 & r_cutoff<1)) r_cutoff <- 0

  if(circles) shape <- 'circle'
  else shape <- 'square'

  theme <- match.arg(theme)
  if(theme=='modern') {
    color.low <- '#00daf2'
    color.center <- '#363636'
    color.high <- '#ed4d3b'
    outline.color <- '#b8b8b8'
  } else {
    color.low <- '#003cff'
    color.center <- 'white'
    color.high <- '#ff2f00'
    outline.color <- '#383838'
  }

  cordata <- ftcor(dataset1=dataset1,dataset2=dataset2,
                   method=method,
                   features1=features1,features2=features2)

  cordata <- sliceCormats(cordata,
                          sigsOnly=sigsOnly, alpha=alpha, adjustp=adjustp,
                          r_cutoff=r_cutoff,
                          matchFts=matchFts)

  paramtxt <- ''
  if(length(features1) | length(features2)) paramtxt <- paste0('_specificfts')
  if(sigsOnly) paramtxt <- paste0(paramtxt,'_alpha_',alpha)
  if(r_cutoff!=0) paramtxt <- paste0(paramtxt,'_r2cutoff_',r_cutoff)

  for(grp in names(cordata)) {
    grp_txt <- gsub('_',' ',grp)
    cormat <- cordata[[grp]]$r
    if(adjustp) pmat <- cordata[[grp]]$q
    else pmat <- cordata[[grp]]$P

    if(sigsOnly) cormat[pmat>alpha] <- 0
    if(!showlowrs) cormat[abs(cormat)<r_cutoff] <- 0

    p <- ggcorrplot(cormat,
                    title=grp_txt,legend.title='R',
                    outline.color=outline.color,
                    method=shape,
                    colors = c(color.low,color.center,color.high))
    p <- p+theme(plot.title = element_text(size=25,hjust = 0.5),
                 axis.text = element_text(size=20),
                 axis.text.x = element_text(angle = 75))

    if(!exists('save_one_all',inherits = F)) save_one_all <- NULL
    save_one_all <- multisave(save_one_all)

    if(save_one_all %in% c('Yes','Yes to all figures')) saveFig <- T
    else saveFig <- F

    show(p)
    if(saveFig) {
      save_directory <- saveResults(dataset1$results_path,foldername = 'Correlations',
                                    filename = paste0(filename,paramtxt,'_',grp_txt),
                                    width = 7, height = 10,
                                    forcesave = T,
                                    verbose = F)
    }
  }
  if(exists('save_directory')) cat('Figures and any associated statistics saved to:\n ',save_directory)

  return(list(dataset1,dataset2))
}

#' Plot correlation trendlines between different features
#'
#' @param dataset1 MicroVis dataset. Defaults to the active dataset
#' @param dataset2 (Optional) MicroVis dataset with features to correlate with
#'     the dataset1
#' @param method Correlation method. One of either "pearson", "spearman", or
#'     "kendall". Defaults to "spearman"
#' @param factor Factor to group samples by
#' @param features1  (Optional) Specific features to consider from first dataset.
#'     If only one dataset is being used, then these features will be correlated
#'     with themselves
#' @param features2 (Optional) Specific features to consider from second dataset,
#'     If no second dataset is provided, then `features1` will be correlated with
#'     `features2` from `dataset1`
#' @param sigsOnly Only plot significant correlations (after adjusting for multiple
#'     comparisons)? Defaults to TRUE
#' @param alpha Significance threshold. Defaults to 0.05
#' @param r_cutoff Correlation coefficient threshold. Defaults to 0 (all correlation
#'     values will be accepted when choosing which correlations to plot)
#'
#' @return NULL
#' @export
#'
plotFtCorlines <- function(dataset1=NULL,dataset2=NULL,
                           method=c('pearson','kendall','spearman'),
                           factor=NULL,
                           features1=NULL,features2=NULL,
                           sigsOnly=T,alpha=0.05,
                           r_cutoff=0.5) {
  if(is.null(dataset1)) dataset1 <- get('active_dataset',envir = mvEnv)

  #TODO:Should we force dataset1 and dataset2 into relative abundance?

  cordata <- ftcor(dataset1,dataset2,
                   method=method,
                   factor=factor,
                   features1=features1,features2=features2)

  r_cutoff <- abs(r_cutoff)
  if(!(r_cutoff>0 & r_cutoff<1)) r_cutoff <- 0

  if(!sigsOnly) alpha <- 1

  plottab <- data.frame(Feature1=0,Feature2=0)[0,]
  for(grp in names(cordata)) {
    if(grp=='All_Groups') next
    cortab <- cordata[[grp]]$r
    ptab <- cordata[[grp]]$q

    if(sigsOnly) cortab[ptab>alpha] <- 0
    if(r_cutoff!=0) cortab[abs(cortab)<r_cutoff] <- 0

    cortab$Feature1 <- rownames(cortab)
    cortab <- cortab %>% pivot_longer(cols=colnames(cortab[1:(ncol(cortab)-1)]),
                                      names_to='Feature2',values_to='r2')

    ptab$Feature1 <- rownames(ptab)
    ptab <- ptab %>% pivot_longer(cols=colnames(ptab[1:(ncol(ptab)-1)]),
                                  names_to='Feature2',values_to='p')

    cortab <- merge(cortab,ptab,by=c('Feature1','Feature2'))
    plottab.tmp <- cortab[cortab$p<=alpha & abs(cortab$r2)>=r_cutoff,][c('Feature1','Feature2')]

    plottab <- unique(rbind(plottab,plottab.tmp))
  }

  if(!is.null(dataset2)) dataset <- mvmerge(dataset1,dataset2)
  else dataset <- dataset1

  factor <- setFVar(dataset)
  abun <- mvmelt(dataset)

  if(nrow(plottab)>25) {
    # message('\nMore than 25 correlations identified with specified thresholds\n Would you like to plot them all?')
    # to_plot <- select.list(c('All',as.list(data.frame(t(plottab)))),multiple = T,
    #                        title = 'Select correlations to plot')
    # if('All' %in% to_plot) plottab <- plottab
    # else {
    #
    # }
    message('\nMore than 25 correlations identified with specified thresholds\n Would you like to plot them all or exit and rerun with more stringent alpha/r2 cutoffs?')
    to_plot <- select.list(c('Plot all','Exit'), title = paste('\nPlot all',nrow(plottab),'or Exit?'))
    if(to_plot=='Exit') return(list(dataset1,dataset2))
  }

  for(pair in 1:nrow(plottab)) {
    ft1 <- plottab[pair,'Feature1']
    ft2 <- plottab[pair,'Feature2']

    labelypos <- rep(1.1*max(abun[c(ft1,ft2)]),length(factor$subset))
    labelypos <- sapply(1:length(labelypos),function(x) labelypos[x]*(1+0.07*(x-1)))

    p <- ggscatter(abun,x=ft1,y=ft2,color=factor$name,shape=factor$name,size=4,
                   add='reg.line',conf.int=T,add.params=list(size=1.2))+
      stat_cor(aes(color=get(factor$name)),size=7,
               label.x.npc=c('left'),label.y=labelypos,
               show.legend=F)+
      scale_color_manual(values=dataset$colors)+
      scale_fill_manual(values=dataset$colors)+
      theme(axis.title = element_text(size=25),
            axis.text = element_text(size=20),
            legend.title = element_blank(),
            legend.key.size = unit(1,'cm'),
            legend.text = element_text(size=20))
    show(p)

    ### Save the Results ###
    #----------------------#
    if(!exists('save_one_all',inherits = F)) save_one_all <- NULL
    save_one_all <- multisave(save_one_all)

    if(save_one_all %in% c('Yes','Yes to all figures')) saveFig <- T
    else saveFig <- F

    if(saveFig) {
      save_directory <- saveResults(dataset$results_path,
                                    foldername = file.path('Correlations',
                                                           paste0(dataset1_name,
                                                                  '_',
                                                                  dataset2_name)),
                                    filename = paste0(ft1,'-',ft2),
                                    factors = dataset$factors,
                                    active_factor = factor$name,
                                    width = 8, height = 6,
                                    forcesave = T,
                                    verbose = F)
    }
  }
}
