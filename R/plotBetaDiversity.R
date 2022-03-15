#' Plot Beta Diversity
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param method Method for dissimilarity calculation. One of either "bray",
#'     "euclidean", "jaccard", "unifrac", "manhattan", "canberra", "clark",
#'     "kulczynski", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#'     "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord"
#' @param weighted (Optional) If method is set to "unifrac", whether to perform
#'     weighted or unweighted unifrac. Defaults to FALSE
#' @param factor Factor along which to group samples for analysis. Defaults to
#'     the active factor
#' @param stratify Another factor to stratify groups
#' @param spokes Draw spokes from group centroids to each sample? Defaults to TRUE
#' @param ellipse Draw a confidence interval ellipse for each group? Defaults to
#'     TRUE
#' @param ci Confidence interval level for ellipses (value between 0 and 1).
#'     Defaults to 0.95
#' @param showStats Show statistical results in the plot. Defaults to TRUE
#' @param labelSamples Label the samples. Defaults to FALSE
#' @param separateLegend Separate the lenged from the figure. Defaults to FALSE
#'
#' @return Scatter plot of the dissimilarity indices of the samples of each group
#' @export
#'
plotBetaDiv <- function(dataset=NULL,
                        method='bray',
                        weighted=F,
                        factor=NULL,
                        stratify=F,
                        spokes=T,
                        ellipse=T,
                        ci=0.95,
                        showStats=T,
                        labelSamples=F,
                        separateLegend=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  factor <- setFVar(dataset)
  colors <- dataset$colors
  colors <- colors[names(colors) %in% factor$subset]
  metadata <- dataset$metadata

  if(stratify & length(names(dataset$factors))>1) {
    cat('\nPlease select a factor to stratify by')
    stratifier <- select.list(names(dataset$factors)[!(names(dataset$factors) %in% factor$name)],
                              title = 'Stratifier',graphics = T)
  }
  if(!exists('stratifier',inherits = F)) stratifier <- NULL

  # Calculate beta-diversity
  bdiv_results <- bdiv(dataset,stratifier=stratifier,method=method)
  coord_tab <- bdiv_results$dist$coord_tab
  stats <- list()
  stats$stats <- bdiv_results$pnova
  if(length(stats)) {
    show(stats$stats)
    cat('\n')
  }

  # Merge coordinate data with metadata
  coord_data <- cleanData(merge(dataset$metadata, coord_tab), factor)

  if(is.null(coord_data)) return(message('ERROR: Something went wrong with beta-diversity analysis'))

  suffix <- paste0('_bdiv_',method)

  if(is.null(stratifier)) {
    p <- ggscatter(coord_data,x='Axis.1',y='Axis.2',
                   color=factor$name,size = 3,
                   ellipse = ellipse,ellipse.level = ci,ellipse.alpha = 0,
                   star.plot = spokes,
                   add.params = list(lwd=1))+
      scale_color_manual(values=colors)+
      labs(colour=factor$name_text)
  }
  else {
    stratified_factor <- paste(factor$name,stratifier,sep = '_')
    coord_data[[stratified_factor]] <- interaction(coord_data[[factor$name]],
                                                   coord_data[[stratifier]],
                                                   sep = ' ')
    reordered <- levels(coord_data[[stratified_factor]])[order(as.character(levels(coord_data[[stratified_factor]])))]

    coord_data[[stratified_factor]] <- factor(coord_data[[stratified_factor]],
                                              levels = reordered)

    p <- ggscatter(coord_data, x='Axis.1',y='Axis.2',
                   color=stratified_factor,shape=stratified_factor,size = 3,
                   ellipse = ellipse,ellipse.level = ci,ellipse.alpha = 0,
                   star.plot = spokes, star.plot.lwd = 1)+
      labs(colour=gsub('_',' and ',stratified_factor))+
      guides(shape='none')

    suffix <- paste0(suffix,'_by_',stratifier)
  }

  suffix <- paste0(suffix,'_',dataset$data$proc$active_rank)

  if(showStats) {
    r2 <- paste('R2:', round(stats$stats$R2[1],3))
    pval <- paste('Pr(>F):',stats$stats$`Pr(>F)`[1])

    x_range <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
    y_range <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range

    p <- p+coord_cartesian(xlim=c(x_range[[1]],
                                  x_range[[2]]+0.2*(x_range[[2]]-x_range[[1]])),
                           ylim=c(y_range[[1]],
                                  y_range[[2]]+0.2*(y_range[[2]]-y_range[[1]])),
                           clip='off')+
      annotation_custom(grob = textGrob(label=r2,hjust=0.5,gp=gpar(cex=2)),
                        ymin=1.2*y_range[[2]],
                        ymax=1.2*y_range[[2]],
                        xmin=x_range[[2]],
                        xmax=x_range[[2]])+
      annotation_custom(grob = textGrob(label=pval,hjust=0.5,gp=gpar(cex=2)),
                        ymin=y_range[[2]],
                        ymax=y_range[[2]],
                        xmin=x_range[[2]],
                        xmax=x_range[[2]])
  }

  p <- p+theme(axis.title = element_text(size = 20),
               axis.text = element_text(size = 20),
               legend.key.size = unit(1,'cm'),
               legend.text = element_text(size = 20,margin = margin(r=20)),
               legend.title = element_blank())+
    guides(fill='none')

  if(labelSamples) {
    p<-p+geom_label_repel(aes(label=sample))
    suffix <- paste0(suffix,'_labeled')
  }

  if(separateLegend) {
    if(!exists('p_legend',inherits = F)) p_legend <- as_ggplot(get_legend(p))
    p <- p+theme(legend.position = 'none')
    suffix <- paste0(suffix,'_nolegend')
    legend_output_location <- paste0(dataset$results_path,'/Results_',Sys.Date(),'/Beta Diversity/')
    if(exists('p_legend')) ggsave(legend_output_location,
                                  filename="Legend.png",
                                  plot=p_legend,
                                  device = 'png',
                                  width = 16,
                                  height = 6)
  }

  show(p)

  saveResults(dataset,foldername = 'Beta Diversity',
              factors = dataset$factors,
              active_factor = factor$name,
              stat_results = stats,
              width = 10, height = 8,
              suffix = suffix)

  activate(dataset)

  return(p)
}
