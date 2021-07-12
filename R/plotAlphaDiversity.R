#' Plot Alpha Diversity
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param method Method for alpha diversity calculation. One of either "chao1",
#'     "shannon", "simpson", "invsimpson", or "pd". Defaults to "chao1"
#' @param factor Factor along which to compare alpha diversity between groups.
#'     Defaults to the active factor
#' @param stratify Whether or not to stratify. Defaults to FALSE
#' @param facet.x (Optional) A factor to stratify by horizontally
#' @param facet.y (Optional) A factor to stratify by vertically
#' @param flattenFactors Combine multiple factors
#' @param param Whether to perform parametrized or nonparametrized univariate
#'     analysis of alpha diversity. Defaults to FALSE (nonparametrized)
#' @param showStats Whether or not to show significance labels. Defaults to TRUE
#' @param add_xaxis Whether to add x-axis labels. Defaults to FALSE
#' @param separateLegend Whether to separate the legend from the plot. Defaults
#'     to FALSE
#'
#' @return Boxplot of alpha diversity of the groups
#' @export
#'
plotAlphaDiv <- function(dataset=NULL,
                         method='Chao1',
                         factor=NULL,
                         stratify=FALSE, facet.x=NULL, facet.y=NULL,
                         flattenFactors=FALSE,
                         param=FALSE,
                         showStats=TRUE,
                         add_xaxis=F, separateLegend=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  factor <- setFVar(dataset)
  colors <- dataset$colors
  colors <- colors[names(colors) %in% factor$subset]

  # Calculate alpha diversity
  div_tab <- adiv(dataset,method=method)
  div_data <- cleanData(merge(dataset$metadata, div_tab), factor)

  # If alpha diversity was successfully calculated, the "chosen feature"
  #   variable is set to the calculation method
  chosen_ft <- tolower(method)

  # Make the initial plot. Significant difference markers will be added later
  #   if "showStats" is true
  p <- ggboxplot(div_data,x=factor$name,y=tolower(method),color=factor$name,size=1)+
    scale_color_manual(values=colors)+
    scale_y_continuous(expand=expansion(mult=c(.1,.1)))+
    labs(y=paste(capitalize(method),'Index'),
         x=factor$name,title='Alpha Diversity',
         colour=factor$name_text)+
    theme(plot.title = element_text(hjust = 0.5, size=25),
          axis.title.y = element_text(size=25),
          axis.text.y = element_text(size=22),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'top',
          legend.title = element_blank(),
          legend.text = element_text(size=20),
          legend.key.size = unit(3,'line')
    )+
    expand_limits(y=0)

  suffix <- '_adiv_'

  if(add_xaxis) {
    p <- p+theme(axis.text.x = element_text(size=20),
                 axis.title.x = element_blank())
  }

  if(separateLegend) {
    if(!exists('p_legend',inherits = F)) p_legend <- as_ggplot(get_legend(p))
    p <- p+theme(legend.position = 'none')
    suffix <- paste0(suffix,'_nolegend')
  }

  # Figure out if there are any factors the user wants to facet by
  #   If so, then facet the plot by those factor(s)
  facets <- parseStratifiers(factor$name, dataset$factors, stratify, facet.x, facet.y)

  # Perform statistical test on alpha diversity data
  stats <- univar(data=div_data,
                  factor=factor$name,
                  stratifiers=c(facets$x,facets$y),
                  features=tolower(method),
                  param=param,
                  dataset_name=dataset_name)

  if(!is.null(facets$x)) p <- facet(p,facet.by=facets$x,nrow=1)+
    theme(strip.text.x = element_text(size=18))
  if(!is.null(facets$y)) p <- facet(p,facet.by=facets$y,ncol=1)+
    theme(strip.text.y = element_text(size=18))

  # Plot the figure
  if(showStats) {
    if(!is.null(stats$pw_stats)) {
      pw_stats <- stats$pw_stats
      pw_stats$y.position <- 1.1*pw_stats$y.position
      p <- p+stat_pvalue_manual(pw_stats,
                                label='p.adj.signif',
                                step.increase=0.03,
                                bracket.size=1,
                                tip.length=0,
                                hide.ns=T)
    } else if(!is.null(stats$stats)) {
      tot_stats <- stats$stats
      tot_stats$y.position <- 1.1*tot_stats$y.position
      p <- p+stat_pvalue_manual(tot_stats,
                                label='p.adj.signif',
                                label.size = 9,
                                bracket.size = 1,
                                tip.length = 0,
                                hide.ns = T)
    }
  }

  show(p)

  saveResults(dataset$results_path,foldername = 'Alpha Diversity',
              factors = dataset$factors,
              active_factor = factor$name,
              facets = facets,
              stat_results = stats,
              width = 8, height = 6,
              suffix = paste0('_adiv_',method))

  legend_output_location <- paste0(dataset$results_path,'/Results_',Sys.Date(),'/Alpha Diversity/')
  if(exists('p_legend')) ggsave(legend_output_location,
                                filename="Legend.png",
                                plot=p_legend,
                                device = 'png',
                                width = 16,
                                height = 6)

  activate(dataset)
  return(p)
}
