#' Plot Cladogram
#'
#' @param dataset MicroVis dataset. Defaults to the activedataset
#' @param factor Factor to analyze by. Defaults to the active factor
#'
#' @return Cladogram plot of LEfSe significant taxa
#' @export
#'
plotClad <- function(dataset=NULL, factor=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  on.exit(cat('There is an incompatibility between ggtree and dplyr 1.0.6 as of May 2021'))


  factor <- factor[factor %in% names(dataset$factor)]
  if(is.null(factor)) factor <- dataset$active_factor
  factor <- setFVar(dataset,factor_name = factor)
  colors <- dataset$colors
  colors <- colors[names(colors) %in% factor$subset]

  if(is.null(dataset$stats[[factor$name]]$lefse)) dataset <- mvlefse(dataset,
                                                                     dataset_name)

  mm <- dataset$stats[[factor$name]]$lefse

  clad <- plot_cladogram(mm, color = clrlist,
                         clade_label_level = 2,
                         annotation_shape = 22,
                         annotation_shape_size = 3,
                         node_size_offset = 2)

  on.exit()

  clad_legend <- as_ggplot(get_legend(clad))
  clad <- clad+theme(legend.position = 'none')

  show(clad)
  show(clad_legend)

  save_directory <- saveResults(dataset$results_path,foldername = 'Cladogram',
                                factors = dataset$factors,
                                active_factor = dataset$active_factor,
                                verbose = F)

  if(!is.null(save_directory)) {
    ggsave(file.path(save_directory,paste0(nameAnalysis(dataset$factors,dataset$active_factor),
                                           'legend.png')),clad_legend,device='png',
           width = 18,height = 12,units = 'in',dpi = 600)
    cat('\nFigures and any associated statistics saved to:\n ',save_directory)
  }

  activate(dataset)

  return(clad)
}
