#' Plot Cladogram
#'
#' @param dataset MicroVis dataset. Defaults to the activedataset
#' @param factor Factor to analyze by. Defaults to the active factor
#'
#' @return Cladogram plot of LEfSe significant taxa
#' @export
#'
plotClad <- function(dataset=NULL, factor=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  clrs <- dataset$colors

  factor <- factor[factor %in% names(dataset$factor)]
  if(is.null(factor)) factor <- dataset$active_factor
  factor <- setFVar(dataset,factor_name = factor)

  if(is.null(dataset$stats[[factor$name]]$lefse)) dataset <- mvlefse(dataset,
                                                                     dataset_name)

  mm <- dataset$stats[[factor$name]]$lefse

  clrlist <- clrs[names(clrs) %in% factor$subset]
  clad <- plot_cladogram(mm, color = clrlist,
                         clade_label_level = 2,
                         annotation_shape = 22,
                         annotation_shape_size = 3,
                         node_size_offset = 2)

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

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))
  print(dataset)

  return(clad)
}
