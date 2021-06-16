#' Plot random forest important features
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor to analyze by
#' @param hideShadow Hide Boruta shadow features? Defaults to TRUE
#' @param confirmedOnly Only plot confirmed features? Defaults to FALSE
#' @param top How many top features to show. Defaults to 15
#' @param max_runs Max Boruta runs. Defaults to 100
#' @param roughfix Roughfix indeterminate Boruta features? Defaults to FALSE
#' @param alpha Significance threshold. Defaults to 0,01
#'
#' @return Plot of important features determined by the Boruta algorithm on
#'     a random forest model
#' @export
#'
plotRFImp <- function(dataset=NULL,
                      factor=NULL,
                      hideShadow=T,confirmedOnly=F,top=15,
                      max_runs=100,roughfix=F,
                      alpha=0.01) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- setFVar(dataset)
  colors <- dataset$colors
  colors <- colors[names(colors) %in% factor$subset]
  colors <- c('gray',colors)
  names(colors)[1] <- 'OOB'
  outlines <- c('#32a83c', '#deb42c', '#cc1f08', '#11007d')
  fills <- c('#32a83c66', '#deb42c66', '#cc1f0866', '#11007d66')
  names(fills) <- names(outlines) <- c('Confirmed','Tentative','Rejected','Shadow')

  rf <- rfboruta(dataset,
                 factor=factor,
                 max_runs=max_runs,
                 roughfix=roughfix,
                 alpha=alpha)

  # First get error rates for each group and plot them
  class_errs <- rf$ClassErrors
  errs <- rf$Errors

  errs$Trees <- as.numeric(rownames(errs))
  errs <- errs %>% pivot_longer(cols=c('OOB',factor$subset),names_to='Class',values_to='Error')

  p_err <- ggline(errs,x='Trees',y='Error',color='Class',plot_type = 'l',size = 1,numeric.x.axis=T)+
    scale_color_manual(values=colors)+
    theme(plot.title = element_text(size=25,hjust = 0.5),
          axis.text = element_text(size=20),
          axis.title = element_text(size = 20),
          axis.title.y = element_text(margin = ggplot2::margin(r=15)),
          legend.title = element_blank(),
          legend.text = element_text(size = 20, margin = ggplot2::margin(r=20)),
          legend.key.size = unit(1,'cm'))

  show(p_err)

  # Now get the important features as identified by the Boruta algorithm
  boruta <- rf$Boruta

  importance <- aggregate(Importance~Feature,boruta,function(x) median(x))
  importance <- unique(merge(importance,boruta[c('Feature','Decision')]))
  importance <- importance[order(-importance$Importance),]

  fts <- importance$Feature
  suffix <- paste0('_rf_alpha_',alpha)
  if(hideShadow) fts <- fts[!(fts %in% c('shadowMin','shadowMean','shadowMax'))]
  if(confirmedOnly) {
    fts <- fts[importance$Decision=='Confirmed']
    titletxt <- paste('Significantly Important Features')
    suffix <- paste0(suffix,'_confirmed')
  } else if(top>3 & top<nrow(importance)) {
    fts <- fts[1:top]
    titletxt <- paste('Top',top,'Important Features')
    suffix <- paste0(suffix,'_top_',top)
  } else titletxt <- paste('Random Forest Importance of All Features')
  boruta <- boruta[boruta$Feature %in% fts,]

  outlines <- outlines[names(outlines) %in% as.character(unique(boruta$Decision))]
  fills <- fills[names(fills) %in% as.character(unique(boruta$Decision))]

  p_imp <- ggboxplot(boruta,x='Feature',y='Importance',fill='Decision',color='Decision',size=1)+
    labs(title = titletxt)+
    scale_fill_manual(values = fills)+
    scale_color_manual(values = outlines)+
    theme(plot.title = element_text(size = 25,hjust = 0.5),
          axis.text = element_text(size=20),
          axis.text.x = element_text(angle = 75,hjust = 1),
          axis.title = element_text(size = 20),
          axis.title.y = element_text(margin = ggplot2::margin(r=15)),
          legend.title = element_blank(),
          legend.text = element_text(size = 20, margin = ggplot2::margin(r=20)),
          legend.key.size = unit(1,'cm'))

  if(confirmedOnly) p_imp <- p_imp+theme(legend.position = 'none')

  show(p_imp)

  savedirectory <- saveResults(dataset$results_path,foldername = 'Random Forest',
                               factors = dataset$factors,
                               active_factor = factor$name,
                               suffix = suffix,
                               width = 15, height = 10)

  if(!is.null(savedirectory)) write.csv(x=importance,
                                        file=file.path(savedirectory,'Importance Values.csv'),
                                        row.names=F)

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))
  print(dataset)

  return(p_imp)
}
