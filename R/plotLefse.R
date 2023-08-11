#' Plot LEfSe analysis (not cladogram plot)
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor to analyze data by. Defaults to the active factor
#' @param alpha Significance threshold. Defaults to 0.05
#' @param lda_cutoff Effect size cutoff for features to plot. Defaults to 2
#' @param top Top number of features by effect size to plot. Defaults to 20
#' @param byrank Facet plot by rank? Defaults to FALSE
#'
#' @return Plot of effect sizes of significant features from LEfSe analysis
#' @export
#'
plotLEFSE <- function(dataset=NULL,factor=NULL,alpha=0.05,lda_cutoff=2,top=20,byrank=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  factor <- setFVar(dataset)
  colors <- dataset$colors
  colors <- colors[names(colors) %in% factor$subset]

  if(is.null(dataset$stats[[factor$name]]$lefse)) dataset <- mvlefse(dataset,dataset_name)

  mm <- dataset$stats[[factor$name]]$lefse
  lefse_results <- data.frame(Feature=mm@marker_table$feature,
                              Enriched=mm@marker_table$enrich_group,
                              LDA_Score=mm@marker_table$ef_lda,
                              p=mm@marker_table$pvalue,
                              q=mm@marker_table$padj)

  lefse_results$Feature <- gsub('.*\\|','',lefse_results$Feature)
  plottab <- subset(lefse_results,p<=alpha & LDA_Score>=lda_cutoff)
  if(top>0) plottab <- plottab %>% dplyr::slice_max(LDA_Score,n = top)

  ranknames <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
  names(ranknames) <- c('k','p','c','o','f','g','s')
  plottab$Rank <- sapply(plottab$Feature, function(ft) ranknames[substr(ft,1,1)])
  plottab$Rank <- factor(plottab$Rank, levels=unname(ranknames))
  if(byrank) {
    taxanames <- plottab$Feature
    plottab$Feature <- gsub('.*__','',as.character(plottab$Feature))
  }

  p <- ggpubr::ggdotchart(plottab,x='Feature',y='LDA_Score',color='Enriched',size=6,
                          group='Enriched',
                          add = 'segment',
                          add.params = list(size=2),
                          sorting = 'descending',
                          rotate = T)+
    scale_color_manual(values=colors)+
    labs(y='LDA Score')+
    theme(axis.title = element_text(size=20),
          axis.text = element_text(size=18),
          legend.text = element_text(size = 20),
          legend.title = element_blank())

  suffix <- paste0('_alpha_',alpha,'_lda_',lda_cutoff)
  if(byrank) {
    p <- p+facet_grid(rows = vars(Rank),scales = 'free_y',space = 'free_y')+
      theme(strip.text.y = element_text(size=18))
    suffix <- paste0(suffix,'_byrank')
  }

  show(p)

  saveResults(dataset,foldername = 'LEfSe Dot Chart',
              factors = dataset$factors,
              active_factor = factor$name,
              stat_results = lefse_results,
              width = 12, height = 10,
              suffix = suffix)

  activate(dataset)

  return(p)
}
