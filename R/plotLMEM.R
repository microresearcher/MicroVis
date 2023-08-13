#' Plot Volcano Plots of LMEM Log2FC Results
#'
#' @param dataset MicroVis dataset. Defaults to active dataset.
#' @param formula Formula for linear model. Defaults to simple linear model of all covariates.
#' @param exclude Factors/covariates to exclude in linear model.
#' @param rank Rank at which to conduct analysis.
#' @param zero.handling (From linda function) A character string of 'pseudo-count' or 'imputation' indicating the zero handling method used when feature.dat is 'count'. If 'pseudo-count', apseudo.cnt will be added to each value in feature.dat. If 'imputation', then we use the imputation approach using the formula in the referenced paper. Basically, zeros are imputed with values proportional to the sequencing depth. When feature.dat is 'proportion', this parameter will be ignored and zeros will be imputed by half of the minimum for each feature.
#' @param alpha Significance threshold. Defaults to 0.05
#' @param fc Fold-change threshold.
#' @param is.winsor Whether to replace outliers (using winsorization)
#'
#' @return Volcano plots of log2fc values for each feature for groups in each covariate
#' @export
#'
plotLMEM <- function(dataset=NULL,
                     formula=NULL,
                     exclude=NULL,
                     rank=NULL,
                     zero.handling='pseudo-count',
                     alpha=0.05,
                     fc=1.5,
                     is.winsor=T) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(!is.null(formula)) formula <- checkFormula(dataset=dataset, formula=formula)
  if(is.null(formula)) {
    message('Creating basic formula from factors in dataset:')
    vars <- names(dataset$factors)[!(names(dataset$factors) %in% exclude)]
    vars <- vars[sapply(vars, function(x) length(dataset$factors[[x]]$subset) >= 2)]
    formula <- paste0('~',paste(vars, collapse = '+'))
  } else vars <- getFormulaVars(formula)

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  if(is.null(dataset$stats$LMEM[[rank]][[formula]])) {
    dataset <- mvLMEM(dataset = dataset,
                      formula = formula,
                      exclude = exclude,
                      rank = rank,
                      zero.handling = zero.handling,
                      alpha = alpha,
                      is.winsor = is.winsor)
  }

  res <- dataset$stats$LMEM[[rank]][[formula]]

  colors <- dataset$colors
  colors <- colors[names(colors) %in% unlist(lapply(res, function(x) unique(x$Contrast)))]
  colors <- c(colors,'Not diff'='gray')

  p <- list()
  for(var in names(res)) {
    res[[var]]$Diff <- ifelse(abs(res[[var]]$log2FoldChange)>=log2(fc) & res[[var]]$padj<alpha,res[[var]]$Contrast,'Not Diff')
    res[[var]]$FeatureLabel <- ifelse(res[[var]]$Diff=='Not Diff','',res[[var]]$Feature)
    p[[var]] <- ggplot(data = res[[var]], aes(x=log2FoldChange, y=-log10(padj),
                                              shape=Contrast, color=Diff, size=-log10(padj), label=FeatureLabel))+
      geom_point(alpha = 0.5)+
      ggpubr::theme_pubr()+
      ggrepel::geom_text_repel(alpha=1, show.legend=F)+
      scale_color_manual(values=colors)+
      geom_vline(xintercept=c(-log2(fc), log2(fc)), col="black", alpha=0.5)+
      geom_hline(yintercept=-log10(alpha), col="black", alpha=0.5)+
      labs(title=var,
           x=paste0('Log2 Fold Change from "',unique(res[[var]]$Reference),'"'),
           y='-Log10(padj)',
           colour=var)+
      guides(size='none', shape='none',
             color=guide_legend(override.aes=list(size=10)))+
      theme(plot.title = element_text(hjust = 0.5, size=25),
            axis.title.y = element_text(size=20),
            axis.text.y = element_text(size=18),
            axis.title.x = element_text(size=20),
            axis.text.x = element_text(size=18),
            legend.title = element_blank(),
            legend.text = element_text(size=18)
      )
    if(all(res[[var]]$Diff=='Not Diff')) p[[var]] <- p[[var]]+guides(shape=guide_legend(override.aes=list(size=10)))

    show(p[[var]])

    analysis_name <- paste0(var,'__',gsub('~','',formula),'__')
    for(f in dataset$factors) {
      excluded_grps <- f$groups[!(f$groups %in% f$subset)]
      if(!length(excluded_grps)) next
      if(length(excluded_grps) < length(f$subset)) {
        # If less than half of the groups are excluded
        temp <- paste0('No ',paste(excluded_grps, collapse = '-'))
      } else {
        # If less than half of the groups are included
        temp <- paste0('Only ',paste(f$subset, collapse = '-'))
      }
      analysis_name <- paste0(analysis_name,paste0('_',temp))
    }
    savedirectory <- saveResults(dataset,foldername = 'LMEM',
                                 analysis_name = analysis_name,
                                 suffix = paste0('_',rank),
                                 width = 10, height = 7)
  }

  activate(dataset)
}
