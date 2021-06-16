#' Plot DESeq Results
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor along which to analyze data
#' @param rank Rank at which to select features
#' @param byGroup Group features by group? Defaults to TRUE
#' @param alpha Significance threshold. Defaults to 0.05
#'
#' @return MicroVis dataset containing DESeq results
#' @export
#'
plotDeseq <- function(dataset=NULL,
                      factor=NULL,
                      rank=NULL,
                      byGroup=T,
                      alpha=0.05) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor
  factor <- setFVar(dataset, factor_name=factor)
  colors <- dataset$colors
  colors <- colors[names(colors) %in% factor$subset]

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  if(is.null(dataset$stats[[factor$name]]$deseq[[rank]])) {
    dataset <- mvdeseq(dataset=dataset,
                       factor=factor$name,
                       rank=rank,
                       dataset_name=dataset_name)
  }
  deseq_res <- dataset$stats[[factor$name]]$deseq[[rank]]

  plottab.all <- deseq_res[!is.na(deseq_res$padj) & deseq_res$padj<=alpha,]
  ftorder <- aggregate(log2FoldChange~Feature,data=plottab.all,mean)
  ftorder <- ftorder[order(ftorder$log2FoldChange),]

  plottab.all$label_ynudge <- sapply(plottab.all$Feature, function(ft) {
    vals <- plottab.all$log2FoldChange[plottab.all$Feature==ft]
    if(all(vals>=0)) {
      -0.5
    } else if(all(vals<0)) {
      0.5
    } else {
      max_abs <- vals[grep(max(abs(vals)),abs(vals))]
      max_abs_side <- max_abs/abs(max_abs)

      max(vals[vals!=max_abs])-max_abs_side*0.5
    }
  })

  temp <- plottab.all
  temp$count <- 1
  temp <- temp %>% pivot_wider(id_cols=c('Reference','Feature'),
                               names_from='Contrast',
                               values_from='count')
  unique_fts <- temp$Feature[complete.cases(temp)]
  other_fts <- temp$Feature[!complete.cases(temp)]

  if(length(unique_fts)) {
    plottab.unique <- plottab.all[plottab.all$Feature %in% unique_fts,]

    p_uniques <- ggbarplot(plottab.unique, x='Feature', y='log2FoldChange',
                           fill='Contrast', color='white',
                           position = position_dodge(), order=ftorder$Feature)+
      geom_hline(yintercept=0, linetype='solid', color='darkgray', size=1)+
      scale_fill_manual(values=colors)+
      labs(y=paste0('Log2FC vs ', plottab.unique$Reference))+
      theme(plot.title = element_text(size=25,hjust = 0.5),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_text(size=20),
            legend.title = element_blank(),
            legend.text = element_text(size=20, margin = margin(r=20)),
            legend.key.size = unit(1,'cm'))+
      coord_flip()+
      geom_text(aes(y=0, label=Feature, size=6,
                    hjust=as.numeric(-plottab.all$label_ynudge/abs(plottab.all$label_ynudge)+1)/2),
                nudge_y=plottab.unique$label_ynudge, show.legend = F,
                check_overlap = T)

    show(p_uniques)
    savedirectory <- saveResults(dataset$results_path,foldername = 'DESeq2',
                                 factors = dataset$factors,
                                 active_factor = factor$name,
                                 suffix = '_uniquesigs',
                                 width = 15, height = 10)
  }

  if(length(other_fts)) {
    plottab.other <- plottab.all[plottab.all$Feature %in% other_fts,]

    p_others <- ggbarplot(plottab.other, x='Feature', y='log2FoldChange',
                         fill='Contrast', color='white',
                         position = position_dodge(), order=ftorder$Feature)+
      geom_hline(yintercept=0, linetype='solid', color='darkgray', size=1)+
      scale_fill_manual(values=colors)+
      labs(y=paste0('Log2FC vs ', plottab.other$Reference))+
      theme(plot.title = element_text(size=25,hjust = 0.5),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_text(size=20),
            legend.title = element_blank(),
            legend.text = element_text(size=20, margin = margin(r=20)),
            legend.key.size = unit(1,'cm'))+
      coord_flip()+
      facet_grid(rows=dplyr::vars(Contrast), scales='free_y', space='free')+
      geom_text(aes(y=0, label=Feature, size=6,
                    hjust=as.numeric(-plottab.all$label_ynudge/abs(plottab.all$label_ynudge)+1)/2),
                nudge_y=plottab.other$label_ynudge, show.legend = F,
                check_overlap = T)+
      theme(strip.text=element_text(size=20))

    show(p_others)
    savedirectory <- saveResults(dataset$results_path,foldername = 'DESeq2',
                                 factors = dataset$factors,
                                 active_factor = factor$name,
                                 suffix = '_othersigs',
                                 width = 15, height = 10)
  }

  if(!is.null(savedirectory)) write.csv(x=deseq_res,
                                        file=file.path(savedirectory,'DESeq2 Results.csv'),
                                        row.names=F)

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))

  return(dataset)
}
