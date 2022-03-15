#' Plot rarefaction curves
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param step_size Step size for drawing the curves
#' @param raw Perform on unfiltered data?
#' @param panels Split groups into separate panels?
#' @param labelSamples Label each curve with the sample it corresponds to?
#' @param getPlot Return the plot object?
#'
#' @return Either the MicroVis dataset or the rarefaction curve plot object
#' @export
#'
plotRareCurves <- function(dataset=NULL,
                           step_size=100,
                           raw=F,
                           panels=TRUE,
                           labelSamples=FALSE,
                           getPlot=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  if(!is.null(dataset$data$proc$normalization)) dataset <- clearNormalization(dataset, temp=T, silent=T)

  rank <- dataset$data$proc$active_rank
  abd <- dataset$data$proc[[rank]]

  metadata <- dataset$metadata
  factor <- setFVar(dataset)
  colors <- dataset$colors
  colors <- colors[names(colors) %in% factor$subset]

  maxrich <- max(rowSums(abd))
  minrich <- min(rowSums(abd))

  p <- rarecurve(abd, step = step_size)
  names(p) <- rownames(abd)

  alldata <- c()
  for(i in 1:length(p)) {
    sname <- names(p)[i]
    sdata <- data.frame(sample_size=as.numeric(str_replace(names(p[[i]]),'N','')),
                        richness=p[[i]])
    sdata <- data.frame(sample=rep(sname,nrow(sdata)),sdata)

    alldata <- rbind(alldata, sdata)
  }
  maxpts <- aggregate(sample_size ~ sample, data=alldata, max)
  maxpts <- merge(maxpts, alldata)

  # Merge these two data tables with the metadata
  maxtab <- cleanData(merge(metadata, maxpts), factor)
  rctab <- cleanData(merge(metadata, alldata), factor)

  p <- ggplot(rctab, aes(sample_size,richness,group=sample))+
    geom_line(aes(color=get(factor$name)),size=0.5)+
    scale_color_manual(values=colors)+
    labs(y='Genus Richness',x='Sample Size',colour=factor$txt)
  if(panels) {
    p <- facet(p,factor$name)
  }
  if(labelSamples) {
    p<-p+geom_label_repel(data=maxtab,
                          aes(label=sample))
  }

  if(!getPlot) {
    show(p)
    saveResults(dataset,foldername = 'Rarefaction Curves',
                factors = dataset$factors,
                active_factor = factor$name,
                suffix = '_rc')

    activate(dataset)
  }

  return(p)
}
