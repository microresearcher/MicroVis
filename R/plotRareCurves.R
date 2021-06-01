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
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  if(!is.null(dataset$data$proc$normalization)) dataset <- clearNormalization(dataset, temp=T, silent=T)

  rank <- dataset$data$proc$active_rank
  abd <- dataset$data$proc[[rank]]

  clrs <- dataset$colors
  metadata <- dataset$metadata
  cmpgrp <- setFVar(dataset)

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
  maxtab <- cleanData(merge(metadata, maxpts), cmpgrp)
  rctab <- cleanData(merge(metadata, alldata), cmpgrp)

  p <- ggplot(rctab, aes(sample_size,richness,group=sample))+
    geom_line(aes(color=get(cmpgrp$name)),size=0.5)+
    scale_color_manual(values=clrs)+
    labs(y='Genus Richness',x='Sample Size',colour=cmpgrp$txt)
  if(panels) {
    p <- facet(p,cmpgrp$name)
  }
  if(labelSamples) {
    p<-p+geom_label_repel(data=maxtab,
                          aes(label=sample))
  }

  if(getPlot) return(p)
  else show(p)

  saveResults(dataset$results_path,foldername = 'Rarefaction Curves',
              factors = dataset$factors,
              active_factor = cmpgrp$name,
              suffix = '_rc')

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))
  return(dataset)
}
