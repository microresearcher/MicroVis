#' Plot Feature Network
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor (Optional) Factor to split data along groups
#' @param rank Rank at which to build network
#' @param fts List of specific features to analyze and plot
#' @param method Co-occurence analysis method
#' @param fill How node colors are chosen
#' @param outline (Optional) How node outlines are chosen, if at all. Defaults to
#'     no outline
#' @param labelfts Which nodes are labeled, if any. Defaults to none
#' @param layout Layout for nodes
#'
#' @return Network plot
#' @export
#'
plotNetwork <- function(dataset=NULL,
                        factor=NULL,
                        rank=NULL, fts=NULL,
                        method=c('spieceasi','sparcc'),
                        fill=NULL, outline=NULL, labelFts=NULL, labelAll=F,
                        layout=c('fr','circle')) {

  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- factor[factor %in% names(dataset$factors)]

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  fts <- fts[fts %in% getFeatures(dataset, ranks = rank)]
  if(is.null(fts)) fts <- getFeatures(dataset, ranks = rank)

  if(length(fill)>1) {
    if(length(fill)!=length(fts)) fill <- NULL
    else filllist <- fill
  }

  if(length(fill)<=1) {
    fill <- fill[fill %in% getRanks(dataset)]
    if(!length(fill)) fill <- getRanks(dataset)[2]
    filllist <- getParentTaxa(fts,
                              from=rank,to=fill,
                              dataset$data$taxa_names)
  }

  layout <- match.arg(layout)

  data <- getdata(dataset, rank = rank, metadata = F)
  data$Other <- NULL

  netcoefs <- getNetwork(dataset,
                         factor=factor,
                         rank=rank, fts=fts,
                         method=method,
                         format='igraph')

  netcoefs <- dataset$networks[[rank]]

  totaledges <- length(as.matrix(netcoefs)[as.matrix(netcoefs)!=0])/2
  if(totaledges) cat('\n',totaledges,'significant correlations identified\n')
  else stop('No significant correlations were identified using',method)

  vsize <- glog(apply(data, 2, mean))
  vsize[vsize<(max(vsize)/5)] <- max(vsize)/5

  net.ig <- igraph::graph_from_adjacency_matrix(netcoefs, weighted = T)

  E(net.ig)$polarity <- E(net.ig)$weight/abs(E(net.ig)$weight)
  E(net.ig)$weight <- abs(E(net.ig)$weight)
  E(net.ig)$width <- abs(E(net.ig)$weight)*20
  E(net.ig)$arrow.size <- 0.2

  if(!labelAll) V(net.ig)$label <- rep(NA, length(V(net.ig)$name))
  if(length(labelfts)) V(net.ig)$label[which(V(net.ig)$name %in% labelfts)] <- labelfts

  V(net.ig)$size <- vsize*20
  V(net.ig)$fill <- filllist

  emptyvs <- V(net.ig)[degree(net.ig)==0]
  net.clean <- delete_vertices(net.ig, emptyvs)

  net.sel <- net.clean

  colors <- rainbow(length(unique(V(net.sel)$fill)), alpha = 0.5)
  names(colors) <- unique(V(net.sel)$fill)

  V(net.sel)$color <- colors[V(net.sel)$fill]
  V(net.sel)$frame.color <- NA

  if(layout=='fr') lay <- layout_with_fr(net.sel)
  else if(layout=='circle') lay <- layout_in_circle(net.sel)

  plot(net.sel, layout=lay)
  legend(x=-1.5, y=-0.6, unique(V(net.sel)$fill), pch=21,
         col="#ffffff", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1)

  return(dataset)
}
