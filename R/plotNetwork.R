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
                        fill=NULL, outline=NULL, labelfts=NULL, labelAll=F,
                        deg_cutoff=0, r_cutoff=0, top_r_prop=100,
                        layout=c('fr','circle','sphere','dh','nicely')) {

  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- factor[factor %in% names(dataset$factors)]

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  fts <- fts[fts %in% getFeatures(dataset, ranks = rank)]
  if(!length(fts)) fts <- getFeatures(dataset, ranks = rank)

  r_cutoff <- abs(r_cutoff)

  method <- match.arg(method)

  if(length(fill)>1) {
    if(length(fill)!=length(fts)) fill <- NULL
    else filllist <- fill
  }

  if(length(fill)<=1) {
    fill <- fill[fill %in% getRanks(dataset)]
    if(!length(fill)) fill <- getRanks(dataset)[2]
    filllist <- getParentTaxa(fts, from=rank, to=fill,
                              dataset$data$taxa_names)
  }

  layout <- match.arg(layout)

  data <- getdata(dataset, rank = rank, metadata = F)
  data$Other <- NULL

  net <- getNetwork(dataset,
                    factor=factor,
                    rank=rank, fts=fts,
                    method=method)

  totaledges <- length(as.matrix(net)[as.matrix(net)!=0])/2
  if(totaledges) cat('\n',totaledges,'significant correlations identified\n')
  else stop('No significant correlations were identified using ',method)

  net <- graph_from_adjacency_matrix(net, weighted = T)

  vsize <- glog(apply(data, 2, mean))
  vsize[vsize<(max(vsize)/5)] <- max(vsize)/5

  E(net)$polarity <- E(net)$weight/abs(E(net)$weight)
  E(net)$weight <- abs(E(net)$weight)
  E(net)$width <- abs(E(net)$weight)
  E(net)$arrow.size <- 0

  if(!labelAll) V(net)$label <- rep(NA, length(V(net)$name))
  if(length(labelfts)) V(net)$label[which(V(net)$name %in% labelfts)] <- labelfts

  V(net)$size <- vsize*20
  V(net)$fill <- filllist

  if(r_cutoff>0 & r_cutoff<1) {
    cutoff <- r_cutoff
    net.filt <- delete_edges(net, E(net)[weight<cutoff])
    cat('  Filtered out',length(E(net)[weight<cutoff]),
        'correlations less than', r_cutoff,'\n')

    emptyvs.filt <- V(net.filt)[degree(net.filt)<=deg_cutoff]
    cat('  Removed nodes with',deg_cutoff,'or fewer correlations\n')

    net.clean <- delete_vertices(net.filt, emptyvs.filt)
  } else if(top_r_prop>0 & top_r_prop<100) {
    cutoff <- min(E(net)$weight[order(E(net)$weight,decreasing = T)]
                  [floor((top_r_prop/100)*length(E(net)$weight))+1])
    net.filt <- delete_edges(net, E(net)[weight<cutoff])
    cat('  Filtered out',length(E(net)[weight<cutoff]),
        'correlations in the bottom',100-top_r_prop,'percentile\n')

    emptyvs.filt <- V(net.filt)[degree(net.filt)<=deg_cutoff]
    cat('  Removed nodes with',deg_cutoff,'or fewer correlations\n')

    net.clean <- delete_vertices(net.filt, emptyvs.filt)
  } else {
    emptyvs <- V(net)[degree(net)<=deg_cutoff]
    net.clean <- delete_vertices(net, emptyvs)
  }

  colors <- rainbow(length(unique(V(net.clean)$fill)), alpha = 0.5)
  names(colors) <- unique(V(net.clean)$fill)

  V(net.clean)$color <- colors[V(net.clean)$fill]
  V(net.clean)$frame.color <- NA

  if(layout=='fr') lay <- layout_with_fr(net.clean)
  else if(layout=='circle') lay <- layout_in_circle(net.clean)
  else if(layout=='sphere') lay <- layout_on_sphere(net.clean)
  else if(layout=='dh') lay <- layout_with_dh(net.clean)
  else if(layout=='nicely') lay <- layout_nicely(net.clean)

  plot(net.clean, layout=lay)
  legend(x=-1.5, y=-0.6, unique(V(net.clean)$fill), pch=21,
         col="#ffffff", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1)

  return(dataset)
}
