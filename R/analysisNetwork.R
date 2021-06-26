#' Construct an association network between features
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor (Optional) Factor to split data into groups
#' @param rank Rank of features to analyze. Defaults to the active rank
#' @param fts (Optional) Vector of feature names to analyze. Defaults to all the
#'     features in the filtered dataset
#' @param method Co-occurence analysis method
#'
#' @return MicroVis dataset with network coefficients attached
#' @export
#'
mvNetwork <- function(dataset=NULL,
                      factor=NULL,
                      rank=NULL, fts=NULL,
                      method=c('spieceasi','sparcc')) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- factor[factor %in% names(dataset$factors)]

  rank <- rank[rank %in% getRanks(dataset)]
  if(!length(rank)) rank <- dataset$data$proc$active_rank

  fts <- fts[fts %in% getFeatures(dataset, ranks = rank)]
  if(!length(fts)) fts <- getFeatures(dataset, ranks = rank)

  data <- getdata(selectFeatures(clearNormalization(dataset, temp = T, silent = T),
                                 features = fts, ranks = rank, temp = T, silent = T),
                  metadata = F, rank = rank)

  data$Other <- NULL

  method <- match.arg(method)
  if(method=='spieceasi' & requireNamespace('SpiecEasi',quietly = T)) {
    library(SpiecEasi)
    cat('\nBuilding network with SpiecEasi\n')
    ### Using spiec-easi ###
    # https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html
    net.se <- SpiecEasi::spiec.easi(as.matrix(data), method='mb', icov.select.params=list(rep.num=50))

    optnet <- net.se$refit$stars
    netcoefs <- SpiecEasi::symBeta(getOptBeta(net.se))

    colnames(netcoefs) <- rownames(netcoefs) <- colnames(optnet) <- rownames(optnet) <- colnames(data)
  } else {
    cat('\nBuilding network with SparCC\n')
    netcoefs <- sparcc(as.matrix(data))$Cov

    colnames(netcoefs) <- rownames(netcoefs) <- colnames(data)
  }

  dataset$networks[[rank]] <- netcoefs

  assign('active_dataset',dataset,envir = mvEnv)
  if(!is.null(dataset$name)) assign(dataset$name,dataset,1)

  return(dataset)
}

#' Get the co-occurence network for a dataset
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor (Optional) Factor to split data into groups
#' @param rank Rank of features to analyze. Defaults to the active rank
#' @param fts (Optional) Vector of feature names to analyze. Defaults to all the
#'     features in the filtered dataset
#' @param method Co-occurence analysis method
#' @param format (Optional) Format of the network. One of "none", "igraph", or
#'     "network"
#'
#' @return Co-occurence network for a dataset
#' @export
#'
getNetwork <- function(dataset=NULL,
                       factor=NULL,
                       rank=NULL, fts=NULL,
                       method=c('spieceasi','sparcc'),
                       format=c('none','network','igraph')) {

  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- factor[factor %in% names(dataset$factors)]

  rank <- rank[rank %in% getRanks(dataset)]
  if(!length(rank)) rank <- dataset$data$proc$active_rank

  fts <- fts[fts %in% getFeatures(dataset, ranks = rank)]
  if(!length(fts)) fts <- getFeatures(dataset, ranks = rank)

  method <- match.arg(method)
  if(is.null(dataset$networks[[rank]])) dataset <- mvNetwork(dataset,
                                                             factor=factor,
                                                             rank=rank, fts=fts,
                                                             method=method)

  format <- match.arg(format)

  if(format=='igraph') net <- graph_from_adjacency_matrix(dataset$networks[[rank]],
                                                          weighted = T)
  else net <- dataset$networks[[rank]]

  return(net)
}

#' Get feature clusters determined by network analysis
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor (Optional) Factor to split data into groups
#' @param rank Rank of features to analyze. Defaults to the active rank
#' @param fts (Optional) Vector of feature names to analyze. Defaults to all the
#'     features in the filtered dataset
#' @param method Co-occurence analysis method
#'
#' @return
#' @export
#'
#' @examples
networkClusters <- function(dataset=NULL,
                            factor=NULL,
                            rank=NULL, fts=NULL,
                            method=c('spieceasi','sparcc')) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  factor <- factor[factor %in% names(dataset$factors)]

  rank <- rank[rank %in% getRanks(dataset)]
  if(!length(rank)) rank <- dataset$data$proc$active_rank

  fts <- fts[fts %in% getFeatures(dataset, ranks = rank)]
  if(!length(fts)) fts <- getFeatures(dataset, ranks = rank)

  method <- match.arg(method)
  net <- getNetwork(dataset,
                    factor=factor,
                    rank=rank, fts=fts,
                    method=method,
                    format='igraph')

  clusters <- split(names(V(net)), components(net)$membership)

  return(clusters)
}
