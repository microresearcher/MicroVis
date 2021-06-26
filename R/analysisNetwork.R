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

  rank <- rank[rank %in% getRanks(dataset)]
  if(!length(rank)) rank <- dataset$data$proc$active_rank

  fts <- fts[fts %in% getFeatures(dataset, ranks = rank)]
  if(!length(fts)) fts <- getFeatures(dataset, ranks = rank)

  data <- getdata(selectFeatures(clearNormalization(dataset, temp = T, silent = T),
                                 features = fts, ranks = rank, temp = T, silent = T),
                  metadata = F, rank = rank)

  data$Other <- NULL

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
