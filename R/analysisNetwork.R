#' Construct an association network between features
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor (Optional) Factor to split data into groups
#' @param rank Rank of features to analyze. Defaults to the active rank
#' @param fts (Optional) Vector of feature names to analyze. Defaults to all the
#'     features in the filtered dataset
#' @param fill (Optional) How to choose the fill color of the nodes. Defaults to
#'     the highest rank for taxonomic data, or a single color for other data.
#'
#' @return MicroVis dataset with network coefficients attached
#' @export
#'
mvNetwork <- function(dataset=NULL, factor=NULL, rank=NULL, fts=NULL, fill=NULL) {
  library(SpiecEasi)

  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  rank <- rank[rank %in% getRanks(dataset)]
  if(!length(rank)) rank <- getLowestRank(dataset)

  fts <- fts[fts %in% getFeatures(dataset)]
  if(!length(fts)) fts <- getFeatures(dataset)

  data <- getdata(selectFeatures(clearNormalization(dataset, temp = T, silent = T),
                                 features = fts, temp = T, silent = T),
                  metadata = F, rank = rank)

  data$Other <- NULL

  ### Using spiec-easi ###
  # https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/inference-of-microbial-ecological-networks.html
  net.se <- spiec.easi(as.matrix(data), method='mb', icov.select.params=list(rep.num=50))

  optnet <- net.se$refit$stars
  netcoefs <- symBeta(getOptBeta(net.se))

  colnames(netcoefs) <- rownames(netcoefs) <- colnames(optnet) <- rownames(optnet) <- colnames(data)

  dataset$networks[[rank]] <- netcoefs

  assign('active_dataset',dataset,envir = mvEnv)
  if(!is.null(dataset$name)) assign(dataset$name,dataset,1)

  return(dataset)
}
