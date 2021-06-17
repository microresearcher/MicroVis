#' Cluster samples
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param rank Rank of features to use for similarity calculation
#' @param dist_method Method for distance calculation. One of either "bray",
#'     "euclidean", "jaccard", "unifrac", "spearman", "pearson", "kendall",
#'     "manhattan", "canberra", "clark", "kulczynski", "gower", "altGower",
#'     "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao",
#'     "mahalanobis", "chisq" or "chord"
#' @param weighted If using unifrac distance method, whether to perform weighted
#'     or unweighted unifrac. Defaults to FALSE
#' @param clust_method Method for sample clustering. One of either "ward.D",
#'     "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#'     "median" (= WPGMC) or "centroid" (= UPGMC). Defaults to "ward.D2"
#' @param clust_num Number of clusters to make. Defaults to 2
#' @param dataset_name (Not recommended) Name of the dataset to save clusters
#'     to. This should not need to be used by users since the function can
#'     determine the name of the dataset directly passed to it, but not when
#'     it is called within another function.
#'
#' @return List containing the sample adjacency matrix, metadata with an
#'     additional column assigning each sample to a cluster, and the relative
#'     abundance table (with metadata) used for adjacency calculations
#' @export
#'
clusterSamples <- function(dataset=NULL,
                           rank=NULL,
                           dist_method='bray',
                           weighted=F,
                           clust_method='ward.D2',
                           clust_num=2,
                           dataset_name=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  allfts <- getFeatures(dataset, ranks=rank)

  dataset_rel <- scaleSamples(clearNormalization(dataset,temp = T,silent = T),
                              scaling = 'relative', temp = T, silent = T)

  melted <- mvmelt(dataset_rel, rank=rank)

  abun_data <- melted[allfts]
  rownames(abun_data) <- melted$sample

  if(dist_method %in% c('pearson','spearman')) {
    dst <- as.data.frame(cor(t(abun_data), method = dist_method))
  } else {
    dst <- 1-as.data.frame(as.matrix(vegdist(abun_data, method = dist_method)))
  }

  clust <- hclust(as.dist(1-dst), method = clust_method)

  if(clust_num > (2*length(melted$sample))) clust_num <- 2

  clusters <- cutree(clust, k=clust_num)
  sample_names <- dataset$metadata$sample
  for(s in sample_names[!(sample_names %in% names(clusters))]) {
    clusters[[as.character(s)]] <- NA
  }

  clustered_samples <- data.frame(sample=names(clusters),
                                  cluster=clusters)
  dataset$metadata$cluster <- NULL
  metadata_clustered <- merge(dataset$metadata, clustered_samples, 'sample')

  dataset$metadata <- metadata_clustered
  if(dataset_name=='active_dataset') assign(dataset_name,dataset,envir = mvEnv)
  else assign(dataset_name,dataset,1)

  return(list(dst=dst,clusters=metadata_clustered,data=melted))
}
