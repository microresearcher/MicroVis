#' Plot similarity matrix between samples
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param rank Rank of features to use for similarity calculation
#' @param dist_method Dissimilarity calculation method. One of either "pearson",
#'     "spearman", bray", "euclidean", "jaccard", "unifrac", "manhattan", "canberra",
#'     "clark", "kulczynski", "gower", "altGower", "morisita", "horn", "mountford",
#'     "raup", "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord".
#'     Defaults to "bray"
#' @param weighted If performing unifract distance, whether to use weighted or
#'     unweighted unifrac. Defults to FALSE (unweighted)
#' @param clust_method Clustering method. One of either "ward.D", "ward.D2",
#'     "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#'     "median" (= WPGMC) or "centroid" (= UPGMC). Defaults to "ward.D2"
#' @param clust_num Number of clusters to try to make
#' @param r_cutoff R-values with absolute value below this cutoff will be shaded
#'     white. Defaults to 0
#'
#' @return Similarity matrix heatmap
#' @export
#'
plotSimilarity <- function(dataset=NULL,
                           rank=NULL,
                           dist_method='bray',
                           weighted=F,
                           clust_method='ward.D2',
                           clust_num=2,
                           r_cutoff=0) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  if(r_cutoff==0 | !(abs(r_cutoff)<1)) {
    r_cutoff <- 0
    clrs <- c('red','white','blue')
    clrvals <- c(1,0,-1)
  } else {
    r_cutoff <- abs(r_cutoff)
    clrs <- c('red','white','white','white','blue')
    clrvals <- c(1,r_cutoff,0,-r_cutoff,-1)
  }

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  cluster_res <- clusterSamples(dataset,
                                rank=rank,
                                dist_method=dist_method, weighted=weighted,
                                clust_method=clust_method, clust_num=clust_num,
                                dataset_name=dataset_name)

  dst <- cluster_res$dst
  metadata_clustered <- cluster_res$clusters
  abun_data <- cluster_res$data

  if(min(dst)>=0) {
    corr_clrs <- colorRamp2(c(0,r_cutoff,1),c('white','white','red'))
  } else if(max(dst<=0)) {
    corr_clrs <- colorRamp2(c(-1,-r_cutoff,0),c('white','white','red'))
  } else {
    corr_clrs <- colorRamp2(c(-1,-r_cutoff,0,r_cutoff,1),c('white','white','red'))
  }

  # Make a colors list from "clrs" that HeatmapAnnotation() will accept
  ha_coloring <- list()
  factor_clrs <- dataset$colors
  if(length(dataset$factors)) {
    for(f in names(dataset$factors)) {
      if(length(dataset$factors[[f]]$subset)>1) ha_coloring[[f]] <- factor_clrs
    }
  }

  ha <- HeatmapAnnotation(df = abun_data[names(ha_coloring)],
                          col = ha_coloring)

  # The Heatmap function from ComplexHeatmap is expansive and can also do the
  #   distance calculation for us, but we want the option to modify values
  #   afterwards and also record the clusters, so the function does the calculation
  #   on its own
  hm <- Heatmap(as.matrix(dst),
                col = corr_clrs,
                top_annotation = ha,
                clustering_method_columns = clust_method)
  draw(hm)

  save_directory <- saveResults(dataset$results_path, foldername = 'Similarity Matrices',
                                factors = dataset$factors,
                                active_factor = dataset$active_factor,
                                other_results = list(Clusters=metadata_clustered),
                                figure = hm,
                                suffix = paste0('_sample-similarity_',
                                                dist_method,'_',clust_method))

  activate(dataset)

  return(hm)
}
