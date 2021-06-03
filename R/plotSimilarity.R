#' Plot similarity matrix between samples
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param dist_method Dissimilarity calculation method. One of either "pearson",
#'     "spearman", bray", "euclidean", "jaccard", "unifrac", "manhattan", "canberra",
#'     "clark", "kulczynski", "gower", "altGower", "morisita", "horn", "mountford",
#'     "raup", "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord".
#'     Defaults to "bray"
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
                           dist_method='bray',
                           clust_method='ward.D2',
                           clust_num=2,
                           r_cutoff=0) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  if(r_cutoff==0 | !(abs(r_cutoff)<1)) {
    r_cutoff <- 0
    clrs <- c('red','white','blue')
    clrvals <- c(1,0,-1)
  } else {
    r_cutoff <- abs(r_cutoff)
    clrs <- c('red','white','white','white','blue')
    clrvals <- c(1,r_cutoff,0,-r_cutoff,-1)
  }

  dataset_rel <- scaleSamples(clearNormalization(dataset,temp = T,silent = T),
                              scaling = 'relative', temp = T, silent = T)

  rank <- dataset_rel$data$proc$active_rank
  abd <- dataset_rel$data$proc[[rank]]

  factor <- setFVar(dataset_rel)
  metadata <- dataset_rel$metadata

  if(is.null(abd)) {return()}
  allfts <- colnames(abd)

  # Merge metadata with abundance data
  abd$sample <- rownames(abd)
  abd_data <- cleanData(merge(metadata, abd), factor)
  abd_data <- arrange(abd_data, get(factor$name))

  clust_data <- abd_data[allfts]

  rownames(clust_data) <- abd_data$sample

  if(dist_method %in% c('pearson','spearman')) {
    dst <- as.data.frame(cor(t(clust_data), method = dist_method))
  } else {
    dst <- 1-as.data.frame(as.matrix(vegdist(clust_data, method = dist_method)))
  }

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

  ha <- HeatmapAnnotation(df = abd_data[names(ha_coloring)],
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

  save_directory <- saveResults(dataset$results_path,foldername = 'Similarity Matrices',
                                factors = dataset$factors,
                                figure = hm,
                                active_factor = factor$name,
                                suffix = paste0('_sample-similarity_',
                                                dist_method,'_',clust_method))

  clust <- hclust(as.dist(1-dst), method = clust_method)

  if(clust_num > (2*length(abd_data$sample))) clust_num <- 2

  clusters <- cutree(clust, k=clust_num)
  sample_names <- dataset$metadata$sample
  for(s in sample_names[!(sample_names %in% names(clusters))]) {
    clusters[[as.character(s)]] <- NA
  }

  clustered_samples <- data.frame(sample=names(clusters),
                                  cluster=clusters)
  metadata_clusters <- merge(dataset$metadata, clustered_samples, 'sample')
  metadata_clusters$cluster <- metadata$cluster.y
  metadata_clusters$cluster.x <- NULL
  metadata_clusters$cluster.y <- NULL

  dataset$metadata <- metadata_clusters
  if(dataset_name=='active_dataset') assign(dataset_name,dataset,envir = mvEnv)
  else assign(dataset_name,dataset,1)

  # res_path <- file.path(dataset$results_path,paste0('Results_',Sys.Date()))
  # dir.create(res_path,recursive = T,showWarnings = F)
  if(!is.null(save_directory)) {
    write.csv(metadata,file=file.path(save_directory,'clustered_metadata.csv'))
  }

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))
  print(dataset)

  return(hm)
}
