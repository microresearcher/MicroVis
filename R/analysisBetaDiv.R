#' Beta Diversity Analysis
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param factor Factor along which to conduct beta-diversity analysis
#' @param stratifier (Optional) A second factor along which to stratify groups
#' @param method Method for dissimilarity calculation. One of either "bray",
#'     "euclidean", "jaccard", "unifrac", "manhattan", "canberra", "clark",
#'     "kulczynski", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#'     "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord"
#' @param weighted (Optional) If method is set to "unifrac", whether to perform
#'     weighted or unweighted unifrac. Defaults to FALSE
#'
#' @return List of data tables containing various beta-diversity results
#' @export
#'
bdiv <- function(dataset=NULL, factor=NULL, stratifier=NULL, method='bray', weighted=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  dataset <- transformData(dataset,temp=T,silent=T,transform_method = 'none')

  rank <- dataset$data$proc$active_rank
  abd <- dataset$data$proc[[rank]]
  abd$Other <- NULL
  metadata <- dataset$metadata
  if(is.null(factor)) factor <- setFVar(dataset)
  if(!is.null(stratifier)) stratifier <- setFVar(dataset,factor_name = stratifier)

  abd$sample <- rownames(abd)
  abd.tmp <- cleanData(merge(metadata,abd),factor)

  if(!is.null(stratifier)) abd.tmp <- cleanData(abd.tmp,stratifier)

  abd <- abd.tmp[(ncol(metadata)+1):ncol(abd.tmp)]
  rownames(abd) <- abd.tmp$sample
  sample_names <- rownames(abd)

  results <- list()
  if(tolower(method)=='euclidean' & mga) {
    PCoA <- prcomp(abd, center = F, scale = F)
    colnames(PCoA$x) <- c(paste0('Axis.',1:ncol(PCoA$x)))
    coord_tab <- data.frame(sample=sample_names,
                            PCoA$x)
  } else {
    if(tolower(method)=='unifrac') {
      dst <- mvunifrac(dataset, weighted = weighted, normalized=F)
    } else {
      dst <- vegdist(abd, method=tolower(method))
    }

    attributes(dst)$Labels <- sample_names
    PCoA <- pcoa(dst)
    coord_tab <- data.frame(sample=sample_names,
                            PCoA$vectors)

    grouping1 <- metadata[metadata$sample %in% sample_names,][[factor$name]]
    if(is.null(stratifier)) dst_stats <- adonis(dst ~ grouping1)
    else {
      grouping2 <- metadata[metadata$sample %in% sample_names,][[stratifier$name]]
      dst_stats <- adonis(dst ~ grouping1 + grouping2)
    }

    # show(dst_stats)
    results$stats <- dst_stats

    results$eigbars <- barplot(PCoA$values$Relative_eig[1:10])
    results$biplot <- biplot.pcoa(PCoA,abd)
  }

  results$coord_tab <- coord_tab

  return(results)
}

#' Unifrac Distance Calculation
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param weighted Whether to perform weighted or unweighted unifrac.
#'     Defaults to FALSE
#' @param normalized Normalize dissimilarities to a range of 0 to 1? Defaults to
#'     TRUE
#'
#' @return "dist" object of unifrac dissimilarity values between samples
#' @export
#'
mvunifrac <- function(dataset=NULL,weighted=F,normalized=T) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  ps <- makePS(dataset)

  unifrac_dist <- UniFrac(ps,weighted=weighted,normalized=normalized)

  return(unifrac_dist)
}
