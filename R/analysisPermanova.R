#' Distance Calculation
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param method Method for dissimilarity calculation. One of either "bray",
#'     "euclidean", "jaccard", "unifrac", "manhattan", "canberra", "clark",
#'     "kulczynski", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#'     "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord"
#' @param weighted (Optional) If method is set to "unifrac", whether to perform
#'     weighted or unweighted unifrac. Defaults to FALSE
#' @param allFactors When performing distance calculation, only include samples
#'     in the active subgroups of all factors? Defaults to TRUE. Otherwise, user
#'     can specify which factors using the factor parameter (explained below)
#' @param factors (Optional) Specify factors to filter samples by
#'
#' @return Dissimilarity matrix and other associated details of the calculation
#' @export
#'
mvdist <- function(dataset=NULL, method='bray', weighted=F, allFactors=T, factors=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  dataset <- transData(dataset,temp=T,silent=T,transform_method = 'none')

  rank <- dataset$data$proc$active_rank
  abd <- dataset$data$proc[[rank]]
  abd$Other <- NULL

  factors <- factors[factors %in% names(dataset$factors)]
  if(!length(factors)) {
    if(!allFactors) factors <- select.list(names(dataset$factors),
                                           multiple = T, graphics = T)
    else factors <- names(dataset$factors)
  }

  metadata <- dataset$metadata[c('sample', factors)]

  abd$sample <- rownames(abd)
  abd.tmp <- merge(metadata, abd)

  for(f in factors) abd.tmp <- cleanData(abd.tmp, setFVar(dataset, f))

  abd <- abd.tmp[(ncol(metadata)+1):ncol(abd.tmp)]
  rownames(abd) <- abd.tmp$sample
  sample_names <- rownames(abd)

  dst_results <- list(metadata=abd.tmp[1:ncol(metadata)])
  if(tolower(method)=='unifrac') {
    dst_results$dst <- mvunifrac(dataset, weighted = weighted, normalized=F)
  } else {
    dst_results$dst <- vegdist(abd, method=tolower(method))
  }

  attributes(dst_results$dst)$Labels <- sample_names
  PCoA <- pcoa(dst_results$dst)
  coord_tab <- data.frame(sample=sample_names,
                          PCoA$vectors)

  dst_results$eigbars <- barplot(PCoA$values$Relative_eig[1:10])
  dst_results$biplot <- biplot.pcoa(PCoA,abd)

  return(dst_results)
}

#' PERMANOVA Analysis
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param dist Method for dissimilarity calculation. One of either "bray",
#'     "euclidean", "jaccard", "unifrac", "manhattan", "canberra", "clark",
#'     "kulczynski", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#'     "binomial", "chao", "cao", "mahalanobis", "chisq" or "chord"
#' @param weighted (Optional) If method is set to "unifrac", whether to perform
#'     weighted or unweighted unifrac. Defaults to FALSE
#' @param allFactors Perform PERMANOVA analysis on all factors? Defaults to TRUE
#'     if no factors are specified (below)
#' @param factors (Optional) Specify factors to perform PERMANOVA on. If none are
#'     specified, then all factors will be included
#' @param alpha Significance threshold for PERMANOVA results before performing
#'     PERMDISP analysis. Defaults to 0.05
#'
#' @return PERMANOVA results
#' @export
#'
pnova <- function(dataset=NULL, dist='bray', weighted=F, allFactors=T, factors=NULL, alpha=0.05) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  dataset <- transData(dataset,temp=T,silent=T,transform_method = 'none')

  factors <- factors[factors %in% names(dataset$factors)]
  if(!length(factors)) {
    if(!allFactors) factors <- select.list(names(dataset$factors),
                                           multiple = T, graphics = T)
    else factors <- names(dataset$factors)
  }

  dst_results <- mvdist(dataset, method=dist, weighted=weighted, factors=factors)
  dst <- dst_results$dst

  f <- unlist(lapply(factors, function(x) dst_results$metadata[x]),recursive = F)

  formula <- as.formula(paste('dst ~ ',paste0('f$',names(f),collapse = '+')))

  dst_stats <- list()
  dst_stats$pnova <- adonis(formula)

  dst_stats$summary <- data.frame(Factor='factor',
                                  Permanova=alpha,
                                  Dispersion=alpha,
                                  Significant='')[0,]
  for(grping in names(f)) {
    dst_stats$permdisp[[sub('f\\$','',grping)]] <- anova(betadisper(dst, f[[grping]]))

    pnova_res <- dst_stats$pnova$aov.tab[paste0('f$',grping),'Pr(>F)']
    disp_res <- dst_stats$permdisp[[sub('f\\$','',grping)]]$`Pr(>F)`[[1]]
    overall_res <- ifelse(pnova_res <= alpha & disp_res > alpha, 'Yes', 'No')

    dst_stats$summary <- rbind(dst_stats$summary,
                               data.frame(Factor=grping,
                                          Permanova=pnova_res,
                                          Dispersion=disp_res,
                                          Significant=overall_res))
  }

  return(dst_stats)
}
