#' Plot feature correlation matrix
#'
#' @param dataset1 MicroVis dataset. Defaults to the active dataset
#' @param dataset2 (Optional) MicroVis dataset with features to correlate with
#'     the dataset1
#' @param method Correlation method. One of either "pearson", "spearman", or
#'     "kendall". Defaults to "spearman"
#' @param factor Factor to group samples by
#' @param rank Rank of features to use for calculation
#' @param sigsOnly Select only significant correlations? Defaults to TRUE
#' @param alpha Significance threshold. Defaults to 0.05
#' @param adjustp Use the adjusted p-values? Defaults to TRUE
#' @param r_cutoff Only select correlations with an absolute R-value above a certain
#'     threshold? This value must be between -1 and 1. Defaults to 0
#' @param showlowrs (If r_cutoff is non-zero) Still shade correlations that are
#'     not above the correlation cutoff (r_cutoff)? Defaults to F
#' @param features1 (Optional) Specific features to consider from first dataset.
#'     If only one dataset is being used, then these features will be correlated
#'     with features of the same dataset
#' @param features2 (Optional) Specific features to consider from second dataset.
#' @param matchFts Show the same features in the correlation matrices for all groups?
#'     This can be helpful when comparing correlations between groups. Defaults
#'     to FALSE
#' @param circles Use circles instead of squares to make the matrix? Defaults to
#'     FALSE
#'
#' @return List of the one or two datasets used for this plot
#' @export
#'
plotFtCor <- function(dataset1=NULL, dataset2=NULL,
                      method=c('spearman','pearson','kendall'),
                      factor=NULL,
                      rank=NULL,
                      sigsOnly=T, alpha=0.05, adjustp=T,
                      r_cutoff=0, showlowrs=F,
                      features1=NULL, features2=NULL,
                      matchFts=F,
                      circles=F) {
  if(is.null(dataset)) dataset1 <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset1$name)) dataset1_name <- 'active_dataset'
  else dataset1_name <- dataset1$name

  if(!is.null(dataset2)) {
    dataset2_name <- dataset2$name
    filename <- paste(dataset1_name,dataset2_name,collapse = '_and_')
  } else filename <- dataset1$data$proc$active_rank

  if(length(method) > 1) method <- 'spearman'
  # Function will accept a positive or negative number
  #   as long as its absolute value is between 0 and 1
  r_cutoff <- abs(r_cutoff)
  if(!(r_cutoff>0 & r_cutoff<1)) r_cutoff <- 0
  if(circles) shape <- 'circle'
  else shape <- 'square'

  cordata <- ftcor(dataset1=dataset1,dataset2=dataset2,
                   method=method,
                   features1=features1,features2=features2)

  cordata <- sliceCormats(cordata,
                          sigsOnly=sigsOnly, alpha=alpha, adjustp=adjustp,
                          r_cutoff=r_cutoff,
                          matchFts=matchFts)

  paramtxt <- ''
  if(length(features1) | length(features2)) paramtxt <- paste0('_specificfts')
  if(sigsOnly) paramtxt <- paste0(paramtxt,'_alpha_',alpha)
  if(r_cutoff!=0) paramtxt <- paste0(paramtxt,'_r2cutoff_',r_cutoff)

  for(grp in names(cordata)) {
    grp_txt <- gsub('_',' ',grp)
    cormat <- cordata[[grp]]$r
    if(adjustp) pmat <- cordata[[grp]]$q
    else pmat <- cordata[[grp]]$P

    if(sigsOnly) cormat[pmat>alpha] <- 0
    if(!showlowrs) cormat[abs(cormat)<r_cutoff] <- 0

    p <- ggcorrplot(cormat,
                    title=grp_txt,legend.title='R',
                    outline.color='#383838',
                    method=shape)
    p <- p+theme(plot.title = element_text(size=25,hjust = 0.5),
                 axis.text = element_text(size=20),
                 axis.text.x = element_text(angle = 75))

    if(!exists('save_one_all',inherits = F)) save_one_all <- NULL
    save_one_all <- multisave(save_one_all)

    if(save_one_all %in% c('Yes','Yes to all figures')) saveFig <- T
    else saveFig <- F

    show(p)
    if(saveFig) {
      save_directory <- saveResults(dataset1$results_path,foldername = 'Correlations',
                                    filename = paste0(filename,paramtxt,'_',grp_txt),
                                    width = 7, height = 10,
                                    forcesave = T,
                                    verbose = F)
    }
  }
  if(exists('save_directory')) cat('Figures and any associated statistics saved to:\n ',save_directory)

  return(list(dataset1,dataset2))
}
