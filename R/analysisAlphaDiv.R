#' Alpha Diversity Calculation
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param method Method for alpha diversity calculation. One of either "chao1",
#'     "shannon", "simpson", "invsimpson", or "pd". Defaults to "chao1"
#' @param rooted (Optional) If method is set to "pd", calculate rooted or unrooted
#'     Faith's phylogenetic diversity index? Defaults to FALSE
#'
#' @return Table with the alpha diversity index for each sample
#' @export
#'
adiv <- function(dataset=NULL, method='chao1', rooted=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  method <- tolower(method)

  abd <- clearProcessing(dataset, temp = T, silent = T)$data$proc$unranked

  sample_names <- rownames(abd)

  if(method=='chao1') divs <- apply(abd, 1, function(x) sum(x>0))
  else if(method=='pd') divs <- phyloDiversity(dataset,rooted)
  else divs <- vegan::diversity(abd,method)

  div_tab <- data.frame(sample_names,
                        divs)

  colnames(div_tab) <- c('sample',method)
  return(div_tab)
}

#' Faith's Phylogenetic Diversity Index
#'
#' @param dataset MicroVis dataset. Defaults to active dataset
#' @param rooted Calculate rooted or unrooted Faith's phylogenetic diversity index?
#'     Defaults to FALSE
#'
#' @return Single data.frame column of Faith's phylogenetic diversity index values
#'     for each sample
#' @export
#'
phyloDiversity <- function(dataset=NULL, rooted=F) {
  if(!requireNamespace('picante', quietly = T)) {
    stop(paste0('You need to install the "picante" package to do this. This package is required for phylogenic distance analysis.\nInstall this package with the following command:\n  install.packages("picante")'))
  }

  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  dataset <- clearFeatureFilt(dataset, temp=T, silent=T)

  abun <- dataset$data$proc[[getLowestRank(dataset)]]
  colnames(abun) <- TaxatoASV(dataset$data,colnames(abun),getLowestRank(dataset))

  tree <- makePS(dataset)@phy_tree

  if(rooted) tree <- picante::prune.sample(abun,tree)

  div <- picante::pd(abun,tree,include.root = rooted)

  return(div['PD'])
}
