#' Merge Two Datasets with Overlapping Samples
#'
#' @param dataset1 First MicroVis dataset
#' @param dataset2 Second MicroVis dataset
#' @param dataset1_name (Optional) What to name the first dataset.
#' @param dataset2_name (Optional) What to name the second dataset.
#' @param features1 (Optional) Features to choose from first dataset.
#' @param features2 (Optional) Features to choose from second dataset.
#' @param keepFiltered Whether or not to keep the filtered-out taxa that are lumped
#'     into the "Other" feature by runFeatureFilter().
#' @param sigsOnly1 Whether or not to only consider the significant features of
#'     the first dataset. Default is FALSE.
#' @param sigsOnly2 Whether or not to only consider the significant features of
#'     the second dataset. Default is FALSE.
#' @param alpha Significance threshold. Default is 0.05.
#'
#' @return Merged dataset with the samples that were shared across both datasets
#'     and all the features selected from both datasets.
#' @export
#'
mvmerge <- function(dataset1,dataset2,
                    dataset1_name=NULL,dataset2_name=NULL,
                    features1=NULL,features2=NULL,keepFiltered=F,
                    sigsOnly1=F,sigsOnly2=F,alpha=0.05) {
  # Load both datasets
  dataset1_name <- dataset1$name
  dataset2_name <- dataset2$name
  if(is.null(dataset1_name)) dataset1_name <- readline(cat('Provide a short name for the first dataset:\n'))
  if(is.null(dataset2_name)) dataset2_name <- readline(cat('Provide a short name for the second dataset:\n'))

  ranks1 <- getRanks(dataset1)
  rank1 <- dataset1$data$proc$active_rank
  abun1 <- dataset1$data$proc[[rank1]]
  fts1 <- colnames(abun1)
  metadata1 <- dataset1$metadata
  factors1 <- dataset1$factors
  clrs1 <- dataset1$colors

  ranks2 <- getRanks(dataset2)
  rank2 <- dataset2$data$proc$active_rank
  abun2 <- dataset2$data$proc[[rank2]]
  fts2 <- colnames(abun2)
  metadata2 <- dataset2$metadata
  factors2 <- dataset2$factors
  clrs2 <- dataset2$colors

  if(!keepFiltered) {
    abun1$Other <- NULL
    abun2$Other <- NULL
    fts1 <- fts1[fts1!='Other']
    fts2 <- fts2[fts2!='Other']
  }

  # Select features if any are specified
  if(length(features1 %in% fts1)) fts1 <- fts1[fts1 %in% features1]
  if(sigsOnly1) {
    dataset1 <- univar(dataset=dataset1,
                       rank=rank1,
                       features=fts1,
                       factor=dataset1$active_factor)$overall
    stats1 <- dataset1$stats[[dataset1$active_factor]]$univar[[rank]]$stats
    sigfts1 <- fts1[fts1 %in% stats1$.y.[stats1$p.adj < alpha]]
    if(length(sigfts1)) fts1 <- union(fts1,sigfts1)
    else message('\nNo significant features were found in ',dataset1_name,' that overlapped with already selected features\n')
  }

  if(length(features2 %in% fts2)) fts2 <- fts2[fts2 %in% features2]
  if(sigsOnly2) {
    stats2 <- univar(dataset=dataset2,
                     rank=rank2,
                     features=fts2,
                     factor=dataset2$active_factor)$overall
    stats2 <- dataset2$stats[[dataset2$active_factor]]$univar[[rank]]$stats
    sigfts2 <- fts2[fts2 %in% stats2$.y.[stats2$p.adj < alpha]]
    if(length(sigfts2)) fts2 <- union(fts2,sigfts2)
    else message('\nNo significant features were found in ',dataset2_name,' that overlapped with already selected features\n')
  }

  # If any features share the same name in the two datasets, append all features
  #   with the shorthand of the dataset they came from
  if(any(colnames(abun1) %in% colnames(abun2))) {
    colnames(abun1) <- paste(colnames(abun1),dataset1_name,sep = '_')
    fts1 <- paste(fts1,dataset1_name,sep = '_')

    colnames(abun2) <- paste(colnames(abun2),dataset2_name,sep = '_')
    fts2 <- paste(fts2,dataset2_name,sep = '_')
  }

  merged_abun <- merge(abun1, abun2, by=0)
  shared_samples <- merged_abun$Row.names
  merged_abun$Row.names <- NULL

  if(!length(shared_samples)) return(message('\nERROR: ',
                                             dataset1_name,' and ',dataset2_name,
                                             ' do not share any samples. Check the sample names to make sure they match.\n'))

  lost_samples1 <- rownames(abun1)[!(rownames(abun1) %in% shared_samples)]
  if(length(lost_samples1)) message('\nThe following samples in "',
                                    dataset1_name,'" were not in "',dataset2_name,
                                    '":\n ',paste0(lost_samples1,collapse = '\t'),'\n')

  lost_samples2 <- rownames(abun2)[!(rownames(abun2) %in% shared_samples)]
  if(length(lost_samples2)) message('\nThe following samples in "',
                                    dataset2_name,'" were not in "',dataset1_name,
                                    '":\n ',paste0(lost_samples2,collapse = '\t'),'\n')

  # Turn the sample names into row names
  rownames(merged_abun) <- shared_samples

  # Merge the metadatas
  metadata <- merge(dataset1$metadata, dataset2$metadata)

  # Check to make sure all the shared samples are in the metadata
  if(!all(shared_samples %in% metadata$sample)) message('\nWARNING: Not all samples are in the merged metadata... for some reason')

  factors <- c(factors1,factors2[!(factors2 %in% factors1)])
  clrs <- c(clrs1,clrs2[!(clrs2 %in% clrs1)])

  features <- list(fts1,fts2)
  names(features) <- c(dataset1_name,dataset2_name)

  if(dataset1$features=='taxa') {
    ftnames1 <- dataset1$data$taxa_names
    ftnames1 <- ftnames1[!duplicated(ftnames1[[rank1]]),
                         1:grep(rank1,colnames(ftnames1))]
  }

  if(dataset2$features=='taxa') {
    ftnames2 <- dataset2$data$taxa_names
    ftnames2 <- ftnames2[!duplicated(ftnames2[[rank2]]),
                         1:grep(rank2,colnames(ftnames2))]
  }

  if(exists('ftnames1') & exists('ftnames2')) {
    rownames(ftnames1) <- paste0(dataset1_name,'_',rownames(ftnames1))
    rownames(ftnames2) <- paste0(dataset2_name,'_',rownames(ftnames2))
    ftnames <- bind_rows(ftnames1,ftnames2)
  } else if(exists('ftnames1')) {
    ftnames <- ftnames1
  } else if(exists('ftnames2')) {
    ftnames <- ftanmes2
  } else ftnames <- NULL

  merged_dataset <- list(metadata=metadata,
                         data=list(taxa_names=ftnames,
                                   proc=list(single_rank=merged_abun,
                                             active_rank='single_rank'),
                                   features=features),
                         features='merged',
                         factors=factors,
                         active_factor=dataset1$active_factor,
                         colors=clrs,
                         results_path=dataset1$results_path,
                         name=paste0('mvmerged_',dataset1_name,'_and_',dataset2_name))

  class(merged_dataset) <- 'mvmerged'

  assign(merged_dataset$name, merged_dataset, 1)

  return(merged_dataset)
}
