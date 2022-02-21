#' Compute feature correlation matrices for each group and all groups together
#'
#' @param dataset1 MicroVis dataset. Defaults to the active dataset
#' @param dataset2 (Optional) MicroVis dataset with features to correlate with
#'     the `dataset1`
#' @param method Correlation method. One of either "pearson", "spearman", or
#'     "kendall". Defaults to "spearman"
#' @param factor Factor to group samples by
#' @param features1 (Optional) Specific features to consider from first dataset.
#'     If only one dataset is being used, then these features will be correlated
#'     with features of the same dataset
#' @param features2 (Optional) Specific features to consider from second dataset.
#'     If no second dataset is provided, then `features1` will be correlated with
#'     `features2` from `dataset1`
#'
#' @return List of correlation matrices for each group with additional matrices
#'     containing p-values and adjusted p-values
#' @export
#'
ftcor <- function(dataset1=NULL,dataset2=NULL,
                  method=c('pearson','kendall','spearman'),
                  factor=NULL,
                  features1=NULL, features2=NULL) {
  if(is.null(dataset1)) dataset1 <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset1$name)) dataset1_name <- 'active_dataset'
  else dataset1_name <- dataset1$name

  if(length(method) > 1) method <- 'spearman'

  if(is.null(dataset2)) {
    clrs <- dataset1$colors
    rank <- dataset1$data$proc$active_rank
    abun_data <- mvmelt(dataset1)
    fts <- colnames(dataset1$data$proc[[rank]])
    fts <- fts[!(fts %in% c('Other','Unknown'))]

    if(length(features1 %in% fts)) fts1 <- fts[fts %in% features1]
    else fts1 <- fts
    if(length(features2 %in% fts)) fts2 <- fts[fts %in% features1]
    else fts2 <- fts

    factor <- dataset1$factors[[dataset1$active_factor]]
  } else {
    dataset2_name <- dataset2$name
    merged_dataset <- mvmerge(dataset1,dataset2,
                              features1=features1,features2=features2)
    abun_data <- mvmelt(merged_dataset)
    # In case any specified features weren't found during merging
    fts1 <- merged_dataset$data$features[[1]]
    fts2 <- merged_dataset$data$features[[2]]
    fts1 <- fts1[!(fts1 %in% c('Other','Unknown'))]
    fts2 <- fts2[!(fts2 %in% c('Other','Unknown'))]
    clrs <- merged_dataset$colors
    factor <- merged_dataset$factors[[merged_dataset$active_factor]]

    filename <- paste(dataset1_name,dataset2_name,collapse = '_and_')
  }
  allfts <- union(fts1,fts2)

  cormats_bygrp <- list()

  cormat <- rcorr(as.matrix(abun_data[allfts]),type=method)
  cormat$n <- NULL
  if(is.null(dataset2)) diag(cormat$r) <- 0
  for(mat in names(cormat)) {
    cormat[[mat]] <- data.frame(cormat[[mat]])
    cormat[[mat]] <- cormat[[mat]][rownames(cormat[[mat]]) %in% fts1,
                                   colnames(cormat[[mat]]) %in% fts2]
  }
  cormats_bygrp[['All_Groups']] <- cormat

  for(grp in factor$subset) {
    abun.temp <- abun_data[abun_data[[factor$name]]==grp,][allfts]
    cormat <- rcorr(as.matrix(abun.temp[allfts],type=method))
    cormat$n <- NULL
    if(is.null(dataset2)) diag(cormat$r) <- 0
    for(mat in names(cormat)) {
      cormat[[mat]] <- data.frame(cormat[[mat]])
      cormat[[mat]] <- cormat[[mat]][rownames(cormat[[mat]]) %in% fts1,
                                     colnames(cormat[[mat]]) %in% fts2]
    }
    cormats_bygrp[[grp]] <- cormat
  }

  # Create a matrix of adjusted p-values for each group
  for(grp in names(cormats_bygrp)) {
    pmat <- cormats_bygrp[[grp]]$P
    rn <- rownames(pmat)
    cn <- colnames(pmat)

    qmat <- data.frame(matrix(p.adjust(as.matrix(pmat),method = 'BH'),
                              nrow = length(rn)))

    rownames(qmat) <- rn
    colnames(qmat) <- cn

    cormats_bygrp[[grp]]$q <- qmat
  }

  return(cormats_bygrp)
}


#' Select only significant and/or high-correlation values from correlation matrices
#'
#' @param cormats List of correlation matrices (output of ftcor)
#' @param sigsOnly Select only significant correlations? Defaults to TRUE
#' @param alpha Significance threshold. Defaults to 0.05
#' @param adjustp Use the adjusted p-values? Defaults to TRUE
#' @param r_cutoff Only select correlations with an absolute R-value above a certain
#'     threshold? This value must be between -1 and 1. Defaults to 0
#' @param matchFts Show the same features in the correlation matrices for all groups?
#'     This can be helpful when comparing correlations between groups. Defaults
#'     to FALSE
#'
#' @return List of correlation matrices with only the desired correlations remaining
#' @export
#'
sliceCormats <- function(cormats,
                         sigsOnly=T, alpha=0.05, adjustp=T,
                         r_cutoff=0, matchFts=F) {
  # Function will accept a positive or negative number
  #   as long as its absolute value is between 0 and 1
  r_cutoff <- abs(r_cutoff)
  if(!(r_cutoff>0 & r_cutoff<1)) r_cutoff <- 0

  fts1 <- list()
  fts2 <- list()
  for(grp in names(cormats)) {
    cormat <- cormats[[grp]]

    fts1[[grp]] <- rownames(cormat$r)
    fts2[[grp]] <- colnames(cormat$r)

    if(sigsOnly) {
      # Get the significant correlations
      if(adjustp) pvals <- cormat$q
      else pvals <- cormat$P
      pvals$Feature1 <- rownames(pvals)
      pvals <- pvals %>% pivot_longer(cols=colnames(cormat$q),names_to='Feature2',values_to='q')
      sigs1 <- pvals[pvals$q<=alpha,]$Feature1
      sigs2 <- pvals[pvals$q<=alpha,]$Feature2

      # Warn user if they wanted significant features but none were found
      if(!length(c(sigs1,sigs2))) message('\nNo correlations found at a significance level of ',alpha,
                                          ' for ',grp)
      else {
        # If there is at least one feature in sigs1, there will be at least one
        #   in sigs2, and vice-versa
        fts1[[grp]] <- fts1[[grp]][fts1[[grp]] %in% sigs1]
        fts2[[grp]] <- fts2[[grp]][fts2[[grp]] %in% sigs2]
      }
    }

    if(r_cutoff!=0) {
      # Get the correlations above an r cutoff
      rvals <- cormat$r
      rvals$Feature1 <- rownames(rvals)
      rvals <- rvals %>% pivot_longer(cols=colnames(cormat$r),names_to='Feature2',values_to='r')
      highrs1 <- rvals[abs(rvals$r)>=r_cutoff,]$Feature1
      highrs2 <- rvals[abs(rvals$r)>=r_cutoff,]$Feature2

      if(!length(c(highrs1,highrs2))) {
        message('\nNo correlations found with absolute R value above ', r_cutoff, appendLF = F)
        if(sigsOnly) message(' at a significance level of ', alpha, appendLF = F)
        message(' for ', grp)
      } else {
        fts1[[grp]] <- fts1[[grp]][fts1[[grp]] %in% highrs1]
        fts2[[grp]] <- fts2[[grp]][fts2[[grp]] %in% highrs2]
      }
    }
  }

  for(grp in names(cormats)) {
    for(m in names(cormats[[grp]])) {
      if(matchFts) {
        fts1[[grp]] <- unique(unlist(fts1))
        fts2[[grp]] <- unique(unlist(fts2))
      }
      cormats[[grp]][[m]] <- cormats[[grp]][[m]][rownames(cormats[[grp]][[m]]) %in% fts1[[grp]],
                                                 colnames(cormats[[grp]][[m]]) %in% fts2[[grp]]]
    }
  }

  return(cormats)
}
