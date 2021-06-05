#' Zero-Inflated Mixed Model Analysis
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor to analyze by. Defaults to the active factor
#' @param rank Rank to select features from. Defaults to active rank
#' @param dataset_name (Not Recommended) Name of the dataset to save results
#'     to. This should not need to be used by users since the function can
#'     determine the name of the dataset directly passed to it, but not when
#'     it is called within another function.
#'
#' @return MicroVis dataset containing results of LEfSe analysis
#' @export
#'
mvZimm <- function(dataset=NULL,
                   factor=NULL,
                   rank=NULL,
                   dataset_name=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else if(is.null(dataset_name)) {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  ngrps <- length(dataset$factors[[factor]]$subset)
  if(ngrps>2) {
    if(is.null(dataset$stats[[factor]][[rank]]$zig)) {
      cat('\nPerforming zero-inflated gaussian fitting\n')
      res <- mvZig(dataset,factor,rank)
      dataset$stats[[factor]][[rank]]$zig <- res
    }
  } else {
    if(is.null(dataset$stats[[factor]][[rank]]$ziln)) {
      cat('\nPerforming zero-inflated log-normal fitting\n')
      res <- mvZiln(dataset,factor,rank)
      dataset$stats[[factor]][[rank]]$ziln <- res
    }
  }

  if(dataset_name=='active_dataset') assign(dataset_name, dataset,envir = mvEnv)
  else assign(dataset_name,dataset,1)

  return(dataset)
}

#' Zero-inflated Log-Normal Fitting and Statistical Analysis of Differences
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor to analyze by. Defaults to the active factor
#' @param rank Rank to select features from. Defaults to active rank
#'
#' @return metagenomeSeq fitZig results
#' @export
#'
mvZiln <- function(dataset=NULL,factor=NULL,rank=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else if(is.null(dataset_name)) {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  ngrps <- length(dataset$factors[[factor]]$subset)
  if(ngrps>2) stop('Zero-inflated log-normal fitting currently only works for 2 groups')

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  data <- mvmelt(dataset, rank=rank)
  sample_names <- data$sample

  metadata <- data[2:ncol(dataset$metadata)]
  abun <- data[(ncol(dataset$metadata)+1):ncol(data)]
  abun$Other <- NULL
  abun <- data.frame(t(abun))

  colnames(abun) <- rownames(metadata) <- sample_names

  abun <- as.matrix(abun)

  formula <- as.formula(paste('~1+',factor))
  mod <- model.matrix(formula,data=metadata)
  fit <- mvfitFeatureModel(abun,mod)

  return(mvMRcoefs(fit,number=Inf))
}

#' Zero-inflated Gaussian Fitting and Statistical Analysis of Differences
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor to analyze by. Defaults to the active factor
#' @param rank Rank to select features from. Defaults to active rank
#'
#' @return metagenomeSeq fitFeatureModel results
#' @export
#'
mvZig <- function(dataset=NULL,factor=NULL,rank=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else if(is.null(dataset_name)) {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- factor[factor %in% names(dataset$factors)]
  if(is.null(factor)) factor <- dataset$active_factor

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  data <- mvmelt(dataset, rank=rank)
  sample_names <- data$sample

  metadata <- data[2:ncol(dataset$metadata)]
  abun <- data[(ncol(dataset$metadata)+1):ncol(data)]
  abun$Other <- NULL
  abun <- data.frame(t(abun))

  colnames(abun) <- rownames(metadata) <- sample_names

  abun <- as.matrix(abun)

  factor <- dataset$active_factor

  mod <- model.matrix(~metadata[[factor]])
  fit <- mvfitZig(abun,mod)

  return(mvMRcoefs(fit,number=Inf))
}
