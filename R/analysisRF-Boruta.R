#' Perform random forest Boruta analysis
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param factor Factor to analyze by. Defaults to the active factor
#' @param max_runs Maximum number of runs for Boruta algorithm. Default is 100
#' @param roughfix Perform roughfix of indeterminate features? Defaults to FALSE
#' @param alpha Significance threshold. Defaults to 0.01
#'
#' @return List of 3 tables with the class error, errors, and the Boruta results
#' @export
#'
rfboruta <- function(dataset=NULL,
                     factor=NULL,
                     max_runs=100,roughfix=F,
                     alpha=0.01) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  factor <- dataset$active_factor
  mdcolnum <- ncol(dataset$metadata)

  data <- mvmelt(dataset)
  abun <- data[(mdcolnum+1):ncol(data)]
  abun$Other <- NULL
  abun$Unknown <- NULL
  groups <- data[[factor]]

  invalid_chars <- c(' '='__1__','-'='__2__','\\+'='__3__')
  groups <- factor(str_replace_all(groups,invalid_chars))

  rferrs <- randomForest(abun, groups)
  errs <- data.frame(rferrs$err.rate)
  class_err <- data.frame(rferrs$confusion)

  res <- Boruta(abun, groups, maxRuns=max_runs, pValue=alpha)
  if(roughfix) res <- TentativeRoughFix(res)

  decisions <- data.frame(Feature=names(res$finalDecision),
                          Decision=res$finalDecision)
  decisions <- rbind(decisions,data.frame(Feature=c('shadowMin','shadowMean','shadowMax'),
                                          Decision=rep('Shadow',3),
                                          row.names = c('shadowMin','shadowMean','shadowMax')))
  imptab <- apply(res$ImpHistory,2,function(x) x[is.finite(x)])
  imptab <- imptab[order(sapply(imptab,function(x) median(x)))]
  imptab <- stack(imptab)
  colnames(imptab) <- c('Importance','Feature')
  boruta <- full_join(imptab,decisions)

  replace_chars <- c('__1__'=' ','__2__'='-','__3__'='\\+')
  colnames(errs) <- str_replace_all(colnames(errs),replace_chars)
  rownames(class_err) <- str_replace_all(rownames(class_err),replace_chars)
  colnames(class_err) <- str_replace_all(colnames(class_err),replace_chars)
  for(i in 1:length(levels(groups))) {
    if(is.numeric(type.convert(substr(levels(groups)[i],1,1)))) {
      colnames(errs)[i+1] <- sub('X','',colnames(errs)[i+1])
      colnames(class_err)[i] <- sub('X','',colnames(class_err)[i])
    }
  }

  return(list(ClassError=class_err,Errors=errs,Boruta=boruta))
}
