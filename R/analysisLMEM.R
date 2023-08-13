#' Linear Mixed Effects Model Discriminant Analysis
#'
#' @param dataset MicroVis dataset. Defaults to active dataset.
#' @param formula Formula for linear model. Defaults to simple linear model of all covariates.
#' @param exclude Factors/covariates to exclude in linear model.
#' @param rank Rank at which to conduct analysis.
#' @param alpha Significance threshold. Defaults to 0.05
#' @param is.winsor Whether to replace outliers (using winsorization)
#' @param zero.handling (From linda function) A character string of 'pseudo-count' or 'imputation' indicating the zero handling method used when feature.dat is 'count'. If 'pseudo-count', apseudo.cnt will be added to each value in feature.dat. If 'imputation', then we use the imputation approach using the formula in the referenced paper. Basically, zeros are imputed with values proportional to the sequencing depth. When feature.dat is 'proportion', this parameter will be ignored and zeros will be imputed by half of the minimum for each feature.
#'
#' @return Results of linear modeling with log2fc, p-values, adjusted p-values for each covariate.
#' @export
#'
mvLMEM <- function(dataset=NULL,
                   formula=NULL,
                   exclude=NULL,
                   rank=NULL,
                   zero.handling='pseudo-count',
                   alpha=0.05,
                   is.winsor=T) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(!is.null(formula)) formula <- checkFormula(dataset=dataset, formula=formula)
  if(is.null(formula)) {
    message('Creating basic formula from factors in dataset:')
    vars <- names(dataset$factors)[!(names(dataset$factors) %in% exclude)]
    vars <- vars[sapply(vars, function(x) length(dataset$factors[[x]]$subset) >= 2)]
    formula <- paste0('~',paste(vars, collapse = '+'))
  } else {
    vars <- getFormulaVars(formula)
    cat('\nFormula provided for analysis:\n')
  }

  cat(' ',formula,'\n\n')

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  dataset.temp <- clearNormalization(dataset, temp = T, silent = T)

  mdcolnum <- ncol(dataset.temp$metadata)

  data <- mvmelt(dataset.temp, rank = rank)
  md <- data[1:mdcolnum]
  abun <- data[(mdcolnum+1):ncol(data)]
  abun$Other <- NULL
  abun$Unknown <- NULL
  abun <- t(abun)

  linda.res <- MicrobiomeStat::linda(feature.dat = abun,
                                     meta.dat = md,
                                     formula = formula,
                                     zero.handling = zero.handling,
                                     is.winsor = is.winsor,
                                     alpha = 0.1)

  res <- lapply(linda.res$output, function(x) cbind(data.frame('Feature'=rownames(x)),x))

  res <- lapply(vars, function(x) {
    if(!any(grepl(x, names(res)))) stop('Some variables were not included in the linear model, possibly because one or more of the variables have too many groups.\n')

    if(x %in% names(dataset$factors)) grps <- dataset$factors[[x]]$groups
    else grps <- levels(factor(md[[x]]))

    ref_grp <- grps[[1]]
    comp_grps <- grps[!(grps %in% ref_grp)]

    # Determine which dataframes in the linDA results correspond to this variable
    # linDA will automatically take the first group in each factor as the reference group
    i <- match(paste0(x,grps[[2]]),names(res))
    j <- i+length(grps)-2

    temp <- lapply(names(res[i:j]), function(y) cbind(data.frame('Reference'=rep(ref_grp), times=nrow(res[[y]]),
                                                                 'Contrast'=rep(sub(x,'',y), times=nrow(res[[y]]))),
                                                      res[[y]]))
    names(temp) <- names(res[i:j])

    # This won't change anything if there is only one comparison group
    temp <- dplyr::bind_rows(temp)

    rownames(temp) <- NULL

    temp
  })
  names(res) <- vars

  if(is.null(dataset$stats$LMEM)) dataset$stats$LMEM <- list()
  if(is.null(dataset$stats$LMEM[[rank]])) dataset$stats$LMEM[[rank]] <- list()
  dataset$stats$LMEM[[rank]][[formula]] <- res

  assign('active_dataset',dataset,envir = mvEnv)
  if(!is.null(dataset$name)) assign(dataset$name,dataset,1)

  return(dataset)
}

#' Check that a formula is valid for a dataset
#'
#' @param dataset MicroVis dataset.
#' @param formula Formula for linear model.
#'
#' @return Formula if it is valid, or NULL if formula is not valid.
#'
checkFormula <- function(dataset, formula) {
  if(!substr(formula,1,1)=='~') {
    warning(' Formula must start with a tilde (~).\n')
    return(NULL)
  }

  md <- dataset$metadata
  # factors <- dataset$factors
  allvariables <- names(md)[2:ncol(md)]

  vars <- getFormulaVars(formula)

  if(!all(vars %in% allvariables)) {
    warning(' The provided formula contains variables that are not part of this dataset.\n')
    return(NULL)
  }

  # valid_factors <- sapply(vars, function(x) length(factors[[x]]$subset) > 1)
  valid_vars <- sapply(vars, function(x) length(unique(md[[x]])) > 1)
  if(!all(valid_vars)) {
    warning(' Fewer than 2 groups remaining in "', names(valid_vars)[!valid_vars], '".\n', immediate. = T)
    return(NULL)
  }

  return(formula)
}

#' Extract Variables from a Formula
#'
#' @param formula Formula as a string.
#'
#' @return Variables extracted from the formula in a vector
#'
getFormulaVars <- function(formula) {
  operation_chars <- c('~' = '',
                       '\\+' = ',',
                       '\\-' = ',',
                       '\\*' = ',',
                       '\\/' = ',',
                       '\\^' = ',',
                       '\\(' = ',',
                       '\\)' = ',',
                       '\\|' = ',')

  vars <- unlist(strsplit(stringr::str_replace_all(formula,operation_chars),','))
  vars <- vars[!(vars=='')]
  vars <- vars[is.na(suppressWarnings(as.numeric(vars)))]

  return(vars)
}
