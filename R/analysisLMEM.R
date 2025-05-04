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
  if(!requireNamespace('MicrobiomeStat', quietly = T)) {
    stop('You need to install "MicrobiomeStat" package in order to perform mixed linear model analysis.')
  }

  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(!is.null(formula)) formula <- checkFormula(dataset=dataset, formula=formula)
  if(is.null(formula)) {
    message('Creating basic formula from factors in dataset:')
    vars <- names(dataset$factors)[!(names(dataset$factors) %in% exclude)]
    vars <- vars[sapply(vars, function(x) length(dataset$factors[[x]]$subset) >= 2)]
    formula <- paste0('~',paste(vars, collapse = '+'))
  } else {
    vars <- getFormulaVars(formula)
    na_terms <- vars[!(vars %in% colnames(dataset$metadata))]
    if(length(na_terms)) stop(paste('The following terms are not variables (column names in metadata) in the dataset:',
                                    paste(na_terms, sep = ',')))
    terms <- parseFormula(formula)
    cat('\nFormula provided for analysis:\n')
  }

  randomEffects <- grepl('\\|', formula)

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
                                     alpha = alpha)

  # if(randomEffects) temp <- lapply(rownames(abun), function(ft) {
  #   fit <- suppressMessages(lme4::lmer(formula = as.formula(paste(ft,formula)),
  #                                      data = data))
  #   list(estimate = coef(summary(fit)), covar = stats::vcov(fit))
  # }) else fit <- lapply(rownames(abun), function(ft) {
  #   stats::lm(formula = as.formula(paste(ft,formula)),
  #             data = data)
  # })
  #
  # names(temp) <- rownames(abun)
  #
  # res <- do.call(rbind, lapply(temp, `[[`, 1))
  # res.cov <- do.call(rbind, lapply(temp, `[[`, 2))
  #
  # estimate <- lapply(unique(rownames(res)), function(v) {
  #   res.voi <- res[which(rownames(res) == x), ]
  #   rownames(res.voi) <- NULL
  #
  #   if(random.effect) {
  #     df <- res.voi[, 3]
  #   }
  #
  #   log2FoldChange <- res.voi[, 1]
  #   lfcSE <- res.voi[, 2]
  #
  #   bias <- suppressMessages(modeest::mlv(sqrt(nrow(md)) * log2FoldChange,
  #                                         method = 'meanshift',
  #                                         kernel = 'gaussian') / sqrt(nrow(md)))
  #
  #   log2FoldChange <- log2FoldChange - bias
  #   stat <- log2FoldChange / lfcSE
  #
  #   pvalue <- 2 * pt(-abs(stat), df)
  #   padj <- p.adjust(pvalue, method = p.adj.method)
  #   reject <- padj <= alpha
  #
  #   output <- cbind.data.frame(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, reject, df)
  #   rownames(output) <- taxa.name
  #
  #   return(list(bias = bias, output = output))
  # })

  res <- lapply(linda.res$output, function(x) cbind(data.frame('Feature'=rownames(x)),x))

  # Loop over only the variables that are also their own terms in the formula
  #  Do not include those denoted as intercepts with vertical bar
  res <- lapply(intersect(vars, terms), function(v) {
    if(!any(grepl(v, names(res)))) stop('Some variables were not included in the linear model, possibly because one or more of the variables have too many groups.\n')

    # If the variable was continuous then
    if(v %in% names(dataset$factors)) grps <- dataset$factors[[v]]$groups
    else if(is.character(md[[v]])) grps <- levels(factor(md[[v]]))
    else grps <- c(v, '')

    ref_grp <- grps[[1]]
    comp_grps <- setdiff(grps, ref_grp)

    # Determine which dataframes in the linDA results correspond to this variable
    # linDA will automatically take the first group in each factor as the reference group
    i <- match(paste0(v,grps[[2]]),names(res))
    j <- i + length(grps)-2

    temp <- lapply(names(res[i:j]), function(y) cbind(data.frame('Reference' = rep(ref_grp, times=nrow(res[[y]])),
                                                                 'Contrast' = rep(sub(v,'',y), times=nrow(res[[y]]))),
                                                      res[[y]]))
    names(temp) <- names(res[i:j])

    # This won't change anything if there is only one comparison group
    temp <- dplyr::bind_rows(temp)

    # Get the plain fold-change values
    temp$FC <- 2 ^ (temp$log2FoldChange)

    rownames(temp) <- NULL

    temp
  })
  names(res) <- intersect(vars, terms)

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
  valid_vars <- sapply(vars, function(var) length(unique(md[[var]])) > 1)
  if(!all(valid_vars)) {
    warning(' Fewer than 2 groups remaining in "', names(valid_vars)[!valid_vars], '".\n', immediate. = T)
    return(NULL)
  }

  return(formula)
}

#' Extract Terms from a Formula
#'
#' @param formula Formula as a string.
#'
#' @return Variables extracted from the formula in a vector
#'
parseFormula <- function(formula) {
  operation_chars <- c('~' = '',
                       '\\+' = ',',
                       '\\-' = ',',
                       '\\*' = ',',
                       '\\/' = ',',
                       '\\^' = ',',
                       '\\(' = ',',
                       '\\)' = ',')

  terms <- unlist(strsplit(stringr::str_replace_all(formula, operation_chars),','))
  terms <- terms[!(terms=='')]

  return(terms)
}

#' Extract Variables from a Formula
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

  vars <- unlist(strsplit(stringr::str_replace_all(formula, operation_chars),','))
  vars <- vars[!(vars=='')]
  vars <- vars[is.na(suppressWarnings(as.numeric(vars)))]

  return(vars)
}
