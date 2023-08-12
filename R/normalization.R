#' @title Data Normalization
#'
#' @description Normalizes unranked data based on parameters set by the
#'     scaleSamples, scaleFeatures, and transData. This function is one of the
#'     processing functions called by processDataset().
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp This parameter has no use in this function and can be removed
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with normalized data
#'
runNormalization <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  if(!is.null(dataset$data$proc$rarefied)) return(dataset)

  if(is.null(dataset$data$proc$unranked)) dataset <- runSampleFilter(dataset, silent = T)

  ft_data <- dataset$data
  ranks <- names(ft_data$proc)[names(ft_data$proc) %in% c(get('taxaRanks',envir = mvEnv),
                                                          'single_rank', 'functional', 'unranked')]
  normalization <- ft_data$proc$normalization

  ### Check Normalization Parameters ###
  #------------------------------------#
  # Make sure that transformation is last in the normalization list
  #   by making it NULL and then putting it back
  if(!is.null(normalization$transformation)) {
    # If there are transformation parameters, make sure transformation is
    #   the LAST element in the normalization list
    if(grep('transformation',names(normalization))!=length(normalization)) {
      transformation <- normalization$transformation
      normalization$transformation <- NULL
      normalization$transformation <- transformation
    }
  }

  if(!is.null(normalization) & !silent) {
    cat(paste0('\n\n|~~~~~~~~~~~~~~  NORMALIZING DATA  ~~~~~~~~~~~~~~|\n'))
  }

  ### Normalize ###
  #---------------#
  # This loop goes through the recorded normalizations IN ORDER
  #   so that it performs feature scaling before sample scaling or vice versa
  #   depending on which one was performed first
  for(norm in names(normalization)) {
    if(norm=='sample_scale') {
      if(!silent) cat('\n  Scaling samples by',normalization[[norm]]$scaling)
      rn <- rownames(ft_data$proc$unranked)
      cn <- colnames(ft_data$proc$unranked)
      # Transpose the dataframe to perform the column-wise functions on samples
      normalized <- data.frame(t(normalizer(data.frame(t(ft_data$proc$unranked)),
                                            norm_type = normalization[[norm]]$scaling,
                                            sum_scale = normalization[[norm]]$sum_scale)))
      rownames(normalized) <- rn
      colnames(normalized) <- cn
      ft_data$proc$unranked <- normalized
    } else if(norm=='feature_scale') {
      if(!silent) cat('\n  Scaling features by',normalization[[norm]]$scaling)
      rn <- rownames(ft_data$proc$unranked)
      cn <- colnames(ft_data$proc$unranked)
      # Perform the column-wise functions on features
      normalized <- normalizer(ft_data$proc$unranked,
                               norm_type = normalization[[norm]]$scaling,
                               sum_scale = normalization[[norm]]$sum_scale)
      rownames(normalized) <- rn
      colnames(normalized) <- cn
      ft_data$proc$unranked <- normalized
    } else if(norm=='transformation') {
      if(!silent) cat('\n  Transforming data by',normalization[[norm]]$trans_method)
      # If performing CLR transformation:
      #   1) Remove any filtering options that no longer make sense with CLR
      #   2) Remove any features with all zeros
      if(normalization[[norm]]$trans_method=='clr') {
        if(!silent) cat('\n    Note: We do not currently recommended using filtering with CLR normalization')
        # if(!silent) cat('\n    Removing feature filters except NA, top prevalence, and variance filters, if they exist')
        # ft_data$proc$filtering[!(names(ft_data$proc$filtering) %in% c('NAfilter',
        #                                                               'top_prevalence',
        #                                                               'low_var_percentile',
        #                                                               'top_var'))] <- NULL

        zerofts <- colnames(ft_data$proc$unranked)[colSums(ft_data$proc$unranked)==0]
        if(!silent) cat('\n    Removing',length(zerofts),'features absent in all subsetted samples')
        ft_data$proc$unranked <- ft_data$proc$unranked[!(colnames(ft_data$proc$unranked) %in% zerofts)]
        ft_data$proc$filtering$filterlist$zeros <- ASVtoTaxa(ft_data,
                                                             zerofts,
                                                             getLowestRank(dataset))
      }

      rn <- rownames(ft_data$proc$unranked)
      cn <- colnames(ft_data$proc$unranked)

      # Perform a function to all values
      # This one is transposed too because GMPR and RLE transformation
      #   are performed with samples as columns
      normalized <- data.frame(t(normalizer(data.frame(t(ft_data$proc$unranked)),
                                            norm_type = normalization[[norm]]$trans_method,
                                            log_base = normalization[[norm]]$log_base)))
      rownames(normalized) <- rn
      colnames(normalized) <- cn
      ft_data$proc$unranked <- normalized

      # Ensure that, if CLR transformation was performed, features with all zeros
      #   are forcibly filtered out
      ft_data$proc
    }
  }

  if(length(normalization)) {
    # Reset normalization in case transformation needed to be moved to the end
    ft_data$proc$normalization <- normalization
    dataset$data <- ft_data
    if(!silent) cat('\n\n>>> DATA NORMALIZATION SUCCESSFUL <<<\n')
  } else if(!silent) cat('\n~~~ No normalization performed ~~~\n')

  return(dataset)
}

#' Clear Normalization
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp If set to TRUE, it tells processDataset() to NOT update the active
#'     dataset.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with all normalization parameters
#'     cleared. This function then calls processDataset() which will call the
#'     3 "run___" functions and therefore reset all abundance values to raw values.
#'
#' @export
#'
clearNormalization <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  dataset$data$proc$normalization <- NULL

  dataset <- processDataset(dataset, temp=temp, silent=silent)

  return(dataset)
}


#' Scale Samples
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp If set to TRUE, it tells processDataset() to NOT update the active
#'     dataset.
#' @param scaling Method for scaling samples. Choices are "sum", "none", "relative",
#'     "median", "range", "centered", and "standardized". Default is "sum".
#' @param sum_scale Number of reads to scale each sample to. Default is 1 million.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return Microvis dataset (mvdata object) with samples scaled.
#' @export
#'
scaleSamples <- function(dataset=NULL, temp=F,
                         scaling='sum',
                         sum_scale=1e6,
                         silent=F) {

  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  dataset$data$proc$rarefied <- NULL

  ### Process Input Arguments ###
  #-----------------------------#
  # First get the recorded sample scaling parameters (NULL if they don't exist)
  sample_scale <- dataset$data$proc$normalization$sample_scale
  scaletypes <- c('sum', 'median', 'range', 'centered', 'standardized', 'relative', 'none')
  # If no valid input parameters were provided:
  #   1) default to recorded sample scaling parameters, if they exist
  #     or
  #   2) have user select a sample scaling option
  if(!(scaling[[1]] %in% scaletypes)) {
    if(!is.null(sample_scale)) {
      scaling <- sample_scale$scaling
      sum_scale <- sample_scale$sum_scale
    } else {
      message('\nPlease choose a method to scale samples by\n')
      Sys.sleep(0.1)
      scaling <- select.list(scaletypes, title = 'Sample scaling method:', graphics = T)
    }
  }

  ### Format Sample Scaling Parameters ###
  #--------------------------------------#
  if(scaling=='none' | scaling=='') {
    sample_scale <- NULL
  } else {
    if(scaling!='sum') {
      sum_scale <- NULL
    }
    sample_scale <- list(scaling=scaling, sum_scale=sum_scale)
  }

  ### Record Sample Scaling Parameters in Dataset ###
  #-------------------------------------------------#
  # Store the new sample scaling parameters in the normalization parent variable
  dataset$data$proc$normalization$sample_scale <- sample_scale

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, temp = temp, silent = silent)

  return(dataset)
}

#' Scale Features
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp If set to TRUE, it tells processDataset() to NOT update the active
#'     dataset.
#' @param scaling Method for scaling features. Choices are "sum", "none", "relative",
#'     "median", "range", "centered", and "standardized". Default is "sum".
#' @param sum_scale Number of reads to scale each feature to. Default is 1 million.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return Microvis dataset (mvdata object) with features scaled.
#' @export
#'
scaleFeatures <- function(dataset=NULL, temp=F,
                          scaling=c('sum',
                                    'median',
                                    'range',
                                    'centered',
                                    'standardized',
                                    'relative',
                                    'none'),
                          sum_scale=1e6,
                          silent=F) {

  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  dataset$data$proc$rarefied <- NULL

  ### Process Input Arguments ###
  #-----------------------------#
  # First get the recorded feature scaling parameters (NULL if they don't exist)
  feature_scale <- dataset$data$proc$normalization$feature_scale
  scaletypes <- c('sum', 'median', 'range', 'centered', 'standardized', 'relative', 'none')
  # If no valid input parameters were provided:
  #   1) default to recorded feature scaling parameters, if they exist
  #     or
  #   2) have user select a feature scaling option
  if(length(scaling)>1 | !(scaling[[1]] %in% scaletypes)) {
    if(!is.null(feature_scale)) {
      scaling <- feature_scale$scaling
      sum_scale <- feature_scale$sum_scale
    } else {
      message('\nPlease choose a method to scale features by\n')
      Sys.sleep(0.1)
      scaling <- select.list(scaletypes, title = 'Feature scaling method:', graphics = T)
    }
  }

  ### Format Feature Scaling Parameters ###
  #---------------------------------------#
  if(scaling=='none' | scaling=='') {
    feature_scale <- NULL
  } else {
    if(scaling!='sum') {
      feature_scale <- NULL
    }
    feature_scale <- list(scaling=scaling, sum_scale=sum_scale)
  }

  ### Record Feature Scaling Parameters in Dataset ###
  #--------------------------------------------------#
  # Store the new feature scaling parameters in the normalization parent variable
  dataset$data$proc$normalization$feature_scale <- feature_scale

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, temp, silent=silent)

  return(dataset)
}

#' Data Transformation
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp If set to TRUE, it tells processDataset() to NOT update the active
#'     dataset.
#' @param trans_method Method for transforming data. Choices are "glog"
#'     (generalized log), "none", "pseudolog", and "log".
#' @param log_base Base of log transformation.
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#' @param impute_method Method to impute zeros by if using zCompositions approach
#'
#' @return Microvis dataset (mvdata object) with samples scaled.
#' @export
#'
transData <- function(dataset=NULL, temp=F,
                      trans_method=c('clr', 'glog', 'pseudolog', 'log', 'gmpr', 'none'),
                      log_base=10,
                      impute_method=c('GBM','SQ','BL','CZM'),
                      silent=F) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  dataset$data$proc$rarefied <- NULL

  ### Process Input Arguments ###
  #-----------------------------#
  # First get the recorded data transformation parameters (NULL if they don't exist)
  transformation <- dataset$data$proc$normalization$transformation
  transtypes <- c('clr', 'glog', 'pseudolog', 'log', 'gmpr', 'none')
  # If no valid input parameters were provided:
  #   1) default to recorded data transformation parameters, if they exist
  #     or
  #   2) have user select a data transformation option
  if(length(trans_method)>1 | !(trans_method[[1]] %in% transtypes)) {
    if(!is.null(transformation)) {
      trans_method <- transformation$trans_method
      log_base <- transformation$log_base
    } else {
      message('\nPlease choose a method to transform data by\n')
      Sys.sleep(0.1)
      trans_method <- select.list(transtypes, title = 'Data transformation method:', graphics = T)
    }
  }

  ### Format Data Transformation Parameters ###
  #-------------------------------------------#
  if(trans_method=='none' | trans_method=='') {
    transformation <- NULL
  } else {
    if(!(trans_method %in% c('clr', 'glog', 'pseudolog', 'log'))) {
      log_base <- NULL
    }
    if(trans_method %in% c('clr','gmpr','rle')) {
      if(!is.null(dataset$data$proc$normalization$sample_scale)) {
        dataset$data$proc$normalization$sample_scale <- NULL
        message('\nWARNING: ',trans_method,' overrides data scaling - no sample scaling will be done')
      }
      if(!is.null(dataset$data$proc$normalization$feature_scale)) {
        dataset$data$proc$normalization$feature_scale <- NULL
        message('\nWARNING: ',trans_method,' overrides data scaling - no feature scaling will be done')
      }
    }
    transformation <- list(trans_method=trans_method, log_base=log_base)
  }

  ### Record Data Transformation Parameters in Dataset ###
  #------------------------------------------------------#
  # Store the new sample_scaling parameters in the normalization parent variable
  dataset$data$proc$normalization$transformation <- transformation

  if(!get('.loading',envir = mvEnv)) dataset <- processDataset(dataset, temp = temp, silent = silent)

  return(dataset)
}

#' @title Normalizing Function
#'
#' @description Actual function that normalizes data, one method at a time. In
#'     other words, it is called three times by runNormalization() if the data
#'     is to be scaled by samples, scaled by features, and transformed.
#'
#' @param data Abundance table that will be normalized by in columns. In other
#'     words, if scaling by samples, then samples should be in columns. Does not
#'     matter for data transformation.
#' @param norm_type Normalization type: "sum", "none", "relative", "standardized",
#'     "centered", "median", "range", "glog" (generalized log), "pseudolog", "log",
#' @param sum_scale Number to scale samples/features to.
#' @param log_base Base for log transformations.
#' @param impute_method Method to impute zeros by if using zCompositions approach
#'
#' @return Normalized abundance table.
#'
normalizer <- function(data,
                       norm_type,
                       sum_scale=NA,
                       log_base=NA,
                       impute_method=c('GBM','SQ','BL','CZM')) {
  norm_type <- tolower(norm_type)

  rn <- rownames(data)
  cn <- colnames(data)

  if(is.null(norm_type)) {
    return(data)

  } else if(norm_type=='relative') {
    # x/sum(x)
    normalized <- data.frame(sapply(data, function(x) x/sum(x,na.rm = T)) )
    normalized[is.na(normalized)] <- 0

  } else if (norm_type=='sum') {
    # 1000*x/sum(x)
    normalized <- data.frame(sapply(data, function(x) sum_scale*x/sum(x,na.rm = T)) )
    normalized[is.na(normalized)] <- 0

  } else if(norm_type=='standardized') {
    # (x-mean)/sd
    normalized <- data.frame(sapply(data, function(x) (x-mean(x,na.rm = T))/sd(x,na.rm = T) ))
    normalized[is.na(normalized)] <- 0

  } else if(norm_type=='centered') {
    # x-mean
    normalized <- data.frame(sapply(data, function(x) x-mean(x,na.rm = T) ))

  } else if(norm_type=='median') {
    # x-median
    normalized <- data.frame(sapply(data, function(x) x-median(x,na.rm = T) ))

  } else if(norm_type=='range') {
    # (x-mean)/range
    normalized <- data.frame(sapply(data, function(x) (x-mean(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)) ))
    normalized[is.na(normalized)] <- 0

  } else if(norm_type=='clr') {
    # Centered-log ratio
    #   First, zeros are replaced
    #   Second, centered-log ratio is performed
    #   Finally, those features with all zeros are re-attached to the data.frame
    #     These features will be forcibly removed during the filtering step

    data.nozeros <- zeroReplace(data)

    normalized <- data.frame(t(apply(data.nozeros, 1,
                                     function(x) log(x/exp(mean(log(x))), base=log_base) )))

    normalized <- cbind(normalized)

  } else if(norm_type=='glog') {
    # Generalized log formula whose domain includes ALL real numbers,
    #   positive numbers map to positive numbers,
    #   negative numbers map to negative numbers,
    #   and zero maps to zero! :)
    normalized <- data.frame(sapply(data, function(x) log((x+sqrt(x^2+4))/2, base=log_base) ))

  } else if(norm_type=='pseudolog') {
    # Take the log10 of all values
    normalized <- data.frame(sapply(data, function(x) log(x+1, base=log_base) ))

  } else if(norm_type=='log') {
    # Compute half of the minimum non-zero abundance value and
    # add this 'halfmin' to all abundance values
    halfmin <- min(data[data>0],na.rm = T)/2
    data.temp <- data.frame(sapply(data, function(x) x+halfmin))
    # Take the log10 of all values
    normalized <- data.frame(sapply(data.temp, function(x) log(x, base=log_base) ))

  } else if(norm_type=='gmpr') {
    # ### GMPR normalization ###
    # gmpr.intersect <- floor(nrow(data)/10)
    # if(gmpr.intersect<1) gmpr.intersect <- 1
    # gmpr <- calcGMPR(counts=data, intersect.no=gmpr.intersect)
    # normalized <- data.frame(data/gmpr)

    # } else if(norm_type=='rle') {
    #   ### RLE Normalization ###
    #   # First, create a DGE list
    #   dges <- DGEList(counts=data, remove.zeros = TRUE)
    #   # Then perform RLE normalization
    #   normfactors <- calcNormFactors(dges, method='RLE')$samples$norm.factors
    #   normalized <- data.frame(data*normfactors)

  } else if(norm_type=='none') {
    return(data)

  } else {
    message('\nWARNING: Normalization type not recognized. No normalization will be performed at this step\n')
    warning_list <- c(warning_list, paste0('WARNING: Normalization type not recognized. No normalization will be performed at this step'))
    return(data.frame(data))

  }

  rownames(normalized) <- rn
  colnames(normalized) <- cn

  return(normalized)
}

#' Normalize a Standalone Abundance Table
#'
#' @param abundance_table Abundance table with samples as rows and features as
#'     columns
#' @param sample_scale Sample scaling method to be used. Choices are "sum", "none",
#'     "relative", "median", "range", "centered", and "standardized". Default is
#'     "sum".
#' @param feature_scale Feature scaling method to be used. Choices are "sum", "none",
#'     "relative", "median", "range", "centered", and "standardized". Default is
#'     "sum".
#' @param transformation Transformation method to be used. Choices are "glog"
#'     (generalized log), "none", "pseudolog", and "log".
#' @param sum_scale Number to scale samples/features by. Default is 1 million.
#' @param log_base Base for log transformations. Default is 10
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return Normalized abundance table.
#' @export
#'
normalizeTable <- function(abundance_table,
                           sample_scale=NULL,
                           feature_scale=NULL,
                           transformation=NULL,
                           sum_scale=1e6, log_base=10,
                           silent=F) {
  rn <- rownames(abundance_table)
  cn <- colnames(abundance_table)

  if(!is.null(sample_scale)) {
    abundance_table <- data.frame(t(normalizer(data.frame(t(abundance_table)),
                                               norm_type = sample_scale, sum_scale = sum_scale)))
  }
  if(!is.null(feature_scale)) {
    abundance_table <- normalizer(data.frame(abundance_table),
                                  norm_type = feature_scale, sum_scale = sum_scale)
  }
  if(!is.null(transformation)) {
    abundance_table <- data.frame(t(normalizer(data.frame(t(abundance_table)),
                                               norm_type = transformation, log_base = log_base)))
  }

  rownames(abundance_table) <- rn
  colnames(abundance_table) <- cn

  return(abundance_table)
}

#' Generalized Log
#'
#' @param x Any real number
#' @param base Log base. Defaults to 10
#'
#' @return Generalized log value of the input
#' @export
#'
glog <- function(x, base=10) {
  return(log((x+sqrt(x^2+4))/2, base=base))
}

#' Replace Zeros with Non-Zero Value
#'
#' @param x Numerical matrix
#'
#' @return Matrix with no zeroes
#'
zeroReplace <- function(x) {
  if(get('zeroReplaceMethod', envir=mvEnv)$method=='replace') {
    div <- get('zeroReplaceMethod', envir=mvEnv)$div
    rdist <- get('zeroReplaceMethod', envir = mvEnv)$rdist

    r <- nrow(x)
    c <- ncol(x)

    min.nonzero <- min(x[x>0])/div

    if(rdist=='runif') {
      x.nozeros <- x + ((x==0)*matrix(runif(r*c, max=min.nonzero),
                                      nrow = r))
    } else if(rdist=='point') {
      x.nozeros <- x + ((x==0)*matrix(rep(min.nonzero,r*c),
                                      nrow = r))
    }

  } else if(get('zeroReplaceMethod', envir=mvEnv)$method=='impute') {
    x.nozeros <- data.frame(zCompositions::cmultRepl(x,
                                                     output = 'p-counts',
                                                     suppress.print = T))
  }


  return(x.nozeros)
}
