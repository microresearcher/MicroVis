#' Generate a table of feature ratios
#'
#' @param melted_data Melted metadata and abundance table
#' @param feature1list List of numerator features
#' @param feature2list List of denominator features
#' @param factor (Optional) Factor along which to group samples to calculate an
#'     aggregated ratio
#' @param aggregateBefore Aggregate samples before calculating ratios? Stats will
#'     not be calculated between groups in this case.
#'
#' @return Table of feature ratios between specified features
#' @export
#'
ftRatiostab <- function(melted_data, feature1list, feature2list,
                        factor=NA, aggregateBefore=F) {
  if(length(feature1list)!=length(feature2list)) {
    message('\n"feature1list" and "feature2list" need to be the same length')
    return(melted_data)
  }

  if(!is.na(factor)) {
    if(!(factor %in% colnames(melted_data))) {
      message(factor,' was not found in the column names of the data')
      factor <- NA
    }
  }

  if(!is.na(factor) & aggregateBefore) {
    form <- formula(paste('.~',factor))
    melted_data <- aggregate(formula=form,
                             data=melted_data[c(factor,union(feature1list,feature2list))],
                             FUN=function(x) mean(x))
    allsamples <- NULL
    se <- NULL
  }
  for(i in 1:length(feature1list)) {
    ft1zeros <- sum(melted_data[[feature1list[i]]]==0)
    ft2zeros <- sum(melted_data[[feature2list[i]]]==0)
    if(ft2zeros>ft1zeros) {
      melted_data <- ftRatio(melted_data,feature2list[i],feature1list[i])
    }
    else melted_data <- ftRatio(melted_data,feature1list[i],feature2list[i])
  }
  ratiocols <- colnames(melted_data[grepl('/',colnames(melted_data))])
  melted_data[ratiocols][abs(melted_data[ratiocols])==Inf] <- NA
  if(!is.na(factor) & !aggregateBefore) {
    form <- formula(paste('.~',factor))
    allsamples <- melted_data
    se <- aggregate(formula=form,
                    data=melted_data[c(factor,ratiocols)],
                    FUN=function(x) sd(x,na.rm=T)/sqrt(length(x[!is.na(x)])))
    melted_data <- aggregate(formula=form,
                             data=melted_data[c(factor,ratiocols)],
                             FUN=function(x) mean(x), na.action = na.exclude)
  }


  if(is.na(factor)) {
    ratios <- list(ratios=ratiocols,data=melted_data[c('sample',ratiocols)])
  } else {
    means <- melted_data %>% pivot_longer(cols=ratiocols,
                                          names_to='ratio',
                                          values_to='Mean')
    SEs <- se %>% pivot_longer(cols=ratiocols,
                               names_to='ratio',
                               values_to='SE')
    mean.se <- merge(means,SEs)

    ratios <- list(ratios=ratiocols,mean.se=mean.se)
    if(!aggregateBefore) {
      ratios$data <- allsamples[c('sample',factor,ratiocols)]

      rationames <- gsub('/','__',ratios$ratios)

      ratiodata <- ratios$data
      colnames(ratiodata) <- gsub('/','__',colnames(ratiodata))

      stats <- univar(data=ratiodata,
                      factor=factor,
                      features=rationames)

      stats$stattab$Feature <- gsub('__','/',stats$stattab$Feature)
      names(stats$stats) <- gsub('__','/',names(stats$stats))

      allstats <- stats$stats[[1]][0,]
      for(ratio in names(stats$stats)) {
        stats$stats[[ratio]]$.y. <- gsub('__','/',stats$stats[[ratio]]$.y.)
        stats$stats[[ratio]] <- stats$stats[[ratio]] %>% add_xy_position(fun = 'mean_se')
        allstats <- full_join(allstats,stats$stats[[ratio]])
      }
      stats$stats <- allstats
      ratios$stats <- stats
    }
  }
  return(ratios)
}

#' Calculate ratio between two features
#'
#' @param data Abundance table
#' @param feature1 Numerator feature
#' @param feature2 Denominator feature
#'
#' @return Ratio of the two specified features
#' @export
#'
ftRatio <- function(data, feature1, feature2) {
  if(!(feature1 %in% colnames(data))) {
    message('\n',feature1,' could not be found in the column names of the melted_data')
    ft1valid <- F
  } else ft1valid <- T
  if(!(feature2 %in% colnames(data))) {
    message('\n',feature2,' could not be found in the column names of the melted_data')
    ft2valid <- F
  } else ft2valid <- T

  if(!ft1valid | !ft2valid) return(data)
  else data[[paste0(feature1,'/',feature2)]] <- data[[feature1]]/data[[feature2]]

  return(data)
}
