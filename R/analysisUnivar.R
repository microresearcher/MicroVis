#' Perform Univariate Analysis
#'
#' @param data Either a MicroVis dataset or a data table containing metadata and
#'     count values with samples as rows and metadata factors or features as columns
#' @param factor Factor along which to form groups from the samples and perform
#'     the statistical analysis
#' @param stratifiers  (Optional) One or two factors to stratify the groups by
#' @param rank Rank at which to perform the dataset. Defaults to the active rank
#' @param features List of features to analyze. Defaults to all
#' @param pairwise_comparisons For 3 or more groups, list of vectors of pairwise
#'     comparisons to conduct. Defaults to all possible pairwise comparisons
#' @param param Run parametrized or nonparametrized univariate analysis? Default
#'     is FALSE
#' @param dataset_name (Not recommended) Name of the dataset to save statistics
#'     to. This should not need to be used by users since the function can
#'     determine the name of the dataset directly passed to it, but not when
#'     it is called within another function.
#'
#' @return MicroVis dataset containing the statistics results in its "stats"
#'     attribute
#' @export
#'
# TODO: Figure out how to incorporate Tukey and Games-Howell post-hoc tests as
#       options
univar <- function(data=NULL,
                   factor=NULL,
                   stratifiers=NULL,
                   rank=NULL,
                   features=NULL, pairwise_comparisons=NULL,
                   param=F,
                   dataset_name=NULL) {
  if(is.null(data)) data <- get('active_dataset',envir = mvEnv)

  # Can analyze either a MicroVis dataset or just a melted table of metadata and counts
  if(inherits(data,'mvdata')) {
    dataset <- data

    # Set the factor
    factor <- factor[factor %in% names(dataset$factors)]
    if(!length(factor)) factor <- dataset$active_factor

    # Set any stratifiers if specified
    stratifiers <- stratifiers[stratifiers %in% names(dataset$factors) & !(stratifiers %in% factor)]
    if(length(stratifiers)>2) {
      stratifiers <- stratifiers[1:2]
      message('\n  Note: Only up to 2 stratifiers can be selected. ',
              paste0(stratifiers,collapse = ' and '),' have been automatically selected.')
    }

    # Set the taxonomic rank
    rank <- rank[rank %in% getRanks(dataset)]
    if(is.null(rank)) rank <- dataset$data$proc$active_rank

    # Set the features to be compared
    features <- features[features %in% getFeatures(dataset,ranks=rank)]
    if(is.null(features)) features <- getFeatures(dataset,ranks=rank)
    features <- features[!(features %in% c('Other','Unknown'))]

    # Get the sizes of each group/stratified group
    anysmallgrps <- checkGroups(dataset=dataset,
                                factor=factor,
                                stratifiers=stratifiers,
                                min_n=3,
                                verbose=F)

    if(anysmallgrps) stop('All groups or stratified groups must have at least 3 samples for statistical analysis')
    if(is.null(stratifiers)) melted <- mvmelt(dataset, rank=rank)
    else melted <- mvstratify(dataset, stratifiers=stratifiers, rank=rank)

  } else {
    # Make sure the factor exists in the data
    if(!(factor %in% colnames(data))) stop(paste0(factor,' is not a column name of the data table'))

    # Make sure stratifiers exist in the data, are not the same as the factor,
    #   and that a maximum of 2 are selected
    stratifiers <- stratifiers[stratifiers %in% colnames(data) & !(stratifiers %in% factor)]
    if(length(stratifiers)>2) {
      stratifiers <- stratifiers[1:2]
      message('\n  Note: Only up to 2 stratifiers can be selected. ',
              paste0(stratifiers,collapse = ' and '),' have been automatically selected.')
    }

    # Make sure there are at least 2 groups in the factor
    if(length(unique(data[[factor]]))<2) stop('There is only 1 group in the ',factor,' factor')

    # Make sure each group or stratum has at least 3 samples
    grp_sizes <- countSamples.base(metadata=data[c(factor,stratifiers)],
                                   factors=factor,
                                   stratifiers=stratifiers,
                                   verbose=F)
    if(is.null(stratifiers)) grp_sizes <- grp_sizes[[factor]]
    else grp_sizes <- grp_sizes$stratified[[factor]]
    if(!(min(grp_sizes$Size>2))) {
      stop('There are less than 3 samples in ',
           paste0(names(table(data[[factor]])>2)[(table(data[[factor]])>2)],collapse = ' and '))
    }
    if(!is.null(stratifiers)) {
      for(stratum in stratifiers) {
        data <- data %>% dplyr::group_by(get(stratum),.add=TRUE)
        data <- data %>% dplyr::select(-stratum)
        colnames(data)[ncol(data)] <- stratum
      }
    }
    melted <- data
  }

  stat_results <- calcUniVar(data = melted,
                             factor = factor,
                             stratifiers = stratifiers,
                             features = features,
                             param=param,
                             pairwise_comparisons=pairwise_comparisons)

  # Record the results in the dataset if a dataset was passed to this function
  if(inherits(data,'mvdata')) {
    if(length(stratifiers)) {
      dataset$stats[[factor]][[nameStratification(stratifiers)]]$univar[[rank]] <- stat_results
    } else {
      dataset$stats[[factor]]$univar[[rank]] <- stat_results
    }

    assign('active_dataset',dataset,envir = mvEnv)
    if(!is.null(dataset$name)) assign(dataset$name,dataset,1)

    return(dataset)
  }

  return(stat_results)
}

#' Core Univariate Analysis Function
#'
#' @param data Data table containing metadata and count values with samples as
#'     rows and metadata factors or features as columns
#' @param factor Factor along which to form groups from the samples and perform
#'     the statistical analysis
#' @param stratifiers (Optional) One or two factors to stratify the groups by
#' @param features List of features to analyze. Defaults to all
#' @param param Run parametrized or nonparametrized univariate analysis? Default
#'     is FALSE
#' @param pairwise_comparisons For 3 or more groups, list of vectors of pairwise
#'     comparisons to conduct. Defaults to all possible pairwise comparisons
#'
#' @return List containing a table for overall statistics results and the features
#'     that were skipped due to having the same mean. For analyses with 3 or more
#'     groups, also contains a second table with the pairwise statistics results
#'
calcUniVar <- function(data,factor,stratifiers=NULL,features,
                       param=F,pairwise_comparisons=NULL) {
  # Figure out the number of groups being compared in the main factor
  numgrps <- length(unique(data[[factor]]))
  # Loop over all the features and do the following:
  #   1) Check whether statistical test can be performed on this data
  #   2) If so, perform the appropriate statistical test(s)
  tot_stats <- list()
  pw_stats <- list()
  stats_skipped <- c()
  for(ft in features) {
    ftTab <- data[c(factor,stratifiers,ft)]

    # Fix features that start with a number
    if(is.numeric(type.convert(substr(ft,1,1), as.is=T))) ft <- paste0('`',ft,'`')

    ### Pre-Analaysis Data Check ###
    #------------------------------#
    # Check if any two groups have the same mean in the factor being analyzed
    #   If this is the case, wilcox_test won't work (for now), so we just won't calculate those (for now)
    formula <- as.formula(paste(ft,'~',paste(c(factor,stratifiers),collapse = '+')))
    statdf <- aggregate(formula, ftTab, mean)

    if(is.null(stratifiers)) {
      if(any(duplicated(statdf[[ft]]))) stats_skipped <- c(stats_skipped, ft)
    } else {
      # In case the feature values are integers, convert them to numeric
      #   (which tibbles store as double) so that any value adjustments
      #   that may be made will be retained and not converted back to integer
      ftTab[[gsub('`','',ft)]] <- as.numeric(ftTab[[gsub('`','',ft)]])

      statdf$strata <- interaction(statdf[stratifiers])
      for(stratum_name in unique(statdf$strata)) {
        temp <- statdf[statdf$strata==stratum_name,]
        values <- temp[[ft]]
        duplicates <- values[duplicated(values)]
        for(d in duplicates) {
          dupe_grps <- temp[temp[ft]==d,][[factor]]
          message('\nWARNING: For ',ft,' ',
                  paste(dupe_grps,collapse = ' and '),
                  ' have the same mean value in the "',stratum_name,'" stratum')
          if(length(dupe_grps)<4 & d!=0 & get('forceStats',envir = mvEnv)) {
            value_adjust <- get('value_adjst',envir = mvEnv)
            message('  Adjusting these values by ', value_adjust)
            value_adj <- c(1+value_adjust, 1-value_adjust, 1)
            i <- 1
            for(grp in dupe_grps) {
              ftTab[interaction(ftTab[stratifiers])==stratum_name & ftTab[[factor]]==grp,][ft] <-
                ftTab[interaction(ftTab[stratifiers])==stratum_name & ftTab[[factor]]==grp,][[ft]]*value_adj[i]
              i <- i+1
            }
          } else {
            stats_skipped <- c(stats_skipped,ft)
          }
        }
      }
    }

    if(ft %in% stats_skipped) {
      message('\nWARNING: Skipping statistical testing for ',ft,' since it will fail pairwise analysis of variance with the rstatix package')
      next
    }

    ### Perform the Statistical Testing ###
    #-------------------------------------#
    #     Important note: when adjusting grouped data (ie stratified data)
    #       the p-values ARE adjusted WITHIN their own stratum only,
    #       which is what we want I am pretty sure
    formula <- as.formula(paste(ft,'~',factor))

    if(numgrps > 2) {
      if(param) {#### Parametric analysis of 3 or more groups ####
        tot_stats[[ft]] <- ftTab %>% anova_test(formula)
        if(!(ft %in% stats_skipped)) {
          if(is.null(pairwise_comparisons)) {
            pw_stats[[ft]] <- ftTab %>% rstatix::pairwise_t_test(formula, p.adjust.method='BH') %>%
              rstatix::add_xy_position(x=factor)
            if(get('tukey_games',envir = mvEnv)){
              pvals <- pw_stats[[ft]][, !(names(pw_stats[[ft]]) %in% c('p.adj','p.adj.signif'))]
              padj <- ftTab %>% rstatix::tukey_hsd(formula)
              padj <- padj[c(stratifiers,'group1','group2','p.adj','p.adj.signif')]
              pw_stats[[ft]] <- dplyr::inner_join(pvals, padj, c(stratifiers,'group1','group2')) %>%
                rstatix::add_xy_position(x=factor)
            }
          } else {
            pw_stats[[ft]] <- ftTab %>% rstatix::pairwise_t_test(formula, p.adjust.method='BH',
                                                                 comparisons = pairwise_comparisons) %>%
              rstatix::add_xy_position(x=factor)
          }
        }

      } else {#### Non-parametric analysis of 3 or more groups ####
        tot_stats[[ft]] <- ftTab %>% kruskal_test(formula)
        if(!(ft %in% stats_skipped)) {
          if(is.null(pairwise_comparisons)) {
            pw_stats[[ft]] <- ftTab %>% rstatix::pairwise_wilcox_test(formula, p.adjust.method='BH') %>%
              rstatix::add_xy_position(x=factor)
            if(get('tukey_games',envir = mvEnv)){
              pvals <- pw_stats[[ft]][, !(names(pw_stats[[ft]]) %in% c('p.adj','p.adj.signif'))]
              padj <- ftTab %>% rstatix::games_howell_test(formula)
              padj <- padj[c(stratifiers,'group1','group2','p.adj','p.adj.signif')]
              pw_stats[[ft]] <- dplyr::inner_join(pvals, padj, c(stratifiers,'group1','group2')) %>%
                rstatix::add_xy_position(x=factor)
            }
          } else {
            pw_stats[[ft]] <- ftTab %>% rstatix::pairwise_wilcox_test(formula, p.adjust.method='BH',
                                                                      comparisons = pairwise_comparisons) %>%
              rstatix::add_xy_position(x=factor)
          }
        }

      }
    } else {
      if(param) {#### Parametric analysis of 2 groups ####
        tot_stats[[ft]] <- ftTab %>% rstatix::t_test(formula) %>% rstatix::add_xy_position(x=factor)
      } else {#### Non-parametric analysis of 2 groups ####
        tot_stats[[ft]] <- ftTab %>% rstatix::wilcox_test(formula) %>% rstatix::add_xy_position(x=factor)
      }
    }
  }

  tot_stats <- do.call('rbind',tot_stats)
  pw_stats <- do.call('rbind',pw_stats)

  # If no p-values could be calculated, error out
  if(!length(tot_stats)) stop('ERROR: P-values could not be calculated for any of the features')

  ### Add q-Values to the Overall Statistics of each Feature ###
  #------------------------------------------------------------#
  # Initialize the dataframe using the first feature in tot_stats to get the
  #   right columns (need the right # of columns and the right names).
  #   Then remove all the rows and start fresh (keeping the column headers)
  if(length(stratifiers)) {
    tot_stats$p.adj <- rep(0,nrow(tot_stats))
    for(stratum in unique(interaction(tot_stats[stratifiers]))) {
      stratum1 <- stringr::str_split(stratum,pattern = '\\.')[[1]][1]
      stratum2 <- stringr::str_split(stratum,pattern = '\\.')[[1]][2]
      tot_stats$p.adj[interaction(tot_stats[stratifiers])==stratum] <-
        p.adjust(tot_stats$p[interaction(tot_stats[stratifiers])==stratum],method = 'BH')
    }
  } else tot_stats$p.adj <- p.adjust(tot_stats$p,method = 'BH')
  # Add significance labels
  tot_stats <- tot_stats %>% rstatix::add_significance(p.col = 'p.adj')
  if(!is.null(pw_stats)) pw_stats <- pw_stats %>% rstatix::add_significance(p.col = 'p.adj')

  return(list(stats=tot_stats,
              pw_stats=pw_stats,
              stats_skipped=stats_skipped))
}
