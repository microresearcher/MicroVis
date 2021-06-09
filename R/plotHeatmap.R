#' Plot Heatmap of Data
#'
#' @param dataset MicroVis dataset. Defaults to the active dataset
#' @param clustNum Number of clusters for rows (features)
#' @param ftlist (Optional) List of specific features to plot
#' @param plotSigs Plot only significant features? Defaults to TRUE
#' @param plotUniques Plot only features that are uniquely expressed in a certain
#'     group? Defaults to FALSE
#' @param param Perform parametrized or nonparametrized univariate analysis for
#'     significance testing? Defaults to FALSE
#' @param aggregated Aggregate samples in the same group into a single column?
#'     Defaults to FALSE
#' @param labelFeatures Label rows with the feature names? Defaults to TRUE
#' @param labelSamples Label the columns with the sample names (at the bottom)?
#'     Defaults to FALSE
#' @param factor Factor to perform univariate analysis on
#' @param stratify Stratify the groups by another factor? Defaults to FALSE
#' @param scaleCol Color scaling to improve contrast between certain cells. Value
#'     should be 0 or more. Defaults to 0 (no scaling)
#' @param alpha Significance threshold. Defaults to 0.05
#' @param width Width of saved heatmap figure in inches. Defaults to 6 in
#' @param height Height of saved heatmap figure in inches. Defaults to 8 in
#'
#' @return Heatmap of relative abundances of features in each group
#' @export
#'
plotHeatmap <- function(dataset=NULL,
                        clustNum=0,
                        ftlist=NULL, plotSigs=T, plotUniques=F, param=F,
                        aggregated=F,
                        labelFeatures=T, labelSamples=F,
                        factor=NULL, stratify=F,
                        scaleCol=1,
                        alpha=0.05,
                        width=6,height=12) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  factor <- dataset$active_factor
  rank <- dataset$data$proc$active_rank
  features <- getFeatures(dataset)

  if(plotSigs) {
    if(is.null(dataset$stats[[factor]][[rank]]$univar)) {
      dataset <- univar(dataset,
                        rank=rank,
                        param=param,
                        dataset_name=dataset_name)
    }
    stats <- dataset$stats[[factor]][[rank]]$univar
    sigfts <- unique(stats$stats$.y.[stats$stats$p.adj<=alpha])

    if(length(sigfts)) {
      if(plotUniques) sigfts <- unique(unlist(listUniques(dataset,
                                                          dataset_name=dataset_name)))
      features <- sigfts
      suffix <- paste0('_alpha_',alpha)
    } else {
      message('\nNo significant features were found. Plotting all features')
      suffix <- 'allfts'
    }
  }

  stratifier <- NULL
  if(stratify) {
    stratifier <- select.list(names(dataset$factors)[!(names(dataset$factors) %in% factor)],
                              title = 'Stratify by:')
  }

  # Melt the data to ensure only the samples we are analyzing are in the metadata
  #   and abundance table
  if(!is.null(stratifier)) {
    melted <- data.frame(mvstratify(scaleFeatures(clearNormalization(dataset,temp = T,silent = T),
                                                  scaling = 'relative',temp = T,silent = T)))
  } else {
    melted <- mvmelt(scaleFeatures(clearNormalization(dataset,temp = T,silent = T),
                                   scaling = 'relative',temp = T,silent = T))
  }

  if(aggregated) {
    melted <- melted[c(factor,stratifier,features)]
    formula <- as.formula(paste0('. ~',paste0(c(factor,stratifier),collapse = ' + ')))
    melted <- aggregate(formula,melted,mean)

    hm_data <- t(melted[features])
    # Column names were lost, so reassign the group names to the columns
    #   Group names were in rows before under the factor name but hm_data is transposed
    colnames(hm_data) <- interaction(melted[c(factor,stratifier)])

    # Make a dataframe for annotation factors
    factor_anno <- melted[c(factor,stratifier)]
    suffix <- paste0('_heatmap_aggregated',suffix)
  } else {
    hm_data <- t(melted[features])

    # Column names were lost, so reassign the sample names to the columns
    #   Sample names were row names before, but hm_data is transposed
    colnames(hm_data) <- melted$sample

    # Make a dataframe for annotation factors
    factor_anno <- melted[names(dataset$factors)]
    suffix <- paste0('_heatmap',suffix)
  }
  splitby <- melted[c(factor,stratifier)]

  # Make a colors list from "clrs" that HeatmapAnnotation() will accept
  ha_coloring <- list()
  if(length(dataset$factors)) {
    for(f in names(dataset$factors)) {
      ha_coloring[[f]] <- dataset$colors
    }
  }

  ha <- HeatmapAnnotation(df = factor_anno,
                          col = ha_coloring,
                          show_annotation_name = stratify,
                          show_legend = stratify)

  maxval <- max(hm_data)
  minval <- min(hm_data)
  col_center <- median(hm_data)*(1/scaleCol)

  if(col_center==minval) col_center <- mean(c(maxval,minval))

  colscale <- colorRamp2(c(minval,col_center,maxval),c('#ebebeb','#ffea75','red'))

  hm <- Heatmap(hm_data,
                top_annotation = ha, column_gap = unit(3,'mm'),
                show_column_names = labelSamples,
                column_split = splitby,
                row_km = clustNum, row_km_repeats = 100, cluster_row_slices = FALSE,
                show_row_names = labelFeatures,
                col = colscale,
                heatmap_legend_param = list(title='Relative Abundance',
                                            legend_direction='horizontal',
                                            title_position='topcenter'))

  hm <- draw(hm,heatmap_legend_side='bottom')

  if(labelFeatures) {
    suffix <- paste0(suffix,'_ftslabeled')
    width <- width + 4
  }
  if(labelSamples) {
    suffix <- paste0(suffix,'_sampleslabeled')
    height <- height + 2
  }
  saveResults(dataset$results_path,foldername = 'Heatmaps',
              factors = dataset$factors,
              active_factor = factor,
              figure = hm, width = width, height = height,
              stat_results = stats,
              suffix = suffix)

  cat(paste0('\n  <|> Active Dataset: "',dataset_name,'" <|>\n'))
  print(dataset)

  return(hm)
}
