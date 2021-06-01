#' List available functions
#'
#' @param pages (Optional) Specify what function options to see with "filtering",
#'     "normalizing", "analysis", "plotting", or "all"
#'
#' @return Prints a list of available functions with short descriptions
#' @export
#'
mvhelp <- function(pages=c('filtering','normalizing','analysis','plotting')) {
  if(any(c('processing','normalizing') %in% pages)) {
    cat(
      '\nNormalization functions currently available:\n
    - scaleSamples()    scale samples by total sum, median, range, standardized,
                        or centered\n
    - scaleFeatures()   scale features by total sum, median, range, standardized,
                        or centered\n
    - transformData()   transform data by log, pseudolog, or generalized log\n')
  }

  if(any(c('processing','filtering') %in% pages)) {
    cat(
      '\nFiltering functions currently available:\n
    - filterLowPrev       filter features with overall prevalence below a cutoff\n
    - filterLowRelAbun()  filter features with relative abundance below a cutoff\n
    - filterLowAbun()     filter features with abundance lower than a cutoff in
                          more than a certain proportion of features\n
    - filterLowTotAbun()  filter features with overall abundance below a cutoff\n
    - filterLowVar()      filter features with variance falling below a certain
                          percentile\n')
  }

  if('analysis' %in% pages) {
    cat(
      '\nAnalysis functions currently available:\n
    - univar()            perform univariate statistical analysis\n
    - mvlefse()           perform LEfSe analysis\n
    - mvdeseq()           perform DESeq2 analysis\n')
  }

  if('plotting' %in% pages) {
    cat(
      '\nPlotting functions currently available:\n
    - plotAlpha()         plot boxplots of alpha diversity\n
    - plotBeta()          plot a PCoA plot using a user-specified metric\n
    - plotStackedBars()   plot stacked abundance barplot of features\n
    - plotHM()            plot a heatmap of significant features\n
    - plotSimilarity()    plot heatmap that clusters samples by varying
                          distance and clustering methods\n
    - plotFtCor()         plot a correlation heatmap between samples of 1 dataset
                          or between features of 2 datasets\n
    - plotChords()        plot chord diagram of feature correlations\n
    - plotRegression()    plot regression line between two features within 1 dataset
                          or between 2 datasets
    - plotBars()          plot a grouped relative abundance bar plot of significant
                          features in every taxonomic rank\n
    - plotBox()           plot the abundance boxplots of significant features\n
    - plotRF()            plot important features determined by the Boruta algorithm\n
    - plotLEFSE()         plot an effect sizes vertical bar chart of enriched
                          taxa determined by LEfSe analysis\n
    - plotClad()          plot a cladogram of enriched taxa determined by LEfSe
                          analysis\n
    \n')
  }
}
