#' Load MicroVis Default Settings
#'
#' @return Does not return anything
#' @export
#'
mvparamDefaults <- function() {
  autosave <- F
  offerSave <- T
  autoNameResults <- T
  detailed_taxa_names <- F
  taxaLevels <- c('domain','phylum','class','order','family','genus','species')

  #### DEFAULT COLORS ####
  defCols <- c('#6700b5', # Purple
               '#00a6cf', # Baby blue
               '#d12121', # Red
               '#00a60c', # Green
               '#b500a9', # Pink
               '#d1a400', # Yellow
               '#4432e3', # Blue
               '#db8e00', # Orange
               '#00d18f', # Mint
               '#ed594c', # Coral
               '#310057', # Dark purple
               '#006a85', # Teal
               '#6e1212', # Maroon
               '#005c07', # Dark green
               '#6e0066', # Dark pink?
               '#8c6e00', # Golden
               '#0a0063', # Navy
               '#6b4500', # Brown
               '#006948', # Dark... mint?
               '#8a423b' # Burnt... coral?
  )

  # When forceStats is set to true, values will be slightly adjusted for
  #   groups with the same mean prior to statistical testing
  forceStats <- F

  ### DEFAULT SAMPLE FILTERING CUTOFF
  rthresh <- 10000 # Minimum read count threshold for each sample

  #### DEFAULT FEATURE FILTERING PARAMETERS ####
  filtering.defaults <- list()
  filtering.defaults$min_totabun <- 10
  filtering.defaults$low_abun <- list(min_abun=1,
                                      min_prop=20)
  filtering.defaults$min_relabun <- 0.0001
  filtering.defaults$min_prevalence <- 1
  filtering.defaults$low_var_percentile <- 5

  keepSigFisher <- T

  imageType <- 'png'

  ### ONLY FOR MICROVIS DEVELOPMENT ###
  .debug <- F
  mga <- F
}
