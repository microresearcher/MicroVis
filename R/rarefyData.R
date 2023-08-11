#' Rarefy Dataset
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp This isn't needed
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with rarefied data.
#'
runRarefaction <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  if(is.null(dataset$name)) dataset_name <- 'active_dataset'
  else dataset_name <- dataset$name

  if(is.null(dataset$data$proc$rarefied)) {
    if(!silent) cat('~~~ No rarefaction performed ~~~\n')
    return(dataset)
  }

  if(!silent) cat('|~~~~~~~~~~~~~  RARIFYING SAMPLES  ~~~~~~~~~~~~~~|\n')

  ft_data <- dataset$data
  abd <- ft_data$proc$unranked

  # Show the rarefaction curves before rarefying
  p_rcbefore <- plotRareCurves(dataset, getPlot=T)
  if(!silent) show(p_rcbefore)

  # Get the minimum richness of all the samples
  minrich <- min(rowSums(abd))

  # If there are samples with 0 reads (and therefore 0 richness), do not rarefy
  if(minrich==0) {return(message('\nERROR: At least one sample in this dataset has 0 reads\n'))}

  # Rarefy the abundance table
  abd.rarefied <- vegan::rrarefy(abd, minrich)
  # Now make ranked abundance tables with this abundance table
  ft_data$proc$unranked <- abd.rarefied
  ft_data <- makeRankTabs(ft_data)

  # Reset the "unranked" table in "ft_data$proc" so that it is unrarefied
  # ft_data$proc$unranked <- dataset$data$proc$unranked
  dataset$data <- ft_data

  # Show the rarefaction curves after rarefying
  p_rcafter <- plotRareCurves(dataset, getPlot=T)
  if(!silent) show(p_rcafter)

  dataset$data$proc$rarefied <- minrich

  if(!silent) cat('\n>>> Rarefied samples to',minrich,'reads <<<\n')

  return(dataset)
}

#' Data Rarefaction
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param temp This is not needed
#' @param silent Argument that is ultimately passed onto runSampleFilter(),
#'     runNormalization(), and runFeatureFilter(), telling them not to output
#'     any messages.
#'
#' @return MicroVis dataset (mvdata object) with updated rarefaction parameters
#' @export
#'
rarefySamples <- function(dataset=NULL, temp=F, silent=F) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  dataset$data$proc$filtering <- NULL
  dataset$data$proc$normalization <- NULL
  if(is.null(dataset$data$proc$unranked)) dataset <- runSampleFilter(dataset, silent=T)

  abd <- dataset$data$proc$unranked

  # Get the minimum richness of all the samples
  minrich <- min(rowSums(abd))

  # If there are samples with 0 reads (and therefore 0 richness), do not rarefy
  if(minrich==0) {return(message('\nERROR: At least one sample in this dataset has 0 reads\n'))}

  dataset$data$proc$rarefied <- minrich

  # Prepare the dataset for use
  dataset <- processDataset(dataset, silent = silent)

  return(dataset)
}
