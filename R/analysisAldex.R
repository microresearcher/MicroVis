#' ALDEx2 Differential Abundance Analysis
#'
#' @param dataset MicroVis dataset. Default is the active dataset
#'
#' @return Results of ALDEx2 differential abundance analysis
#' @export
#'
mvaldex <- function(dataset=NULL, rank=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  rank <- rank[rank %in% getRanks(dataset)]
  if(is.null(rank)) rank <- dataset$data$proc$active_rank

  aldex_obj <- makeAldex(dataset,rank = rank)

  numgrps <- nrow(countSamples(dataset, getSizes = T, verbose = F)[[dataset$active_factor]])
  if(numgrps==2) {
    res <- ALDEx2::aldex.ttest(aldex_obj)
    res <- cbind(res, ALDEx2::aldex.effect(aldex_obj))
  } else if(numgrps>2) {
    res <- ALDEx2::aldex.kw(aldex_obj)
  } else stop(dataset$active_factor,' currently has only 1 group')

  dataset$stats[[dataset$active_factor]]$aldex[[dataset$data$proc$active_rank]] <- res

  assign('active_dataset',dataset,envir = mvEnv)
  if(!is.null(dataset$name)) assign(dataset$name,dataset,1)

  return(res)
}
