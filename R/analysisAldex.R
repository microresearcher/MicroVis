#' ALDEx2 Differential Abundance Analysis
#'
#' @param dataset MicroVis dataset. Default is the active dataset
#'
#' @return Results of ALDEx2 differential abundance analysis
#' @export
#'
mvaldex <- function(dataset=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  aldex_obj <- makeAldex(dataset)

  numgrps <- nrow(countSamples(dataset, getSizes = T, verbose = F)[[dataset$active_factor]])
  if(numgrps==2) {
    res <- aldex.ttest(aldex_obj)
    res <- cbind(res,aldex.effect(aldex_obj))
  } else if(numgrps>2) {
    res <- aldex.kw(aldex_obj)
  } else stop(dataset$active_factor,' currently has only 1 group')

  dataset$stats[[dataset$active_factor]]$aldex[[dataset$data$proc$active_rank]] <- res

  assign('active_dataset',dataset,envir = mvEnv)
  if(!is.null(dataset$name)) assign(dataset$name,dataset,1)

  return(res)
}
