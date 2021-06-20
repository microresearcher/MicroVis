#' Get samples with certain metadata characteristics from a dataset
#'
#' @param data MicroVis dataset. Defaults to the active dataset
#' @param id_cols Vector of names of metadata columns that uniquely identify subjects
#' @param filter Vector of names of metadata columns that each subject must have all of
#'
#' @return Data table with just the samples identifiedy by "id_cols" that have
#'     all the values in filter
#' @export
#'
getSamples <- function(dataset=NULL, id_cols=NULL, complete=NULL, filter=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  data <- mvmelt(dataset)[1:ncol(dataset$metadata)]

  if(is.null(id_cols)) id_cols <- 'sample'

  id_cols <- colnames(data)[colnames(data) %in% id_cols]
  if(!length(id_cols)) {
    message(paste(id_cols,collapse=', '),' must be column names in the data')
    id_cols <- 'sample'
  }

  samples <- samplesWith(data, id_cols=id_cols,
                         complete=complete,
                         filter=filter)[unique(c('sample',id_cols))]

  return(samples)
}

#' Filter a datatable to only samples with certain metadata characteristics
#'
#' @param data Data table with metadata
#' @param id_cols Vector of names of columns that uniquely identify subjects
#' @param filter (Optional) Vector of names of columns that each subject must have all of
#'
#' @return Data table with just the samples identified by "id_cols" that have
#'     all the values in filter
#' @export
#'
samplesWith <- function(data, id_cols, complete=NULL, filter=NULL) {
  id_cols <- colnames(data)[colnames(data) %in% id_cols]
  if(!length(id_cols)) stop(paste(id_cols,collapse=', '),' must be column names in the data')

  if(!length(c(complete,filter))) return(data[unique(c('sample',id_cols))])

  filter <- filter[names(filter) %in% colnames(data)]

  for(f in names(filter)) {
    filter[[f]] <- filter[[f]][filter[[f]] %in% as.character(unique(data[[f]]))]
    if(!length(filter[[f]])) filter[[f]] <- NULL
    else data <- data[data[[f]] %in% filter[[f]],]
  }
  filter.list <- filter
  for(c in complete) filter.list[[c]] <- 'all'

  if(!length(filter.list)) stop(paste(c(complete,names(filter)),collapse=', '),
                                ' must be column names in the data')

  md <- data[c(id_cols, names(filter.list))]
  md$count <- rep(1,nrow(md))

  pivoted <- md %>% tidyr::pivot_wider(id_cols=id_cols,
                                       names_from=names(filter.list),
                                       values_from='count')

  filtered <- pivoted[complete.cases(pivoted),]

  filtered_data <- merge(data,filtered[id_cols],by=id_cols)

  return(filtered_data)
}
