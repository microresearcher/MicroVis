#' Load metadata file
#'
#' @param path_to_metadata (Optional) full path to metadata csv file
#'
#' @return Loaded metadata dataframe
#'
loadMDFile <- function(path_to_metadata=NA) {
  cat(paste0('\n\n|~~~~~~~~ LOADING metadata FILE ~~~~~~~~|\n'))
  mdfile <- NULL

  if(!is.null(path_to_metadata)) if(!file.exists(as.character(path_to_metadata))) {
    path_to_metadata <- NA

    message('\nSelect metadata table (csv format)')
    Sys.sleep(0.1) # To make sure it displays the above message before opening the dialogue box

    mdfile <- rstudioapi::selectFile(caption='Select metadata table (csv format)',
                                     path=get('project_dir',envir = mvEnv),
                                     filter="Comma-Separated Value (*.csv)")

  } else mdfile <- path_to_metadata

  if(is.null(mdfile)) {
    message('\nERROR: No metadata file selected. No comparative analysis can be run !!!\n')
    return(NULL)
  } else cat('\nmetadata loading from:\n',mdfile,'\n')

  metadata <- read.csv(file.path(mdfile))

  # Sample name must be in first column of metadata file!!
  if(length(unique(metadata[[1]])) != length(metadata[[1]])) {
    assign('.loading',FALSE,envir = mvEnv)
    return(message('\nERROR: First column of metadata must have sample names which should all be unique'))
  }

  colnames(metadata)[1] <- 'sample'

  # Excel writes NA values as "#N/A"
  metadata[metadata == '#N/A'] <- NA

  # Convert columns into numeric/integer if they are continuous data
  vars.numeric <- sapply(type.convert(metadata, as.is = T), is.numeric)
  vars.numeric <- names(vars.numeric[vars.numeric])
  metadata[vars.numeric] <- lapply(vars.numeric, function(v) {
    if(select.list(c('Yes','No'),
                   title = paste('\nIs',v,'a continuous variable?'))=='Yes') {
      type.convert(metadata[[v]], as.is = T)
    } else as.character(metadata[[v]])
  })

  cat(paste0('\n>>> METADATA FILE LOADED SUCCESSFULLY <<<'))

  return(metadata)
}
