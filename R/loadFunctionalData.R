#' Load pathways file
#'
#' @param path_to_fxnl (Optional) full path to pathway csv file
#' @param metadata Metadata dataframe
#'
#' @return Loaded pathway abundance dataframe
#' @export
#'
loadFxnlFile <- function(path_to_fxnl=NA,metadata=NULL) {
  cat(paste0('\n\n|~~~~~~~ LOADING FUNCTIONAL DATA ~~~~~~~|\n'))
  fxnlfile <- NULL

  if(!is.null(path_to_fxnl)) if(!file.exists(as.character(path_to_fxnl))) {
    # If path_to_fxnl is set to "NULL" or a file that doesn't exist, user will
    #   be asked to choose a file in the project directory
    message('\nSelect an additional abundance table (csv format) associated with the same metadata. Press "Cancel" or hit "Esc" to skip')
    Sys.sleep(0.1) # To make sure it displays the above message before opening the dialogue box
    fxnlfile <- selectFile(caption='Select additional abundance table (csv format)',
                           path=get('project_dir',envir = mvEnv),
                           filter='Comma-separated Value (*.csv)')
  } else fxnlfile <- path_to_fxnl

  # If a valid functional abundace table was chosen, load it
  if(is.null(fxnlfile)) {
    message('\n No additional abundance table selected')
    return(NULL)
  } else cat('\nAdditional abundance data loading from:\n',fxnlfile)

  fxnl_data <- data.frame(t(read.csv(file.path(fxnlfile),header=FALSE)))
  fxnl_names <- fxnl_data[1,2:ncol(fxnl_data)]

  # List of invalid characters in taxonomy names to replace with an underscore
  invalid_chars <- c('-'='_',' '='_',':'='_','/'='_',
                     '\\['='','\\]'='','\\('='_','\\)'='')
  fxnl_names <- lapply(fxnl_names, function(x) str_replace_all(x,invalid_chars))
  # for(c in invalid_chars) fxnl_names <- lapply(fxnl_names, function(x) str_replace_all(x,c,'_'))

  colnames(fxnl_data)[1] <- c('sample')

  # Remove the row containing original taxa names
  fxnl_data <- data.frame(fxnl_data[2:nrow(fxnl_data),])
  if(!any(is.na(suppressWarnings(as.numeric(fxnl_data$sample))))) {
    # Make the 'sample' column numeric if it can be coerced
    fxnl_data$sample <- as.numeric(fxnl_data$sample)
  }

  # Reorder the rows by the sample number (since our pipeline yields numerical sample names)
  fxnl_data <- fxnl_data %>% arrange(sample)

  sample_names <- fxnl_data$sample

  # Check to make sure all the sample names are in the metadata
  if(!is.null(metadata)) {
    if(!all(sample_names %in% metadata$sample)) {
      missing_samples <- sample_names[!(sample_names %in% metadata$sample)]
      cat('\nThe following sample names were in the additional abundance file but not in the metadata:\n ',paste(missing_samples,collapse='\n  '))
      message('\nSkipping additional data.\n To analyze this additional data, please fix sample names and run "mvLoad()" again\n')
      return(NULL)
    }
  }

  # Remove the 'sample' column and create a numeric dataframe
  fxnl_data <- data.frame(sapply(fxnl_data[,2:ncol(fxnl_data)], function(x) as.numeric(x)))

  # Turn functional names into column names and sample names into row names
  colnames(fxnl_data) <- fxnl_names
  rownames(fxnl_data) <- sample_names

  fxnl_data_list <- list(orig=fxnl_data)
  fxnl_data_list$proc$unranked <- fxnl_data

  fxnl_ranked <- makeRankTabs(fxnl_data_list)

  cat(paste0('\n\n>>> FUNCTIONAL DATA LOADED SUCCESSFULLY <<<\n'))
  return(fxnl_ranked)
}
