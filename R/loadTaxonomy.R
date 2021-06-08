#' Load Taxonomy File
#'
#' @param path_to_taxa (Optional) full path to taxonomy csv file
#' @param metadata Metadata dataframe
#' @param combineDupes If set to TRUE (default), MicroVis will try to combine
#'     duplicate features.
#'
#' @return Loaded taxonomy data - a list of:
#'     1) orig: A sample-by-ASV (rows-by-columns) original abundance dataframe
#'     that will not be processed.
#'     2) A temporary copy of orig that will be used in downstream processing
#'     3) An ASV-by-rank dataframe with ASVs and their assigned taxa at each
#'     taxonomic rank
#'
loadTaxaFile <- function(path_to_taxa=NA,metadata=NULL,combineDupes=T) {
  cat(paste0('\n\n|~~~~~~~~ LOADING TAXONOMY DATA ~~~~~~~~|\n'))
  taxafile <- NULL

  if(!is.null(path_to_taxa)) if(!file.exists(as.character(path_to_taxa))) {
    # If path_to_taxa is set to "NA" or a file that doesn't exist, user will be asked to choose a file in the project directory
    message('\nSelect taxonomy abundance table (csv format). Press "Cancel" or hit "Esc" to skip')
    Sys.sleep(0.1) # To make sure it displays the above message before opening the dialogue box
    taxafile <- selectFile(caption='Select taxonomy abundance table (csv format)',
                           path=get('project_dir',envir = mvEnv),
                           filter='Comma-Separated Value (*.csv)')
  } else taxafile <- path_to_taxa

  # If a valid taxonomy abundace table was chosen, load it
  if(is.null(taxafile)) {
    message('\n No taxonomy abundance table selected')
    return(NULL)
  } else cat('\nTaxonomy abundance data loading from:\n',taxafile)

  # List of invalid characters in taxonomy names to replace with an underscore
  invalid_chars <- c('-'='_',' '='_',':'='_','/'='_',
                     '\\['='','\\]'='','\\('='_','\\)'='')

  taxa_data <- read.csv(file.path(taxafile),header=FALSE)

  rank_cols <- c()
  for(rank in get('taxaRanks',envir = mvEnv)) rank_cols <- c(rank_cols,grep(rank, taxa_data[1,], ignore.case = T))
  if(length(rank_cols)) {
    taxa_ranks <- tolower(unname(unlist(taxa_data[1,][tolower(taxa_data[1,]) %in%
                                                        tolower(taxa_data[rank_cols][1,])])))
    taxa_names_tab <- taxa_data[2:nrow(taxa_data),rank_cols]
    colnames(taxa_names_tab) <- taxa_ranks
    rownames(taxa_names_tab) <- paste(1:nrow(taxa_names_tab))

    taxa_names_tab <- data.frame(apply(taxa_names_tab,2,function(x) str_replace_all(x,invalid_chars)))
    # for(c in invalid_chars) taxa_names_tab <- data.frame(apply(taxa_names_tab, 2, function(x) str_replace_all(x,c,'_')))

    taxa_names_tab <- cleanASVs(taxa_names_tab)

    # Remove the rank columns from the abundance table
    taxa_data <- taxa_data[,-rank_cols]
  } else {
    taxa_names_tab <- taxa_data[2:nrow(taxa_data),1]

    taxa_names_tab <- makeASVtab(taxa_names_tab)

    # Remove the rank column from the abundance table
    taxa_data <- taxa_data[,2:ncol(taxa_data)]
  }
  # Now that taxa names have been cleaned up and - if necessary - split by taxa
  #   levels, we can assign ASVs to the rows in both the taxa names table and
  #   the abundance table
  rownames(taxa_names_tab) <- paste0('ASV_',1:nrow(taxa_names_tab))

  # First turn the first row of the abundance table into column names (assuming
  #   that these are sample names)
  colnames(taxa_data) <- sample_names <- unname(unlist(taxa_data[1,]))
  taxa_data <- taxa_data[2:nrow(taxa_data),]

  rownames(taxa_data) <- rownames(taxa_names_tab)

  taxa_data <- data.frame(t(taxa_data))
  if(length(unique(sample_names))!=length(sample_names)) {
    dupe_samples <- sample_names[duplicated(sample_names)]
    message('\nThe following sample names were duplicated:\n',
            paste(sample_names[sample_names %in% dupe_samples], collapse='\t'),
            '\n\nThey have been changed to:')
    for(dupe in unique(dupe_samples)) {
      dupe_ind <- which(sample_names==dupe)
      for(i in seq_along(dupe_ind)) sample_names[dupe_ind[i]] <- paste(sample_names[dupe_ind[i]], i, sep='_')
      message(paste(sample_names[dupe_ind], collapse='\t'))
    }
  }
  rownames(taxa_data) <- sample_names

  # Clean up taxonomy table and order the samples (stored in rows) by sample number
  taxa_data$sample <- rownames(taxa_data)
  if(!any(is.na(suppressWarnings(as.numeric(taxa_data$sample))))) {
    # Make the 'sample' column numeric if it can be coerced
    taxa_data$sample <- as.numeric(taxa_data$sample)
  }
  # Reorder the rows by the sample number (since our pipeline yields numerical sample names)
  taxa_data <- taxa_data %>% arrange(sample)

  # If there is a metadata, check to make sure sample names in abundance table
  #   are the same
  if(!is.null(metadata)) {
    if(!all(taxa_data$sample %in% metadata$sample)) {
      missing_samples <- taxa_data$sample[!(taxa_data$sample %in% metadata$sample)]
      cat('\nThe following sample names were in the taxonomic abundance file but not in the metadata:\n',paste(missing_samples,collapse='\t'))
      message('\nSkipping taxonomic data.\n To analyze taxonomic data, please fix sample names and run "mvLoad()" again\n')
      return(NULL)
    }
  }
  # Store sample names in their new order
  sample_names <- taxa_data$sample
  # Remove the 'sample' column
  taxa_data$sample <- NULL
  # Create a numeric dataframe
  taxa_data <- data.frame(sapply(taxa_data, function(x) as.numeric(x)))
  rownames(taxa_data) <- sample_names
  colnames(taxa_data) <- rownames(taxa_names_tab)

  taxa_data_list <- list(orig=taxa_data,
                         taxa_names=taxa_names_tab)

  # Now identify identical ASVs and combine them in the taxa table
  taxa_data_list <- combineDupeASVs(taxa_data_list,combineDupes=combineDupes)

  # Make another copy of the abundance table labeled "unranked"
  #   This serves as a placeholder for when the abundance tables are actually
  #   processed (e.g. samples or taxa are filtered)
  taxa_data_list$proc$unranked <- taxa_data_list$orig

  taxa_ranked <- makeRankTabs(taxa_data_list)

  cat(paste0('\n>>> TAXONOMY DATA LOADED SUCCESSFULLY <<<\n'))

  return(taxa_ranked)
}

#' Rename Taxa with NA and at Species Level
#'
#' @param dataset Unprocessed MicroVis dataset (mvdata object) that still contains
#'     the unranked abundance table.
#'
#' @return MicroVis dataset (mvdata object) with cleaned taxonomy names in its
#'     asv-taxonomy reference table
#' @export
#'
cleanUnkTaxa <- function(dataset=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  taxa_names_tab <- dataset$data$taxa_names

  taxa_names_tab <- data.frame(apply(taxa_names_tab,2,function(x) {
    unks <- grep('unidentified_.*_of_',x)
    x[unks] <- paste0(gsub('_[^_]*$','',gsub('unidentified_.*_of_','',x[unks]),perl=T),
                      '_',substr(gsub('.*_','',x[unks]),1,1),
                      gsub('unidentified','',gsub('_of_.*','',x[unks])))
    return(x)
  }))

  dataset$data$taxa_names <- taxa_names_tab

  return(dataset)
}
