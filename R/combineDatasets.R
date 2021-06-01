#' Merge Two Datasets with Different Sample Sets
#'
#' @param path_to_folder Folder with project directory.
#'
#' @return Merged dataset with samples from both datasets kept separate and all
#'     the features from both datasets.
#'
mvcombine <- function(path_to_folder=NA) {
  if(is.null(path_to_folder)) {path_to_folder <- NA}
  if(!dir.exists(as.character(path_to_folder))) {
    message('\nSelect directory with abundance data and metadata. Results will also be saved to this directory\n')
    Sys.sleep(0.1) # To make sure it displays the above message before opening the dialogue box
    project_dir <- selectDirectory('Select directory with abundance data and metadata',
                                   path=dirname(getwd()))
  } else {
    project_dir <- path_to_folder
  }
  metadata <- list()
  taxa <- list()
  rank_cols <- list()
  pwys <- list()
  ds_names <- list()
  combined_paths <- list()

  # Start loading the first dataset
  # First load metadata
  message('\nSelect metadata file for first dataset')
  Sys.sleep(0.1)
  metadata$a <- read.csv(selectFile(path = file.path(project_dir)))
  colnames(metadata$a)[1] <- 'sample'

  # Then load taxonomy tables
  message('\nSelect taxonomy file for first dataset')
  Sys.sleep(0.1)
  taxa.a_path <- selectFile(path = file.path(project_dir))
  if(!is.null(taxa.a_path)) {
    taxa$a <- read.csv(taxa.a_path,header=F)
    rank_cols$a <- c()
    for(rank in get('taxaranks',envir = mvEnv)) rank_cols$a <- c(rank_cols$a,grep(rank, taxa$a[1,],
                                                                                  ignore.case = T))
  }

  # Then load pathways table
  message('\nSelect pathways file for first dataset')
  Sys.sleep(0.1)
  pwys.a_path <- selectFile(path = file.path(project_dir))
  if(!is.null(pwys.a_path)) {
    pwys$a <- read.csv(pwys.a_path,header=F)
  }

  # Now get a name for the first dataset
  ds_names$a <- readline(prompt = cat('Provide a short name for this dataset:\n'))
  metadata$a$Dataset <- rep(ds_names$a,nrow(metadata$a))

  # Start loading the second dataset
  # First load metadata
  message('\nSelect metadata file for second dataset')
  Sys.sleep(0.1)
  metadata$b <- read.csv(selectFile(path = file.path(project_dir)))
  colnames(metadata$b)[1] <- 'sample'

  # Then load taxonomy table if one was loaded for the first dataset
  if(!is.null(taxa.a_path)) {
    message('\nSelect taxonomy file for second dataset')
    Sys.sleep(0.1)
    taxa.b_path <- selectFile(path = file.path(project_dir))
    taxa$b <- read.csv(taxa.b_path,header=FALSE)
    rank_cols$b <- c()
    for(rank in get('taxaranks',envir = mvEnv)) {
      rank_cols$b <- c(rank_cols$b,grep(rank, taxa$b[1,], ignore.case = T))
    }
  }

  # Then load pathways table, if one was loaded for the first dataset
  if(!is.null(pwys.a_path)) {
    message('\nSelect pathways file for second dataset')
    Sys.sleep(0.1)
    pwys.b_path <- selectFile(path = file.path(project_dir))
    pwys$b <- read.csv(pwys.b_path,header=F)
  }

  # Now get a name for the second dataset
  ds_names$b <- readline(prompt = cat('Provide a short name for this dataset:\n'))
  metadata$b$Dataset <- rep(ds_names$b,nrow(metadata$b))

  # Rename the samples in each of the metadata file
  metadata$a$sample <- paste0(metadata$a$sample,'_',ds_names$a)
  metadata$b$sample <- paste0(metadata$b$sample,'_',ds_names$b)

  # Only keep the columns that are shared in the metadatas
  if(!identical(colnames(metadata$a),colnames(metadata$b))) {
    shared_columns <- colnames(metadata$a)[colnames(metadata$a) %in% colnames(metadata$b)]
    metadata$a <- metadata$a[shared_columns]
    metadata$b <- metadata$b[shared_columns]
  }

  project_dir <- file.path(project_dir,paste0(ds_names$a,'_',ds_names$b,'_combined'))
  dir.create(project_dir,showWarnings = F)
  combined_paths$project <- project_dir

  # Combine the metadatas
  metadata <- rbind(metadata$a,metadata$b)
  combined_paths$metadata <- file.path(project_dir,paste0(ds_names$a,'_',ds_names$b,'_combined-metadata.csv'))
  write.csv(metadata,file = combined_paths$metadata,row.names = F)

  invalid_chars <- c('-'='_',' '='_',':'='_','/'='_',
                     '\\['='','\\]'='','\\('='_','\\)'='')

  if(!length(rank_cols$a) & !length(rank_cols$b)) {
    if(length(taxa)) for(tab in 1:length(taxa)) {
      colnames(taxa[[tab]]) <- c('taxa',paste0(taxa[[tab]][1,2:ncol(taxa[[tab]])],'_',ds_names[[tab]]))
      taxa[[tab]] <- taxa[[tab]][2:nrow(taxa[[tab]]),]
    }
    if(length(pwys)) for(tab in 1:length(pwys)) {
      colnames(pwys[[tab]]) <- c('pwys',paste0(pwys[[tab]][1,2:ncol(pwys[[tab]])],'_',ds_names[[tab]]))
      pwys[[tab]] <- pwys[[tab]][2:nrow(pwys[[tab]]),]
    }
  } else if(rank_cols$a != rank_cols$b) {
    message('\nTaxonomy names need to be in the same format in both datasets')
    return(NULL)
  }

  if(length(taxa)) {
    taxa_data <- full_join(taxa$a,taxa$b)
    taxa_data[is.na(taxa_data)] <- 0
    combined_paths$taxa <- file.path(project_dir,paste0(ds_names$a,'_',ds_names$b,'_combined-taxa.csv'))
    write.csv(taxa_data,file = combined_paths$taxa,row.names = F)
  }
  if(length(pwys)) {
    pwys_data <- full_join(pwys$a,pwys$b)
    pwys_data[is.na(pwys_data)] <- 0
    combined_paths$pwys <- file.path(project_dir,paste0(ds_names$a,'_',ds_names$b,'_combined-pwys.csv'))
    write.csv(taxa_data,file = combined_paths$pwys,row.names = F)
  }

  return(combined_paths)
}
