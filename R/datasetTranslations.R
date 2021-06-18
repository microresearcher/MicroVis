#' Make a phyloseq object from a MicroVis dataset
#'
#' @param dataset MicroVis dataset. Default is the active dataset
#'
#' @return Phyloseq object
#' @export
#'
makePS <- function(dataset=NULL) {
  # https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- dataset$name
  }

  if(dataset$features!='taxa') stop('Phyloseq is intended for taxonomic data only')

  taxaranks <- get('taxaRanks',envir = mvEnv)

  if(!(all(colnames(dataset$data$taxa_names) %in% taxaranks))) stop('Taxonomic ranks must be one or more of:\n   ',paste0(taxaranks,collapse = ', '))

  taxa_data <- dataset$data

  rank <- getLowestRank(dataset)
  abd <- taxa_data$proc[[rank]]
  abd$Other <- NULL
  colnames(abd) <- TaxatoASV(taxa_data,taxalist = colnames(abd), taxa_rank = rank)
  metadata <- data.frame(lapply(dataset$metadata,function(x) if(is.factor(x)) unfactor(x) else x))

  sample_names <- rownames(abd)

  taxa_df <- taxa_data$taxa_names

  colnames(taxa_df)[grep('domain',colnames(taxa_df))] <- 'kingdom'
  colnames(taxa_df) <- unlist(lapply(colnames(taxa_df),function(x) capitalize(x)))

  samples_df <- data.frame(metadata[metadata$sample %in% sample_names, 2:ncol(metadata)])
  colnames(samples_df) <- colnames(metadata)[2:ncol(metadata)]
  rownames(samples_df) <- metadata[metadata$sample %in% sample_names, 'sample']

  abun_mat <- as.matrix(abd)
  taxa_mat <- as.matrix(taxa_df)

  abd <- otu_table(abun_mat, taxa_are_rows = FALSE)
  tax <- tax_table(taxa_mat)
  samples <- sample_data(samples_df)

  phylodata <- phyloseq(abd,tax,samples)

  # Make tree and add to phyloseq object
  tree <- rtree(ntaxa(phylodata), rooted=TRUE, tip.label=taxa_names(phylodata))
  phylodata <- merge_phyloseq(phylodata,tree)
  if(exists('tree')) {cat('\nPhylogenetic tree created successfully!\n')}

  return(phylodata)
}

#' Make a DESeq object from a MicroVis dataset
#'
#' @param dataset MicroVis dataset. Default is the active dataset
#' @param baseline A reference group for building the DESeq dataset. If
#'     none is selected the user will be asked to select one during function run.
#'
#' @return DESeq object
#' @export
#'
makeDeseq <- function(dataset=NULL,baseline=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- dataset$name
  }

  dataset <- clearNormalization(dataset, temp = T, silent = T)

  active_rank <- dataset$data$proc$active_rank
  abd <- dataset$data$proc[[active_rank]]
  abd$Other <- NULL

  metadata <- dataset$metadata[as.character(dataset$metadata$sample) %in% rownames(abd),]
  metadata <- metadata %>% arrange('sample')

  samples <- rownames(abd)
  fts <- colnames(abd)
  abd <- data.frame(t(abd))
  rownames(abd) <- fts
  colnames(abd) <- samples
  abd <- abd %>% select(as.character(metadata$sample))

  factor <- dataset$active_factor
  baseline <- baseline[baseline %in% dataset$factors[[factor]]$subset]

  cts <- abd
  coldata <- metadata

  if(is.null(baseline)) {
    cat('\nSelect a baseline group\n')
    baseline <- select.list(unique(as.character(coldata[[factor]])),
                            title = 'Select a baseline',graphics = T)
  }

  if(is.null(baseline)) {
    baseline <- levels(coldata[[factor]])[1]
    message('\nWarning: No baseline chosen so ',baseline,' was automatically selected')
  }
  coldata[[factor]] <- relevel(coldata[[factor]],baseline)

  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = as.formula(paste('~',factor)))
  return(dds)
}

#' Make a taxmap object from a MicroVis dataset
#'
#' @param dataset MicroVis dataset
#' @param unfiltered Use the unfiltered MicroVis dataset?
#' @param ftlist List of features to include into the taxmap object
#'
#' @return taxmap object
#' @export
#'
makeTaxMap <- function(dataset=NULL,unfiltered=F,ftlist=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- dataset$name
  }

  if(dataset$features!='taxa') stop('"taxmap" objects are for taxonomic data only')

  taxa_data <- dataset$data$taxa_names
  taxa_data$ASV <- rownames(taxa_data)

  md <- dataset$metadata

  if(unfiltered) {
    # The original abundance table already has the ASV#'s instead of taxa names
    abun <- dataset$data$orig
  } else {
    lowest_rank <- getLowestRank(dataset)

    ftlist <- ftlist[ftlist %in% getFeatures(dataset)]
    if(is.null(ftlist)) ftlist <- getFeatures(dataset, ranks=lowest_rank)

    abun <- dataset$data$proc[[lowest_rank]][ftlist]
    abun$Other <- NULL

    # Need to convert taxa names in processed abundance tables to ASV#'s
    colnames(abun) <- TaxatoASV(taxa_data = dataset$data,
                                taxalist = colnames(abun),
                                taxa_rank = lowest_rank)

    # Subset the taxa_data to only those taxa being analyzed if filtered is TRUE
    taxa_data <- subset(taxa_data,ASV %in% colnames(abun))
  }
  sample_names <- rownames(abun)

  # Subset the metadata so that only the samples being analyzed are included
  md <- md[md$sample %in% sample_names,]

  abun <- cbind(ASV=colnames(abun),data.frame(t(abun)))
  colnames(abun)[2:ncol(abun)] <- sample_names

  mvtaxmap <- parse_tax_data(taxa_data,class_cols = 1:(ncol(taxa_data)-1),named_by_rank = T,
                             datasets = list(abundance=abun),mappings = c('ASV'='ASV'))
  mvtaxmap$data$abundance <- calc_taxon_abund(mvtaxmap,'abundance',cols=sample_names)
  mvtaxmap$data$metadata <- md

  return(mvtaxmap)
}
