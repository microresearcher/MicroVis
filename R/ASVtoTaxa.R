#' Translate ASV Numbers to Taxa Names
#'
#' @param taxa_data Taxonomy data (list object) from a dataset containing
#'     asv-taxonomy reference table.
#' @param asvlist List of ASV numbers to look up in taxa_data.
#' @param taxa_rank (Optional) The desired rank at which to look up ASV numbers.
#'     Will default to the active rank if not specified.
#' @param uniqueOnly (Optional) If multiple taxa map to the same parent taxon,
#'     only return the name of that parent taxon once. Default is FALSE
#'
#' @return List of taxa corresponding to the ASVs in the same order.
#' @export
#'
ASVtoTaxa <- function(taxa_data, asvlist, taxa_rank=NULL, uniqueOnly=F) {
  taxa_names_tab <- taxa_data$taxa_names
  if(is.null(taxa_rank)) {
    taxa_rank <- taxa_data$proc$active_rank
  }

  if(ncol(taxa_names_tab)==1) taxanames <- taxa_names_tab[rownames(taxa_names_tab) %in% asvlist,]
  else taxanames <- taxa_names_tab[[taxa_rank]][rownames(taxa_names_tab) %in% asvlist]

  if(uniqueOnly) taxanames <- unique(taxanames)

  others_index <- grep('Other',asvlist)
  if(length(others_index)) taxanames[[others_index]] <- 'Other'

  return(taxanames)
}

#' Translate Taxa Names to ASV Numbers
#'
#' @param taxa_data Taxonomy data (list object) from a dataset containing
#'     asv-taxonomy reference table.
#' @param taxalist List of taxa names to look up in taxa_data.
#' @param taxa_rank (Optional) The desired rank at which to look up taxa names.
#'     Will default to the active rank if not specified.
#'
#' @return List of ASVs corresponding to the taxa in the same order.
#' @export
#'
TaxatoASV <- function(taxa_data, taxalist, taxa_rank=NULL) {
  taxa_names_tab <- taxa_data$taxa_names
  if(is.null(taxa_rank)) {
    taxa_rank <- taxa_data$proc$active_rank
  }

  for(i in 1:length(taxalist)) if(startsWith(taxalist[i],'X') & is.numeric(type.convert(substr(taxalist[i],2,2)))) taxalist[i] <- sub('X','',taxalist[i])

  if(ncol(taxa_names_tab)==1) asvs <- rownames(taxa_names_tab)[taxa_names_tab[[taxa_rank]] %in% taxalist]
  else asvs <- rownames(taxa_names_tab)[taxa_names_tab[[taxa_rank]] %in% taxalist]

  if(length(asvs)!=length(taxalist)) {
    message('\nWarning: Some taxa names at the ',taxa_rank,
            ' level mapped to more than one ASV.\n This means there are some duplicate taxa names at the ',
            taxa_rank,' level\n')
  }

  others_index <- grep('Other',taxalist)
  if(length(others_index)) asvs[[others_index]] <- 'Other'

  return(asvs)
}

#' Agglomerate taxa into parent taxa
#'
#' @param taxa_data Taxonomy data with abundance tables and taxa names table
#' @param abundance_table Abundance table
#' @param from_rank Starting rank (must be lower than "to_rank")
#' @param to_rank Destination rank (must be higher than "from_rank")
#'
#' @return List of parent taxa at the "to_rank" rank for each taxon in "abundance_table"
#' @export
#'
agglomTaxa <- function(taxa_data, abundance_table, from_rank, to_rank) {
  taxa_names_tab <- taxa_data$taxa_names

  if(!(tolower(to_rank) %in% colnames(taxa_names_tab))) {
    stop(paste0(to_rank,' is not a rank in this dataset'))
  }
  if(!(tolower(from_rank) %in% c('asv',colnames(taxa_names_tab)))) {
    stop(paste0(from_rank,' is not a rank in this dataset'))
  }

  if(tolower(from_rank)=='asv') {
    colnames(abundance_table) <- ASVtoTaxa(taxa_data,
                                           colnames(abundance_table),
                                           taxa_rank = tolower(to_rank))
  } else {
    colnames(abundance_table) <- taxa_names_tab[[to_rank]][taxa_names_tab[[from_rank]] %in% colnames(abundance_table)]
  }

  abundance_table <- data.frame(t(rowsum(t(abundance_table),group = colnames(abundance_table))))

  for(i in 1:ncol(abundance_table)) {
    if(startsWith(colnames(abundance_table)[i],'X') & is.numeric(type.convert(substr(colnames(abundance_table)[i],2,2)))) {
      colnames(abundance_table)[i] <- sub('X','',colnames(abundance_table)[i])
    }
  }

  return(abundance_table)
}

#' Get parent taxa
#'
#' @param taxalist Vector of current taxa
#' @param from Rank of current taxa
#' @param to Rank of desired parent taxa
#' @param taxanames Table of all taxa names
#'
#' @return Vector of parent taxa for each taxa in "taxalist"
#' @export
#'
getParentTaxa <- function(taxalist, from, to, taxanames) {
  from <- colnames(taxanames)[tolower(colnames(taxanames)) %in% tolower(from)]
  if(!length(from)) stop(from,' is not a taxa rank in the dataset')

  to <- colnames(taxanames)[tolower(colnames(taxanames)) %in% tolower(to)]
  if(!length(to)) stop(to,' is not a taxa rank in the dataset')

  taxalist <- unlist(lapply(taxalist, function(x) {
    if(startsWith(x,'X') & is.numeric(type.convert(substr(x,2,2)))) {
      x <- sub('X','',x)
    } else x
  }))

  notfound <- taxalist[!(taxalist %in% taxanames[[from]])]
  if(length(notfound)) message('\nThe following taxa were not found at the ',from,' rank:\n',
                               paste0(notfound,collapse = '\t'))

  taxalist <- taxalist[taxalist %in% taxanames[[from]]]
  if(!length(taxalist)) stop('None of the taxa could be found in the ',from,' rank')

  taxanames <- taxanames[!duplicated(taxanames[[from]]),]
  parents <- taxanames[[to]][taxanames[[from]] %in% taxalist]

  return(parents)
}
