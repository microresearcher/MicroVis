#' Make Abundance Tables at Each Rank
#'
#' @param ft_data Feature data containing table of feature names - if taxa,
#'     then separated by rank - an original abundance table, and an unranked
#'     abundance table. If features are not taxa, then this function returns
#'     the input unchanged.
#'
#' @return Feature data with table of feature names, an original abundance table
#'     and processed abundance tables at each taxonomy rank.
#'
makeRankTabs <- function(ft_data) {
  new_ft_data.proc <- ft_data$proc
  abd <- new_ft_data.proc$unranked
  sample_names <- rownames(abd)

  if(!is.null(ft_data$taxa_names)) {
    taxaRanks <- get('taxaRanks',envir = mvEnv)
    taxa_names_tab <- ft_data$taxa_names
    if(ncol(taxa_names_tab)>1) {
      # if(colnames(taxa_names_tab) %in% taxaRanks) taxaRanks <- colnames(taxa_names_tab)
      # Make abundance tables for all taxonomy ranks
      for(i in c(1:ncol(taxa_names_tab))) {
        abd_temp <- abd

        # Go through each ASV in the abundance table and assign the taxon name
        #   at the appropriate taxa level. This way, if we are using the
        #   filtered abundance table, we only look for the taxa names of
        #   remaining taxa
        asvs <- colnames(abd_temp)[!(colnames(abd_temp) %in% c('Other','Unknown'))]
        for(asv in asvs) {
          colnames(abd_temp)[which(colnames(abd_temp)==asv)] <- taxa_names_tab[asv,i]
        }
        unique_taxa <- unique(colnames(abd_temp))
        combined <- list()
        for(taxon in unique_taxa) {
          combined[[taxon]] <- rowSums(data.frame(abd_temp[,colnames(abd_temp)==taxon]))
        }
        new_ft_data.proc[[colnames(taxa_names_tab)[i]]] <- data.frame(combined, row.names = sample_names)
        # Un-rename any taxa that were silently renamed by R (eg R adds an 'X' in front
        #   of column names that start with a number, so we want to reverse that)
        colnames(new_ft_data.proc[[colnames(taxa_names_tab)[i]]]) <- sapply(colnames(new_ft_data.proc[[colnames(taxa_names_tab)[i]]]), function(x) {
          if(startsWith(x,'X') & is.numeric(type.convert(substr(x,2,2), as.is=T))) {
            x <- sub('X','',x)
          } else x
        })
      }
      if(is.null(new_ft_data.proc$active_rank)) {
        new_ft_data.proc$active_rank <- rev(taxaRanks[taxaRanks %in% colnames(taxa_names_tab)])[[1]]
      }
    } else {
      asvs <- colnames(abd)[!(colnames(abd) %in% c('Other','Unknown'))]
      for(asv in asvs) {
        colnames(abd)[which(colnames(abd)==asv)] <- taxa_names_tab[rownames(taxa_names_tab)==asv,1]
      }
      new_ft_data.proc[['single_rank']] <- data.frame(abd, row.names = sample_names)
      new_ft_data.proc$active_rank <- 'single_rank'
    }
  } else {
    new_ft_data.proc[['functional']] <- data.frame(abd, row.names = sample_names)
    new_ft_data.proc$active_rank <- 'functional'
  }

  ft_data$proc <- new_ft_data.proc

  return(ft_data)
}
