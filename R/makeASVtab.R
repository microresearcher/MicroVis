#' @title Make ASV Table
#'
#' @description Generates an ASV table based on the taxa present in the taxonomy
#'     csv file. Will attempt to guess how many ranks are represented in the
#'     taxa names if they are not explicitly stated in the file.
#'
#' @param taxa_names_tab List of unprocessed taxa that may or may not already
#'     be split into columns for each rank
#'
#' @return Dataframe with taxa in rows, split into columns for each taxonomic rank.
#'
makeASVtab <- function(taxa_names_tab) {
  taxaranks <- c('domain','phylum','class','order','family','genus','species')
  # List of invalid characters in taxonomy names to replace with an underscore
  invalid_chars <- c('-'='_',' '='_',':'='_','/'='_',
                     '\\['='','\\]'='','\\('='_','\\)'='')
  taxa_names_tab <- sapply(taxa_names_tab, function(x) str_replace_all(x,invalid_chars))
  # for(c in invalid_chars) taxa_names_tab <- sapply(taxa_names_tab, function(x) str_replace_all(x,c,'_'))

  names(taxa_names_tab) <- taxa_names_tab
  if(!identical(names(taxa_names_tab),unname(taxa_names_tab))) {
    message('\nWARNING: Taxa names have been changed after turning them into headers')
    warning_list <- c(warning_list, '\nWARNING: Taxa names have been changed after turning them into headers')
  }

  # Check if the taxa names are prepended lineages or single rank, and identify
  #   the delimiter if they are prepended names
  if(any(grepl(';',taxa_names_tab))) {
    prepended_taxa_names <- TRUE
    delim <- ';'
    cat('\n\n Prepended taxonomy names detected. Abundance tables for each taxonomic rank will be created\n')
  } else if(any(grepl(' | ',taxa_names_tab))) {
    prepended_taxa_names <- TRUE
    delim <- ' | '
    cat('\n\n Prepended taxonomy names detected. Abundance tables for each taxonomic rank will be created\n')
  } else prepended_taxa_names <- FALSE

  if(prepended_taxa_names) {
    # Store original taxa names in a vector
    prepended_names <- as.vector(taxa_names_tab)

    # Make sure each prepended name has the same number of delimiters
    if(length(unique(str_count(prepended_names,';')))>1) {
      stop('Taxonomy names have different number of rank levels: ',
                                                              paste(unique(str_count(prepended_names,
                                                                                     ';')),
                                                                    collapse = ', '))
    }

    # Select the names of all the ranks based on how many sections each
    #   prepended name has
    #   (e.g. 4 delimiters -> 5 sections -> kingdom to family)
    ranks <- taxaranks[1:(max(str_count(prepended_names,delim))+1)]

    # Clean up these names. Silva and GreenGenes refer to the top taxa rank as
    #   "kingdom", so switch to that temporarily
    ranks[1] <- 'kingdom'
    for(rank in ranks) {
      # Silva adds tags to the end of taxonomy names when a genus has the same
      #   name at multiple taxonomic ranks (e.g. "_o" for order)
      prepended_names <- sapply(prepended_names,
                                function(x) gsub(paste0('_',substr(rank,1,2)),'',x))
      # GreenGenes adds a prefix in front of each taxonomy name to denote the
      #   taxonomic rank (e.g. "p__" for phylum)
      prepended_names <- sapply(prepended_names,
                                function(x) gsub(paste0(substr(rank,1,1),'__'),'',x))
    }
    # Switch back to "domain" as the top taxa rank name
    ranks[1] <- 'domain'

    # If there is no text remaining after the last delimiter because the taxa
    #   is unidentified at the lowest rank, put "NA" there so that it doesn't
    #   get chopped off when splitting by the delimiter
    for(i in 1:length(prepended_names)) if(endsWith(prepended_names[[i]],delim)) {
      prepended_names[[i]] <- paste0(prepended_names[[i]],'NA')
    }

    taxa_names_tab <- data.frame(t(sapply(prepended_names,
                                          function(x) strsplit(x,delim)[[1]])))
    colnames(taxa_names_tab) <- ranks

    taxa_names_tab <- cleanASVs(taxa_names_tab)
  } else {
    taxa_names_tab <- data.frame(taxa_names_tab)
    colnames(taxa_names_tab) <- 'single_rank'
  }

  return(taxa_names_tab)
}

#' @title Clean ASVs in ASV Dataframe
#'
#' @description Rename NA taxa to refer to their highest resolution of identification
#'     and make sure species names are prepended with the genus name
#'
#' @param taxa_names_tab Processed ASV dataframe
#'
#' @return Returns an ASV dataframe with taxa names that make more sense and are
#'     easier to deal with.
#'
cleanASVs <- function(taxa_names_tab) {
  ranks <- colnames(taxa_names_tab)

  taxa_names_tab <- apply(taxa_names_tab,1,function(asv) {
    for(rank in 1:length(asv)) {
      if(is.na(asv[rank])) {
        if(rank==1) {
          asv[rank] <- paste0('unidentified_domain')
          highest_res <- asv[rank]
        } else asv[rank] <- paste0(highest_res,'_',ranks[rank])
      } else if(asv[rank]=='NA' | asv[rank]=='' | asv[rank]=='Incertae_Sedis') {
        if(rank==1) {
          asv[rank] <- paste0('unidentified_domain')
          highest_res <- asv[rank]
        } else asv[rank] <- paste0(highest_res,'_',ranks[rank])
      } else {
        highest_res <- paste0(substr(ranks[rank],1,1),'_',asv[rank])
        if(ranks[rank]=='species') if(!(grepl(asv[rank-1],asv[rank]))) {
          asv[rank] <- paste0(asv[rank-1],'_',asv[rank])
        }
      }
    }
    return(asv)
  })
  taxa_names_tab <- data.frame(t(taxa_names_tab))
}

#' Combine Identical ASVs
#'
#' @param taxa_data_list Taxa data with taxa names and un-processed abundance table.
#' @param combineDupes If set to TRUE (default), MicroVis will try to combine
#'     duplicate features.
#'
#' @return Taxa data with identical taxa combined
combineDupeASVs <- function(taxa_data_list,combineDupes=T) {
  taxa_data <- taxa_data_list$orig
  taxa_names_tab <- taxa_data_list$taxa_names

  dupes <- data.frame(unique(taxa_names_tab[duplicated(taxa_names_tab),]))

  # Combine all exact duplicates
  if(nrow(dupes)) {
    temp <- data.frame(taxa_names_tab,ASV=rownames(taxa_names_tab))
    if(combineDupes) {
      for(dupe in 1:nrow(dupes)) {
        dupeasvs <- merge(dupes[dupe,],temp)$ASV
        taxa_data[dupeasvs[1]] <- rowSums(taxa_data[dupeasvs])
        taxa_data[dupeasvs[2:length(dupeasvs)]] <- NULL
        taxa_names_tab <- taxa_names_tab[!(rownames(taxa_names_tab) %in% dupeasvs[2:length(dupeasvs)]),]
      }
    }
  }

  if(ncol(taxa_names_tab)>1) {
    # Check for any duplicates within each rank
    ranks <- colnames(taxa_names_tab)
    for(rank in ranks[2:length(ranks)]) {
      dupes <- taxa_names_tab[rank][duplicated(taxa_names_tab[rank]),]
      for(dupe in unique(dupes)) {
        parents <- apply(data.frame(taxa_names_tab[taxa_names_tab[[rank]]==dupe,
                                                   seq(grep(rank,ranks)-1)]),
                         1,paste0,collapse=';')
        if(length(unique(parents))==1) dupes <- dupes[dupes!=dupe]
      }

      if(length(dupes)) {
        # Identify the next lowest rank that uniquely identifies these taxa
        prepend_ranks <- sapply(unique(dupes),
                                function(dupe) {
                                  tab <- unique(taxa_names_tab[taxa_names_tab[[rank]] %in% dupe,
                                                               seq(grep(rank,ranks))])
                                  max(sapply(rev(seq(grep(rank,ranks)-1)), function(higher_rank) {
                                    ifelse(length(tab[[higher_rank]])==length(unique(tab[[higher_rank]])),
                                           higher_rank, 0)
                                  }))
                                })

        taxa_names_tab <- apply(taxa_names_tab, 1, function(asv) {
          if(asv[[rank]] %in% dupes) {
            prepend_rank <- prepend_ranks[[asv[[rank]]]]

            if(!grepl('^[dkpcofg]_',asv[[prepend_rank]])) {
              prepend <- paste0(substr(ranks[prepend_rank],1,1),'_',asv[[prepend_rank]])
            } else prepend <- asv[[prepend_rank]]

            asv[rank] <- paste0(prepend,'_',gsub('[dkpcofg]_.*_','',asv[[rank]]))
          }
          return(asv)
        })
        taxa_names_tab <- data.frame(t(taxa_names_tab))
      }
    }
  }

  taxa_data_list$orig <- taxa_data
  taxa_data_list$taxa_names <- taxa_names_tab

  return(taxa_data_list)
}
