#' Choose Factors for the Dataset
#'
#' @param dataset MicroVis dataset (mvdata object)
#'
#' @return Microvis dataset with factors data
#' @export
#'
chooseFactors <- function(dataset) {
  if(is.null(dataset$metadata)) return(dataset)
  else metadata <- dataset$metadata

  # Initialize it here in case metadata doesn't have more than 1 column and the
  #   first if-statement block isn't run
  possible_factors <- c(colnames(metadata)[2:(ncol(metadata))])
  chosen_factors <- list()
  if(length(possible_factors) > 1) {
    cat('\n')
    print(head(metadata))
    chosen_factors <- select.list(possible_factors,
                                  multiple=TRUE,
                                  graphics=TRUE,
                                  title='\nSelect all the potential factors you would like to compare data by:')
    cat('\n')
  } else chosen_factors <- possible_factors

  factors <- list()
  # Make sure at least one factor was chosen
  if(length(chosen_factors)) {
    # Record name, groups, and 'cleaned' name of the factor
    #   Also create a 'subset' element that will contain only selected groups
    #   that the user chooses
    for(f in chosen_factors) {
      # if(is.numeric(type.convert(metadata[[f]], as.is=T))) {
      if(is.numeric(metadata[[f]])) {
        cutoffs <- -Inf
        minval <- min(metadata[[f]],na.rm = T)
        maxval <- max(metadata[[f]],na.rm = T)
        prompt <- paste0('For "',f,'", type a space-separated list of numbers between ',minval,' and ',maxval,' for cutoffs: ')
        while(!all(dplyr::between(cutoffs, minval, maxval))) {
          cutoffs <- as.numeric(unlist(strsplit(readline(prompt),' ')))
        }
        prompt <- 'Are these cutoffs the lower or upper limit?'
        lowerlim <- ifelse(select.list(c('Lower limit',
                                         'Upper limit'),
                                       title = prompt)=='Lower limit',yes = T,no = F)
        if(lowerlim) {
          groups <- rangetotext(cut(metadata[[f]],c(cutoffs,minval,Inf),right=F))
        } else {
          groups <- rangetotext(cut(metadata[[f]],c(cutoffs,-Inf,maxval),right=T))
        }
        f <- chosen_factors[grepl(f,chosen_factors)] <- paste0(f,'_Range')
        metadata[[f]] <- groups
        dataset$metadata <- metadata
      }
      factors[[f]] <- list()
      factors[[f]]$name <- f
      factors[[f]]$name_text <- stringr::str_replace(f,'_',' ')
      factors[[f]]$groups <- levels(as.factor(metadata[[f]]))
      factors[[f]]$subset <- factors[[f]]$groups
    }

    if(length(chosen_factors)>1) {
      # Have user choose the primary factor if more than 1 potential factor was
      #   selected
      primary_factor <- select.list(names(factors),title='\nWhat is the primary factor you would like to compare by?',graphics=TRUE)
    } else primary_factor <- chosen_factors
  } else {
    # If no factors were chosen, warn the user but keep going
    message('\nWARNING: No factors selected. Comparative analysis will not be available')
    assign('warning_list',c(get('warning_list',envir = mvEnv),paste0('WARNING: No factors selected. Comparative analysis will not be available')),envir = mvEnv)
  }

  dataset$factors <- factors
  dataset$active_factor <- primary_factor

  if(!get('.loading',envir = mvEnv)) {
    dataset <- orderGroups(dataset)
  }

  return(dataset)
}

#' Translate Range Factor Levels into Interpretable Text
#'
#' @param factor Factor (list object) with its groups to be translated
#'
#' @return New factor (list object) with name appended with "_Range" and group
#'     names changed to the the formats as appropriate: "# or less", "# to #",
#'     or "# or more"
#'
rangetotext <- function(factor) {
  renamedlevels <- levels(factor)
  renamedgroups <- as.character(factor)

  for(range in renamedlevels) {
    rangesplit <- unlist(strsplit(range,','))
    if(any(grepl('Inf',rangesplit))) {
      if(grepl('-Inf',rangesplit[1])) {
        lowertext <- paste0(gsub('\\]','',gsub('\\)','',rangesplit[2])))
        uppertext <- ' or less'
      } else if(grepl('Inf',rangesplit[2])) {
        lowertext <- paste0(gsub('\\[','',gsub('\\(','',rangesplit[1])))
        uppertext <- ' or more'
      }
    } else {
      lowertext <- paste0(gsub('\\[','',gsub('\\(','',rangesplit[1])))
      uppertext <- paste0(gsub('\\]','',gsub('\\)','',rangesplit[2])))
    }

    if(is.numeric(type.convert(c(lowertext,uppertext), as.is=T))) rangetext <- paste(lowertext,'to',uppertext)
    else rangetext <- paste0(lowertext,uppertext)

    renamedlevels[renamedlevels==range] <- rangetext
    renamedgroups[renamedgroups==range] <- rangetext
  }
  renamedfactor <- factor(renamedgroups,levels=renamedlevels)
  return(renamedfactor)
}
