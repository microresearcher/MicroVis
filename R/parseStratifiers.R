#' Parse Stratifiers
#'
#' @param primary_factor The name of the primary factor being analyzed
#' @param factors List of all factors (with their attributes) of the dataset
#' @param stratify Is stratification being performed? If not, this function does
#'     nothing
#' @param facet.x (Optional) Name of a factor for horizontal faceting
#' @param facet.y (Optional) Name of a factor for vertical faceting
#'
#' @return A nested list of the facets and facet text to be used by MicroVis
#'     plotting functions for determining faceting and file autonaming
#'
parseStratifiers <- function(primary_factor, factors, stratify, facet.x, facet.y) {

  facets <- list()

  # "avail_factors" is a list of available factors (any factors that are not
  #   the primary factor) to potentially facet by. This list gets modified
  #   any time a facet is successfully assigned so that its length can be
  #   checked before trying to assign a second facet
  avail_factors <- unlist(names(factors)[names(factors) != primary_factor])
  if(!length(avail_factors)) stratify <- F

  if(is.null(facet.x) & is.null(facet.y)) {
    if(stratify) {
      # If user specified wanting to stratify, but didn't specify either the x or y direction
      facet.x <- TRUE
      facet.y <- TRUE
    } else {
      # User did not call for stratification
      return(facets)
    }
  } else {
    if(!(is.null(facet.x))) {
      if(tolower(facet.x) %in% tolower(avail_factors)) {
        # The factor user specified was found
        facets$x <- factors[tolower(factors) %in% tolower(facet.x)]$name
        avail_factors <- avail_factors[avail_factors != facets$x]

        if(!length(avail_factors)) {
          facets$txt <- nameStratification(c(facets$x,facets$y))
          return(facets)
        }
      } else {
        if(facet.x != TRUE) {message('\nCould not find a secondary factor named ',facet.x,'. Please select an available factor to facet by when prompted')}
        facet.x <- TRUE
      }
    } else {
      facet.x <- FALSE
    }
    if(!(is.null(facet.y))) {
      if(tolower(facet.y) %in% tolower(avail_factors)) {
        # The factor user specified was found
        facets$y <- factors[tolower(factors) %in% tolower(facet.y)]$name
        avail_factors <- avail_factors[avail_factors != facets$y]

        if(!length(avail_factors)) {
          facets$txt <- nameStratification(c(facets$x,facets$y))
          return(facets)
        }
      } else {
        if(facet.y != TRUE) {message('\nCould not find a secondary factor named ',facet.y,'. Please select an available factor to facet by when prompted')}
        facet.y <- TRUE
      }
    } else {
      facet.y <- FALSE
    }
  }

  # If stratify was true and all the above failed to name any of the facets,
  #   have the user select from the list of available secondary factors
  if(facet.x==TRUE & facet.y==TRUE) {
    # If we are trying to assign factors to both a horizontal and vertical facet

    # Horizontal direction first
    facets$x <- select.list(avail_factors,graphics = TRUE, title = 'Horizontal facet')
    avail_factors <- avail_factors[avail_factors != facets$x]

    if(facets$x=='') {facets$x <- NULL}
    if(!length(avail_factors)) {
      facets$txt <- nameStratification(c(facets$x,facets$y))
      return(facets)
    }

    # Now the vertical direction
    facets$y <- select.list(avail_factors,graphics = TRUE, title = 'Vertical facet')
    avail_factors <- avail_factors[avail_factors != facets$y]

    if(facets$y=='') {facets$y <- NULL}
    if(!length(avail_factors)) {
      facets$txt <- nameStratification(c(facets$x,facets$y))
      return(facets)
    }

  } else if(facet.x==TRUE) {
    facets$x <- select.list(avail_factors,graphics = TRUE, title = 'Horizontal facet (enter 0 for none)')
    avail_factors <- avail_factors[avail_factors != facets$x]
    if(!length(avail_factors)) {
      facets$txt <- nameStratification(c(facets$x,facets$y))
      return(facets)
    }

  } else if(facet.y==TRUE) {
    facets$y <- select.list(avail_factors,graphics = TRUE, title = 'Vertical facet (enter 0 for none)')
    avail_factors <- avail_factors[avail_factors != facets$y]
    if(!length(avail_factors)) {
      facets$txt <- nameStratification(c(facets$x,facets$y))
      return(facets)
    }
  }

  if(!length(facets)) {
    message('\nWARNING: No valid factor was selected to facet by -- Proceeding with unstratified analysis and plotting')
  } else {
    # If there is a valid facet that was selected, create a "txt" element of
    #   "facets" that will be printed out in the name of the figure file and
    #   the "stats.all" sub-list
    facets$txt <- nameStratification(c(facets$x,facets$y))
  }

  return(facets)
}
