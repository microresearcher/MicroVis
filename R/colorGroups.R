#' Color Groups
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param new_order.unlisted Single list of any new groups across all factors
#'     with their ordering retained
#' @param clrd_grps Named list of already colored groups and their colors
#'
#' @return MicroVis dataset (mvdata object) with a colors attribute containing
#'     the color assignment for each group
#'
colorGrps <- function(dataset, new_order.unlisted=NULL, clrd_grps=NULL) {
  # If no factors are selected for this analysis, skip this function and assign
  #   default colors to "clrs" variable (for now, maybe we don't need this)
  if(!length(dataset$factors)) {
    dataset$colors <- get('defCols',envir = mvDefaults)
    return(dataset)
  }

  defCols <- get('defCols',envir = mvEnv)

  # Initialize an empty named colors list
  #   This will be the global list with groups as names and their corresponding
  #   colors as their elements in the list
  clrs <- list()

  # Create a list of unused colors
  new_clrs <- defCols[!(defCols %in% unlist(clrd_grps))]

  # For each group in the new list of ordered groups
  for(grp in new_order.unlisted) {
    # If this group was not colored already (i.e. it was not in the save-file)
    #   then give it a color from unused colors
    if(!(grp %in% names(clrd_grps))) {
      clrd_grps[grp] <- new_clrs[1]
      new_clrs <- new_clrs[-1]
      # If the unused colors are all used up, refresh it with the original list
      #   of all colors
      if(!length(new_clrs)) {new_clrs <- defCols}
    }
    clrs[grp] <- clrd_grps[grp]
  }
  # If there were any groups in the save-file that are not in any of the groups
  #   in this analysis, we still want to retain their colors. So we will tack
  #   them onto the end of the "clrs" list and print them out in the save-file
  grps_left <- names(clrd_grps)[!(names(clrd_grps) %in% names(clrs))]
  for(grp in grps_left) {clrs[grp] <- clrd_grps[grp]}

  save_file <- file.path(get('project_dir',envir = mvEnv),'group_data.txt')
  writeLines(c('Groups:',paste(names(clrs),collapse='\t'),'Colors:',paste(clrs,collapse='\t')),save_file)

  dataset$colors <- unlist(clrs)
  return(dataset)
}
