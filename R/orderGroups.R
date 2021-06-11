#' Order Groups
#'
#' @param dataset MicroVis dataset (mvdata object)
#' @param reorder Set this to true if you want to reorder groups after they have
#'     already been ordered in initial loading (during mvload)
#'
#' @return MicroVis dataset with groups in each factor ordered according to user's
#'     input
#' @export
#'
orderGroups <- function(dataset, reorder=F) {
  grps_saved.unlisted <- list()
  clrd_grps <- list()

  metadata <- dataset$metadata
  factors <- dataset$factors

  # If a save-file exists, any groups in this save-file and their order will be loaded
  if(file.exists(file.path(get('project_dir',envir = mvEnv),'group_data.txt'))) {
    save_file <- file.path(get('project_dir',envir = mvEnv),list.files(get('project_dir',envir = mvEnv),
                                                                       pattern=paste0('group_data.txt')))
    save_data <- strsplit(readLines(save_file),split='\t')
    #close(save_file)

    grps_saved.unlisted <- save_data[[grep('Groups:',save_data)+1]]
    grp_clrs <- save_data[[grep('Colors:',save_data)+1]]
    #grps_saved.unlisted <- save_data[[1]]
    #grp_clrs <- save_data[[2]]

    # For the groups that have assigned colors, create a temporary named list of these groups and their colors
    for(i in 1:length(grps_saved.unlisted)) {clrd_grps[grps_saved.unlisted[i]] <- grp_clrs[i]}
  }

  # For each factor, if there were any new groups not in the save-file, list the groups and ask user for reordering
  for(f in factors) {
    # Don't do this for the $active_factor attribute of the factors list
    if(length(f)==1) {next}

    # For this factor, make a list of the groups and their order that was loaded from the save-file
    grps_saved <- grps_saved.unlisted[grps_saved.unlisted %in% f$groups]
    # Any groups in this factor that are not in this ordered list are new and need to be assigned colors
    new_grps <- f$groups[!(f$groups %in% grps_saved)]

    # The initial assumed order will be [order of groups for that factor as listed in save-file, then any new groups]
    #   This list should contain ALL the groups in flist$levels[[i]] (the levels for this factor),
    #   with any groups from the save-file ordered the way they are in the save-file
    current_order <- unlist(c(grps_saved, new_grps))

    # If there are any new groups that weren't in the save-file,
    #   or the user wanted to reorder the groups when they ran this function,
    #   ask if they would like to reorder the groups in *this* factor
    # This way, if they wanted to reorder the groups in one factor,
    #   they don't have to order the groups in *all* the factors
    if((!length(new_grps) | length(new_grps)>10) & !reorder) {
      change_order <- FALSE
    } else {
      cat(paste0('\n\nCurrent group order for ',f$name_text,':'))
      cat(paste0('\t',current_order))
      change_order <- ifelse(select.list(choices = c('Keep','Change'),
                                         title = '\n\nChange the ordering of the groups?')=='Change',TRUE,FALSE)
    }

    # If the order of groups for this factor doesn't need to be changed, the current order is the "new" order
    #   If it does need to be changed, then ask the user to order the groups how they would like them, and that will be the new order
    if(!change_order) {
      new_order <- current_order
      cat('\n\nGroup order for', f$name_text, 'is:\n',paste0(new_order,collapse='\t'))
    } else {
      new_order <- list()

      # User will place groups into this new order list one by one in the order they want
      #   (until one item is remaining, at which point the loop will exit)
      for(j in 1:(length(f$groups)-1)) {
        new_order <- c(new_order,select.list(current_order[!(current_order %in% new_order)],title=paste0('Choose group #',j,':'),graphics=TRUE))
        # Print the new order after each selection so that the user can keep track of their ordering
        cat('\n',paste0('\n',new_order))
      }

      # Once one item is remaining, we don't need to ask the user to 'choose' that item for the last slot
      #   which is why the loop exited when 1 item was remaining
      new_order <- unlist(c(new_order, current_order[!(current_order %in% new_order)]))
      cat('\n\nGroup order for ', f$text, 'is:\n',paste0(new_order,collapse='\t'))
    }

    # Reorder the "levels" list of this factor
    f$groups <- new_order
    # Reorder the "grps" list of this factor
    f$subset <- new_order[new_order %in% f$subset]
    factors[[f$name]] <- f

    # Regardless of new groups added or whether the order is changed, this is where the metadata column is factorized
    #   in the appropriate order, now that the group order is finalized
    metadata[[f$name]] <- factor(metadata[[f$name]], levels=f$groups)
  }

  dataset$metadata <- metadata
  dataset$factors <- factors

  # Finally, color any new groups while maintaining the colors of any old ones
  new_order.unlisted <- c()
  for(f in factors) {
    if(length(f)>1) {
      new_order.unlisted <- c(new_order.unlisted, f$groups)
    }
  }
  dataset <- colorGrps(dataset, new_order.unlisted, clrd_grps)

  # If this is not run during loading, then "cmpgrp" will not be reset without
  #   running "chooseGrps()" which is undesirable
  return(dataset)
}
