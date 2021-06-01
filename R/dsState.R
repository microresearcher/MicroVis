#' Change the active rank of a dataset
#'
#' @param dataset Microvis dataset. Default is the active dataset
#' @param rank (Optional) Specify rank to change to. If none is specified, the
#'     function will ask
#'
#' @return MicroVis dataset with the active rank changed to the selected rank
#' @export
#'
changerank <- function(dataset=NULL,rank=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset', envir=mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  ranks <- getRanks(dataset)
  if(length(ranks)==1) {
    message('\nERROR: Only one rank in this dataset')
    return(dataset)
  }

  if(!length(ranks[ranks %in% rank])) rank <- select.list(ranks,title='Select a rank',graphics = T)

  dataset$data$proc$active_rank <- rank

  assign('active_dataset',dataset,envir = mvEnv)
  return(dataset)
}

#' Pool groups of a factor into fewer, larger groups
#'
#' @param dataset Microvis dataset. Default is the active dataset
#' @param factor Groups of this factor will be grouped
#' @param include_prefix Include the factor name as the prefix for each group?
#'
#' @return Dataset with a new factor of pooled groups
#' @export
#'
poolGroups <- function(dataset=NULL,factor=NULL,include_prefix=T) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset', envir=mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }

  if(is.null(factor)) factor <- setFVar(dataset)

  samples <- rownames(dataset$data$orig)
  md <- dataset$metadata
  factors <- dataset$factors

  # pool_name <- paste0('Pooled_',factor$name)
  # i <- 1
  # while(pool_name %in% colnames(md)) {
  #   pool_name <- paste0(pool_name)
  # }
  for(l in capitalize(letters)) {
    pool_name <- paste0('Pooled_',l)
    if(!(pool_name %in% colnames(md))) break
  }

  md[[pool_name]] <- ''

  md[!(md$sample %in% samples)][[pool_name]] <- NA

  if(include_prefix) prefix <- paste0(factor$name,'_')
  else prefix <- ''

  avail_grps <- factor$groups
  i <- 1
  while(length(avail_grps)>0) {
    if(length(avail_grps)>1) {
      cat(paste('\nChoose one or more groups for pool',i,'(excluded groups are also shown here)'))
      pooled_grps <- select.list(avail_grps,multiple = T,graphics = T,
                                 title = paste('Choose one or more groups for pool',i))
    } else pooled_grps <- avail_grps

    if(!length(pooled_grps)) {
      message('\nNo groups selected. All remaining groups in this factor will be pooled together')
      pooled_grps <- avail_grps
    }
    pool <- paste0(prefix,paste(pooled_grps,collapse = '+'))
    md[[pool_name]][md[[factor$name]] %in% pooled_grps] <- pool
    avail_grps <- avail_grps[!(avail_grps %in% pooled_grps)]

    i <- i+1
  }

  current_order <- unique(md[[pool_name]])
  new_order <- list()

  # User will place groups into this new order list one by one in the order they want
  #   (until one item is remaining, at which point the loop will exit)
  for(i in 1:(length(current_order)-1)) {
    new_order <- c(new_order,select.list(current_order[!(current_order %in% new_order)],
                                         title=paste0('Choose group #',i,':'),graphics=TRUE))
    # Print the new order after each selection so that the user can keep track of their ordering
    cat('\n',paste0('\n',new_order))
  }

  # Once one item is remaining, we don't need to ask the user to 'choose' that item for the last slot
  #   which is why the loop exited when 1 item was remaining
  new_order <- unlist(c(new_order, current_order[!(current_order %in% new_order)]))
  cat('\n\nGroup order for ', pool_name, 'is:\n',paste0(new_order,collapse='\n '))

  md[[pool_name]] <- factor(md[[pool_name]],levels=new_order)

  pool_factor <- list(name=pool_name,
                      name_text=gsub('_',' ',pool_name),
                      groups=levels(md[[pool_name]]))

  pool_factor$subset <- new_order[new_order %in% (md[md[[factor$name]] %in% factor$subset,][[pool_name]])]

  factors[[pool_name]] <- pool_factor

  defCols <- get('defCols',envir = 'mvEnv')

  avail_clrs <- defCols[!(defCols %in% dataset$colors)]
  if(length(avail_clrs)<length(pool_factor$groups)) avail_clrs <- defCols

  newclrs <- avail_clrs[1:length(pool_factor$groups)]
  names(newclrs) <- pool_factor$groups

  dataset$colors <- c(dataset$colors,newclrs)

  dataset$metadata <- md
  dataset$factors <- factors

  dataset <- processDataset(dataset,silent = T)

  return(dataset)
}

#' Remove a factor of pooled groups
#'
#' @param dataset Dataset to remove pooling factor from
#'
#' @return Returns a dataset without the pooling factor that was selected for
#'     for removal
#' @export
#'
removePool <- function(dataset=NULL) {
  if(is.null(dataset)) {
    dataset <- get('active_dataset',envir = mvEnv)
    dataset_name <- 'active_dataset'
  } else {
    dataset_name <- deparse(substitute(dataset))
  }
}
