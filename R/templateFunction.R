#### Template Function ####
# This function is for other developers or those who fork the MicroVis repo.
#    Feel free to use this as a basis for any functions you make!
templateFunction <- function(dataset=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  ### Your code here

  # If this is a plotting function or analysis function, it is recommended
  #   to activate the dataset (uncomment the line below)
  #activate(dataset)
}
