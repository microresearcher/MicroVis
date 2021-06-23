#### Template Function ####
# This function is for other developers or those who fork the MicroVis repo.
#    Feel free to use this as a basis for any functions you make!
templateFunction <- function(dataset=NULL) {
  if(is.null(dataset)) dataset <- get('active_dataset',envir = mvEnv)

  ### Your code here

  # If this is a plotting function, it is recommended
  #   to activate the dataset (uncomment the line below)
  #activate(dataset)

  # If this is an analysis function that takes a while, it is recommended to
  #   store the results in the dataset and rewrite the dataset in the global
  #   and MicroVis environments (uncomment the lines below)
  #assign('active_dataset',dataset,envir = mvEnv)
  #if(!is.null(dataset$name)) assign(dataset$name,dataset,1)
}
