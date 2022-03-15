#' Reset the Results Directory
#'
#' @param save_directory Current results directory path
#'
#' @return New results directory path.
#' @export
#'
resetResDir <- function(save_directory) {
  if(!(getwd()==save_directory)) {
    cwd <- strsplit(getwd(),split = '/')[[1]]
    savedir <- strsplit(save_directory,split = '/')[[1]]

    cwd_intersect <- grep(intersect(cwd,savedir)[1],cwd)[1]
    if(is.na(cwd_intersect)) return(save_directory)
    if(cwd_intersect==1) return(save_directory)

    savedir_intersect <- grep(intersect(cwd,savedir)[1],savedir)[1]
    new_savedir <- c(cwd[1:(cwd_intersect-1)],savedir[savedir_intersect:length(savedir)])

    save_directory <- new_savedir[1]
    for(d in 2:length(new_savedir)) save_directory <- file.path(save_directory,new_savedir[d])
  }
  return(save_directory)
}

#' Save figure and statistics results
#'
#' @param dataset Dataset to associate these results with. Directory path is pulled
#'     from the dataset
#' @param foldername (Optional) Name of the directory within save_directory
#' @param filename (Optional) Name of file to save figure and statistics under
#' @param factors (Optional) Factors of the dataset
#' @param analysis_name (Optional) Name for this analysis
#' @param active_factor (Optional) The active factor of this analysis
#' @param facets (Optional) Facets if there were any facets in this analysis
#' @param figure (Optional) The plot object to be saved
#' @param width (Optional) Width of the image to be saved in inches
#' @param height (Optional) Height of the image to be saved in inches
#' @param stat_results (Optional) Table of statistics results to be saved
#' @param other_results (Optional) Table of other results to be saved
#' @param suffix (Optional) Any suffix to be added to the end of the file name
#' @param forcesave If set to TRUE, the function acts as if "autosave" is TRUE
#' @param verbose If set to TRUE, will print out the location of file saving
#'
#' @return save_directory
#' @export
#'
saveResults <- function(dataset,
                        foldername='',
                        filename=NULL,
                        factors=NULL,
                        active_factor=NULL,
                        analysis_name=NULL,
                        facets=NULL,
                        figure=NULL,
                        width=18,height=12,
                        stat_results=NULL,
                        other_results=NULL,
                        suffix='',
                        forcesave=F,
                        verbose=T) {
  autosave <- get('autosave',envir = mvEnv)
  offerSave <- get('offerSave',envir = mvEnv)
  autoNameResults <- get('autoNameResults',envir = mvEnv)

  if(autosave | forcesave) {
    saveFig <- T
  } else if(offerSave) {
    saveFig <- ifelse(select.list(c('Yes','No'),
                                  title='\nWould you like to save this figure?')=='Yes',
                      TRUE,FALSE)
  } else return()

  if(saveFig) {
    # Create the output directories
    save_directory <- file.path(resetResDir(dataset$results_path),
                                paste0('Results_',Sys.Date()),
                                foldername)
    dir.create(save_directory,recursive = T,showWarnings = F)

    subfoldername <- ''

    if(autoNameResults=='groups') {
      if(!is.null(analysis_name)) filename <- analysis_name
      # No analysis_name, but factors and active_factor is specified
      else if(!is.null(factors) & !is.null(active_factor)) {
        # No filename is specified, in which case make one and don't make a subfolder
        if(is.null(filename)) {
          filename <- nameAnalysis(factors, active_factor, facets)
          # Check if path is too long
          if(nchar(file.path(save_directory,paste0(filename,suffix,'.png'))) > 259) {
            filename <- NULL
            autoNameResults <- 'dataset'
          }
        }
        # A filename is specified, in which case make a subfolder to put that filename file in
        else {
          subfoldername <- nameAnalysis(factors, active_factor, facets)

          # Check if path is too long
          if(nchar(file.path(file.path(save_directory, subfoldername),
                             paste0(filename,suffix,'.png'))) > 259) autoNameResults <- 'dataset'
        }
      }
    }

    if(autoNameResults=='dataset') {
      # Check if dataset has a name. If not, ask user for one
      if(is.null(dataset$name)) {
        charlimit <- 259 - nchar(file.path(file.path(save_directory,subfoldername),paste0(suffix,'.png')))

        cat('\n')
        dataset_name <- readline(paste0('Please type a name for this analysis (up to ',charlimit,' characters): '))

        # Ask user if they would like to use this name for the dataset
        if(select.list(title = "\nWould you like to save this as the name for the dataset?",
                       choices = c('Yes','No'))=='Yes') mvsave(dataset_name)

      } else dataset_name <- dataset$name

      # No filename is specified, in which case make one and don't make a subfolder
      if(is.null(filename)) {
        filename <- paste0(dataset$name,'_',
                           strsplit(nameAnalysis(factors, active_factor, facets),'_')[[1]][1])
        if(!is.null(facets$txt)) filename <- paste0(filename, facets$txt)
        # Check if path is too long
        if(nchar(file.path(save_directory,paste0(filename,suffix,'.png'))) > 259) {
          filename <- NULL
          autoNameResults <- 'manual'
        }
      }
      # A filename is specified, in which case make a subfolder to put that filename file in
      else {
        subfoldername <- paste0(dataset$name,'_',
                                strsplit(nameAnalysis(factors, active_factor, facets),'_')[[1]][1])
        if(!is.null(facets$txt)) subfoldername <- paste0(subfoldername, facets$txt)
        # Check if path is too long
        if(nchar(file.path(file.path(save_directory, subfoldername),
                           paste0(filename,suffix,'.png'))) > 259) autoNameResults <- 'manual'
      }
    }

    if(autoNameResults=='manual') {
      # No filename is specified, in which case make one and don't make a subfolder
      if(is.null(filename)) {
        charlimit <- 259 - nchar(file.path(save_directory,paste0(suffix,'.png')))
        # Ask user for a custom name for the figure
        message('\nPlease type a file name for this figure (up to ',charlimit,' characters)\n')
        filename <- readline('File name: ')
      }
      # A filename is specified, in which case make a subfolder to put that filename file in
      else {
        charlimit <- 258 - nchar(file.path(save_directory,paste0(filename,suffix,'.png')))
        message('\nPlease type a folder name up to ',charlimit,' characters for these figures\n')
        subfoldername <- readline('Folder name: ')
        # save_directory <- file.path(save_directory, subfoldername)
        # dir.create(save_directory,recursive = T,showWarnings = F)
      }
    }

    save_path <- file.path(file.path(save_directory,subfoldername),paste0(filename,suffix,'.png'))

    if(is.null(figure)) ggsave(save_path, device='png',
                               width=width, height=height, units = 'in',
                               dpi=600)
    else {
      png(save_path,
          width=width,
          height=height,
          units='in',
          res=600)
      draw(figure)
      dev.off()
    }

    if(!is.null(stat_results)) {
      dir.create(file.path(save_directory,'Statistics'),showWarnings = F)
      if(!is.null(stat_results$stats)) {
        stats <- apply(stat_results$stats,2,function(x) as.character(x))
        write.csv(stats,
                  file=file.path(save_directory,'Statistics',
                                 paste0(filename,suffix,'.csv')),
                  row.names = F)
      }
      if(!is.null(stat_results$pw_stattab)) {
        pw_stats <- apply(stat_results$pw_stats,2,function(x) as.character(x))
        write.csv(pw_stats,
                  file=file.path(save_directory,'Statistics',
                                 paste0(filename,suffix,'-pairwise.csv')),
                  row.names = F)
      }
    }

    if(!is.null(names(other_results))) {
      for(res in names(other_results)) {
        dir.create(file.path(save_directory,res),showWarnings = F)
        write.csv(other_results[[res]],
                  file=file.path(save_directory,res,
                                 paste0(filename,suffix,'.csv')),
                  row.names = F)
      }
    }

    if(verbose) return(cat('Figure(s) and any associated statistics saved to:\n ',save_path,'\n'))
    return(save_directory)
  }
}

#' Determine how to save multiple figures
#'
#' @param save_one_all Variable to determine how to save figures
#'
#' @return save_one_all
#' @export
#'
multisave <- function(save_one_all=NULL) {
  if(is.null(save_one_all)) {
    if(get('autosave',envir = mvEnv)) save_one_all <- 'Yes to all figures'
    else if(get('offerSave',envir = mvEnv)) {
      save_one_all <- select.list(c('Yes',
                                    'No',
                                    'Yes to all figures',
                                    'No to all figures'),
                                  title = 'Save this figure?',
                                  graphics = T)
    }
    else save_one_all <- 'No to all figures'

  } else if(save_one_all %in% c('Yes','No')) {
    save_one_all <- select.list(c('Yes',
                                  'No',
                                  'Yes to all figures',
                                  'No to all figures'),
                                title = 'Save this figure?',
                                graphics = T)
  }
  return(save_one_all)
}

#' Automatically name the analysis based on factors and the active factor
#'
#' @param factors List of factors of the dataset
#' @param active_factor The name of the active factor of the dataset
#' @param facets (Optional) Any facets used in an analysis
#'
#' @return Name of the analysis as a string
#' @export
#'
nameAnalysis <- function(factors, active_factor, facets=NULL) {
  othertxt <- ''

  for(f in factors) {
    excluded_grps <- f$groups[!(f$groups %in% f$subset)]
    if(!length(excluded_grps)) {
      # If no groups are excluded
      if(f$name==active_factor) {
        # All groups is implied for non-active factors unless specified using
        #   the lines of code below
        temp <- paste0('All ',f$name,'s')
      } else next
    } else if(length(excluded_grps) < length(f$subset)) {
      # If less than half of the groups are excluded
      temp <- paste0('No ',paste(excluded_grps, collapse = '-'))
      # temp <- paste0('No ',paste(excluded_grps, collapse = '-'),' ',f$name,'s')
    } else {
      # If less than half of the groups are included
      temp <- paste0('Only ',paste(f$subset, collapse = '-'))
      # temp <- paste0('Only ',paste(f$subset, collapse = '-'),' ',f$name,'s')
    }
    if(f$name==active_factor) maintxt <- temp
    else othertxt <- paste0(othertxt,paste0('_',temp))
  }

  analysistxt <- paste0(maintxt,othertxt)
  if(!is.null(facets$txt)) analysistxt <- paste0(analysistxt, '_', facets$txt)
  return(analysistxt)
}
