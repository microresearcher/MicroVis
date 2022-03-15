#' @title MicroVis Data Loading
#'
#' @description Initial function for loading taxonomy and/or pathway data with
#'     associated metadata into MicroVis
#'
#' @param path_to_folder (Optional) full path to directory containing csv files
#' @param autoProcess If set to TRUE (default), MicroVis will automatically
#'     process datasets using default normalization and filtering parameters.
#'     Datasets can still be re-processed by the user afterwards.
#' @param combineDupes If set to TRUE (default), MicroVis will try to combine
#'     duplicate features.
#' @param combineDataSets If set to TRUE, MicroVis will ask for a second set of
#'     data and combine both datasets into one.
#' @param path_to_metadata (Optional) Path to metadata csv file
#' @param path_to_taxa (Optional) Path to taxonomy abundance csv file
#' @param path_to_fxnl (Optional) Path to functional abundance csv file
#'
#' @return List of successfully loaded datasets with/without processing.
#' @export
#'
#' @importFrom graphics barplot
#' @importFrom grid textGrob gpar
#' @importFrom picante pd prune.sample
#' @importFrom metacoder calc_taxon_abund heat_tree parse_tax_data
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap draw
#' @importFrom microbiomeMarker run_lefse plot_cladogram
#' @importFrom phyloseq phyloseq otu_table tax_table sample_data merge_phyloseq ntaxa taxa_names UniFrac
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq resultsNames results
#' @import rstudioapi
#' @importFrom stats IQR aggregate formula sd as.formula median na.exclude relevel p.adjust prcomp cor as.dist cutree hclust complete.cases
#' @importFrom utils head read.csv select.list type.convert
#' @importFrom methods show
#' @import stringr
#' @import dplyr
#' @importFrom tidyr pivot_longer pivot_wider
#' @import ggplot2
#' @import ggpubr
#' @importFrom ggrepel geom_label_repel geom_text_repel
#' @import rstatix
#' @importFrom vegan rrarefy rarecurve diversity vegdist adonis adonis2 betadisper
#' @importFrom ape rtree pcoa biplot.pcoa
#' @importFrom randomForest randomForest
#' @importFrom Boruta Boruta TentativeRoughFix
#' @importFrom grDevices dev.off png
#' @importFrom utils View write.csv stack
#' @importFrom S4Vectors unfactor
#' @rawNamespace import(crayon,except=c('%+%'))
#' @importFrom purrr reduce
#' @importFrom scales pretty_breaks
#' @importFrom limma lmFit eBayes
#' @importFrom matrixStats colQuantiles
#' @importFrom Hmisc rcorr
#' @importFrom zCompositions cmultRepl
#' @importFrom ALDEx2 aldex.clr aldex.ttest aldex.kw aldex.effect
#'
mvload <- function(path_to_folder=NULL,path_to_metadata=NA,path_to_taxa=NA,path_to_fxnl=NA,
                   autoProcess=T,combineDupes=T,combineDataSets=F) {
  assign('.loading',TRUE,envir = mvEnv)
  on.exit(assign('.loading',F,envir = mvEnv))

  do_taxa <- F
  do_fxnl <- F
  #### Load Project Directory ####
  # If a directory path is provided, check if it's valid
  if(is.null(path_to_folder)) {path_to_folder <- NA}
  if(!dir.exists(as.character(path_to_folder))) {
    message('\nSelect directory with abundance data and metadata. Results will also be saved to this directory\n')
    Sys.sleep(0.1) # To make sure it displays the above message before opening the dialogue box
    assign('project_dir',selectDirectory('Select directory with abundance data and metadata',
                                       path=dirname(getwd())),
           envir = mvEnv)
  } else {
    assign('project_dir',path_to_folder,envir = mvEnv)
  }

  if(is.null(get('project_dir',envir = mvEnv))) {
    return(message('\nERROR: No project directory selected! Exiting \n'))
  } else {
    cat('\nProject directory set to:\n',get('project_dir',envir = mvEnv),
        '\n\n All results will be stored in:\n ',file.path(get('project_dir',envir = mvEnv),'Results'))
  }

  #### Combine Datasets if Desired ####
  if(combineDataSets) {
    combined_ds_paths <- mvcombine(get('project_dir',envir = mvEnv))
    if(is.null(combined_ds_paths)) return(message('\nSelected data files could not be combined'))
    path_to_metadata <- combined_ds_paths$metadata
    path_to_taxa <- combined_ds_paths$taxa
    path_to_fxnl <- combined_ds_paths$pwys
    assign('project_dir',combined_ds_paths$project,envir = mvEnv)
  }

  #### Load Metadata ####
  if(get('.loading',envir = mvEnv)) {metadata <- loadMDFile(path_to_metadata=path_to_metadata)}
  cat('\n')

  #### Load Taxonomic Data ####
  if(get('.loading',envir = mvEnv) & ifelse(select.list(title='\nWould you like to load a taxonomic dataset?', choices=c('Yes','No'))=='Yes',T,F)) {
    taxadata <- loadTaxaFile(path_to_taxa=path_to_taxa,metadata=metadata,combineDupes=combineDupes)

    if(!is.null(taxadata)) {
      taxa_ds <- list(metadata=metadata,
                      data=taxadata,
                      features='taxa')
      do_taxa <- T
    }
  }

  #### Load Non-Taxonomic Data ####
  add_ds_list <- list()
  while(ifelse(select.list(title = '\nWould you like to add a non-taxonomic dataset?',
                           choices=c('Yes','No'))=='Yes',T,F)) {
    if(get('.loading',envir = mvEnv)) {
      fxnldata <- loadFxnlFile(path_to_fxnl=path_to_fxnl,metadata=metadata)
    }

    if(!is.null(fxnldata)) {
      cat('\n')
      data_name <- str_replace_all(readline(prompt='What kind of data is this? Give a one-word name: '),
                                   '[- +=!@#$%^&*()]','_')
      while(data_name %in% names(add_ds_list)) {
        data_name <- str_replace_all(readline(prompt='Please choose a name that has not been used: '),
                                     '[- +=!@#$%^&*()]','_')
      }
      assign(paste0(data_name,'_ds'),list(metadata=metadata,
                                          data=fxnldata,
                                          features=data_name),-1)
      add_ds_list[[data_name]] <- paste0(data_name,'_ds')
      do_fxnl <- T
    }
  }

  #### Sanity Checks and Data Packaging ####
  if(!exists('taxadata') & !exists('fxnldata')) {
    # If neither taxonomic nor functional data could be loaded, then exit
    return(cat('\n!!!\nERROR: No valid taxonomy or functional abundance data found! Exiting\n!!!'))
  }

  if(do_taxa) {
    taxa_ds <- chooseFactors(taxa_ds)
    taxa_ds <- orderGroups(taxa_ds)
    taxa_ds <- removeLowQuality(taxa_ds)
    taxa_ds <- processDataset(taxa_ds, silent = T)
    taxa_ds$results_path <- file.path(get('project_dir',envir = mvEnv),'Taxonomic Analysis')
  }

  if(do_fxnl) {
    for(ds in add_ds_list) {
      fxnl_ds <- get(ds)
      if(do_taxa) {
        fxnl_ds$metadata <- taxa_ds$metadata
        fxnl_ds$factors <- taxa_ds$factors
        fxnl_ds$active_factor <- taxa_ds$active_factor
        fxnl_ds$colors <- taxa_ds$colors
      } else if(grep(ds,add_ds_list)>1) {
        fxnl_ds$metadata <- get(add_ds_list[[1]])$metadata
        fxnl_ds$factors <- get(add_ds_list[[1]])$factors
        fxnl_ds$active_factor <- get(add_ds_list[[1]])$active_factor
        fxnl_ds$colors <- get(add_ds_list[[1]])$colors
      } else {
        fxnl_ds <- chooseFactors(fxnl_ds)
        fxnl_ds <- orderGroups(fxnl_ds)
      }
      fxnl_ds <- processDataset(fxnl_ds, silent = T)
      fxnl_ds$results_path <- file.path(get('project_dir',envir = mvEnv),
                                        paste0(capitalize(fxnl_ds$features),' Analysis'))
      assign(ds,fxnl_ds,-1)
    }
  }

  if(do_taxa | do_fxnl) cat('\n\n>>>  DATA LOADED SUCCESSFULLY!  <<<\n')

  if(do_taxa) {
    taxa_ds$name <- 'taxa_raw'
    assign('taxa_raw',taxa_ds,1)
    if(autoProcess) {
      taxa_ds$name <- 'taxa_proc'
      taxa_ds <- transData(taxa_ds, trans_method = 'clr', silent = T)
      taxa_ds <- filterNAs(taxa_ds,ranks = 'domain',silent = T)
      taxa_ds <- processDataset(taxa_ds,temp = T)

      assign('taxa_proc',taxa_ds,1)
    }
    assign('active_dataset',taxa_ds,envir = mvEnv)
    print(taxa_ds)
  }

  if(do_fxnl) {
    for(ds in add_ds_list) {
      fxnl_ds <- get(ds,inherits = F)
      fxnl_ds$name <- paste0(fxnl_ds$features,'_raw')
      assign(paste0(fxnl_ds$features,'_raw'),fxnl_ds,1)
      if(autoProcess) {
        fxnl_ds$name <- paste0(fxnl_ds$features,'_proc')
        fxnl_ds <- filterLowPrev(fxnl_ds, silent = T)
        fxnl_ds <- filterLowRelAbun(fxnl_ds, silent = T)
        fxnl_ds <- scaleSamples(fxnl_ds, scaling = 'sum', silent = T)
        fxnl_ds <- processDataset(fxnl_ds,temp = T)

        assign(paste0(fxnl_ds$features,'_proc'),fxnl_ds,1)
        print(fxnl_ds)
      }
    }
    if(!do_taxa) {
      if(autoProcess) {
        assign('active_dataset',
               get(paste0(sub('_([^_]*)$','',add_ds_list[[1]]),'_proc')),
               envir = mvEnv)
        print(get(paste0(sub('_([^_]*)$','',add_ds_list[[1]]),'_proc')))
      } else {
        assign('active_dataset',
               get(paste0(sub('_([^_]*)$','',add_ds_list[[1]]),'_raw')),
               envir = mvEnv)
        print(get(paste0(sub('_([^_]*)$','',add_ds_list[[1]]),'_raw')))
      }
    }
  }

  #### Output Messages ####
  if(do_taxa) cat('\nRaw taxonomy dataset is stored in "taxa_raw"')
  if(do_fxnl) for(ds in add_ds_list) cat(paste0('\nRaw functional dataset is stored in "',
                                                sub('_([^_]*)$','',ds),'_raw"'))
  if(autoProcess) {
    if(do_taxa) cat('\nAuto-processed taxonomy dataset is stored in "taxa_proc"')
    if(do_fxnl) for(ds in add_ds_list) cat(paste0('\nAuto-processed functional dataset is stored in "',
                                                  sub('_([^_]*)$','',ds),'_proc"'))
    if(do_taxa) cat('\n\n  <|> Active Dataset: "taxa_proc" <|>\n\n')
    else if(do_fxnl) cat(paste0('\n\n  <|> Active Dataset: "',
                                sub('_([^_]*)$','',add_ds_list[[1]]),'_proc" <|>\n\n'))
  } else {
    if(do_taxa) cat('\n\n  <|> Active Dataset: "taxa_raw" <|>\n\n')
    else if(do_fxnl) cat('\n\n  <|> Active Dataset: "fxnl_raw" <|>\n\n')
  }

  return(cat('\nType "mvhelp()" then press Enter to see the available functions!\n\n'))
}

#' MicroVis ggplot theme
#'
#' @return NULL
#' @export
#'
theme_mv <- function() {
  theme_pubr()+
    theme(plot.title = element_text(hjust = 0.5, size=25),
          axis.title = element_text(size=25),
          axis.text = element_text(size=22),
          legend.title = element_blank(),
          legend.text = element_text(size=20, margin=margin(r=10)),
          legend.key.size = unit('5','mm'))
}
