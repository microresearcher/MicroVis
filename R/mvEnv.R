mvEnv <- new.env()
mvDefaults <- new.env()

assign('.loading',TRUE,envir = mvEnv)
assign('warning_list',list(),envir = mvEnv)
assign('autosave',F,envir = mvEnv)
assign('offerSave',T,envir = mvEnv)
assign('autoNameResults',T,envir = mvEnv)
assign('detailed_taxa_names',F,envir = mvEnv)
assign('taxaRanks',c('domain','phylum','class','order','family','genus','species'),envir = mvEnv)
assign('keepSigFisher',T,envir = mvEnv)
assign('imageType','png',envir = mvEnv)
# When forceStats is set to true, values will be slightly adjusted for
#   groups with the same mean prior to statistical testing
assign('forceStats',F,envir = mvEnv)
assign('adjust_value',0.001,envir = mvEnv)

#### DEFAULT COLORS ####
assign('defCols',c('#6700b5', # Purple
                 '#00a6cf', # Baby blue
                 '#d12121', # Red
                 '#00a60c', # Green
                 '#b500a9', # Pink
                 '#d1a400', # Yellow
                 '#4432e3', # Blue
                 '#db8e00', # Orange
                 '#00d18f', # Mint
                 '#ed594c', # Coral
                 '#310057', # Dark purple
                 '#006a85', # Teal
                 '#6e1212', # Maroon
                 '#005c07', # Dark green
                 '#6e0066', # Dark pink?
                 '#8c6e00', # Golden
                 '#0a0063', # Navy
                 '#6b4500', # Brown
                 '#006948', # Dark... mint?
                 '#8a423b' # Burnt... coral?
                 ),envir = mvEnv)

### DEFAULT SAMPLE FILTERING CUTOFF
assign('rthresh',10000,envir = mvDefaults) # Minimum read count threshold for each sample

#### DEFAULT FEATURE FILTERING PARAMETERS ####
assign('filtering.defaults',list(min_totabun=10,
                                 low_abun=list(min_abun=1,min_prop=20),
                                 min_relabun=0.0001,
                                 min_prevalence=1,
                                 low_var_percentile=5),
       envir = mvDefaults)

### ONLY FOR MICROVIS DEVELOPMENT ###
assign('.debug',F,envir = mvEnv)
assign('mga',F,envir = mvEnv)
