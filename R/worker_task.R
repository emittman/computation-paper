# i = identify the scenario 
# code below produce data, analysis and summary for ONE scenario

ss  = Sys.getenv("SLURM_ARRAY_TASK_ID") 
i = as.numeric(ss) + 1

# 0 ) libraries and functions
pkgs <- c('dplyr', 'cudarpackage')
lapply(pkgs, library,  character.only = TRUE, quietly=TRUE )

# inputs
load('targets.RData')
source('workflow_fun.R')

# ==================================================
# 1) Create a dataset, run the analysis and its summary
# 1.1 Create data set and save it
eval( parse(text= paste(my_plan$target[i], my_plan$command[i], sep=' <- ')  ) )
saveRDS(get(my_plan$target[i]) , file = paste('arrayOutput/', my_plan$target[i], '.rds', sep='') )

# 1.2 Run all analyses that use this dataset
id <- grep(my_plan$target[i], my_plan$command) 
for (j in id) {
  eval( parse(text= paste(my_plan$target[j], my_plan$command[j], sep=' <- ')  ) )
  saveRDS(get(my_plan$target[j]), file = paste('arrayOutput/', my_plan$target[j], '.rds', sep='') )
}