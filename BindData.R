
#################
# Author: Thomas Hart
# Contact: thomas.hart@columbia.edu
#################

################# START updates list 
# 2016/09/06 @ 18:10
#   1. file created
################# END updatse list 

# rm(list=ls())

# ################# START source 
# 
# setwd("F:\\Script Rebuild\\Transform Script\\")
# source(file = "ConstructPatientFilePathTree.R")
# 
# ################# END source 

################# START announcer 
cat(
  "BindData.R contains:\n",
  "\tBindAllDataInDir( file_path )\n"
)
################# END announcer 

BindAllDataInDir = function(file_path) {
  file_paths = dir(path = file_path,full.names = TRUE)
  data_list = list()
  for ( file_path in file_paths ) {
    cat("\treading in [ ",file_path," ] ... ")
    data = read.csv(file = file_path,sep = "\t",check.names=FALSE)
    data_list[[file_path]] = data
    cat("done\n")
  }
  cat("\tbinding data_list[[",1,"]] ... ")
  bound_data = data_list[[1]]
  cat("done\n")
  for ( i in 2:length(data_list) ) {
    cat("\tbinding data_list[[",i,"]] ... ")
    bound_data = rbind(bound_data,data_list[[i]])
    cat("done\n")
  }
  return( bound_data )
}
