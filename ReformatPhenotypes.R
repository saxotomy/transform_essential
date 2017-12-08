#################
# Author: Thomas Hart
# Contact: thomas.hart@columbia.edu
#################

################# START updates list 
# 08/22/12016 @ 1846
#   1. file created
################# END updatse list 

################# START source 

#source(file = "ConstructPatientFilePathTree.R")
#source(file = "ReformatPhenotypes.R")

################# END source 

################# START announcer 
cat(
  "ReformatPhenotypes contains:\n",
  "\tFormatPhenotypes(data_names)\n",
  "\tFormatPhenotypesinPatientFilePathTree(patient_file_path_tree)\n",sep=""
)
################# END announcer 

FormatPhenotypes = function(phenotype,label_corrections) {
  phenotype = as.character(phenotype)
  for ( label_correction in label_corrections ) {
    selection = which(phenotype==label_correction[1])
    phenotype[ selection ] = label_correction[2]
  }
  return(phenotype)
}

# # START FixPhenotypeLabels check
# label_corrections = list(
#   c("CD3","CD3+"),
#   c("CD68","CD68+"),
#   c("Ki67","Ki67+"),
#   c("HLA-DR","HLA-DR+"),
#   c("CD8","CD3-CD8+"),
#   c("other","Other")
# )
# Phenotypes = c("CD8","CD3","CD68","CD28","Ki67","Ki67","HLA-DR","HLA-DR","other","other")
# FormatPhenotypes(phenotype = Phenotypes,label_corrections = label_corrections)
# check = c("CD8+","CD3+","CD68+","CD28","Ki67+","Ki67+","HLA-DR+","HLA-DR+","Other","Other")
# # END FixPhenotypeLabels check

FormatPhenotypesinPatientFilePathTree = function(file_path_tree,label_corrections) {
  # read in, correct and rewrite each file
  for ( patient_id in names(file_path_tree) ) {
    cat("\tpatient_id: ",patient_id,"\n")
    for ( patient_image_id in names(file_path_tree[[patient_id]]) ) {
      cat("\t\t image_id: ",patient_image_id,"\n")
      
      # format names for cell seg data
      cat("\t\t\tcorrecting phenotype labels : _cell_seg_data ... ")
      cell_seg_data_file_path = GetCellSegDataFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      cell_seg_data = read.csv(file = cell_seg_data_file_path,sep = "\t")
      cell_seg_data$Phenotype = FormatPhenotypes(phenotype = cell_seg_data$Phenotype,label_corrections = label_corrections)
      write.table(x = cell_seg_data,file = cell_seg_data_file_path,sep = "\t",row.names = FALSE)
      cat("done\n")
      
    }
  }
}

# ### START FormatPhenotypesinPatientFilePathTree check\
# rm(list=ls())
# setwd(dir = "G:\\Script Rebuild\\Transform Script/")
# # label corrections 
# label_corrections = list(
#   c("CD3","CD3+"),
#   c("CD68","CD68+"),
#   c("Ki67","Ki67+"),
#   c("HLA-DR","HLA-DR+"),
#   c("CD8","CD3-CD8+"),
#   c("other","Other")
# )
# # sink dir
# sink_directory = "G:\\Script Rebuild/Training_Data_P/"
# # construct patient file tree
# file_path_tree = ConstructPatientFilePathTree(dir_name = sink_directory)
# FormatPhenotypesinPatientFilePathTree(file_path_tree = file_path_tree,label_corrections = label_corrections)
# ### END FormatPhenotypesinPatientFilePathTree check





