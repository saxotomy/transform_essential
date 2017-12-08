#################
# Author: Thomas Hart
# Contact: thomas.hart@columbia.edu
#################

################# START updates list 
# 08/22/12016 @ 1442
#   1. file created
################# END updatse list 

################# START source 

#source(file = "ConstructPatientFilePathTree.R")

################# END source 

################# START announcer 
cat(
  "ReformatFileNames.R contains:\n",
  "\tFormatDataNames(data_names)\n",
  "\tFormatVariableNamesinPatientFilePathTree(patient_file_path_tree)\n",sep=""
)
################# END announcer 

FormatVariableNames = function(data_names) {
  updated_data_names = sub(pattern = "[.]$",replacement = "",x = data_names)
  return( gsub(pattern = "[.]+",replacement = "_",x = updated_data_names) )
}

# ### START check FormatVariableNames
# data_names = c("Path", "Sample.Name", "Tissue.Category", "Cell.ID", "Total.Cells", "Tissue.Category.Area..pixels.")
# FormatVariableNames(data_names = data_names)
# check = c("Path", "Sample_Name", "Tissue_Category", "Cell_ID", "Total_Cells", "Tissue_Category_Area_pixels")
# ### END check FormatVariableNames

FormatVariableNamesinPatientFilePathTree = function(file_path_tree) {
  
  # read in, correct and rewrite each file
  for ( patient_id in names(file_path_tree) ) {
    cat("\tpatient_id: ",patient_id,"\n")
    for ( patient_image_id in names(file_path_tree[[patient_id]]) ) {
      cat("\t\timage_id: ",patient_image_id,"\n")
      
      # format names for cell seg data
      cat("\t\t\tcorrecting names : _cell_seg_data ... ")
      cell_seg_data_file_path = GetCellSegDataFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      cell_seg_data = read.csv(file = cell_seg_data_file_path,sep = "\t")
      names(cell_seg_data) = FormatVariableNames(data_names = names(cell_seg_data))
      write.table(x = cell_seg_data,file = cell_seg_data_file_path,sep = "\t",row.names = FALSE)
      cat("done\n")
      
      # format names for cell seg data summary
      cat("\t\t\tcorrecting names : _cell_seg_data_summary  ... ")
      cell_seg_data_summary_file_path = GetCellSegDataSummaryFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      cell_seg_data_summary = read.csv(file = cell_seg_data_summary_file_path,sep = "\t")
      names(cell_seg_data_summary) = FormatVariableNames(data_names = names(cell_seg_data_summary))
      write.table(x = cell_seg_data_summary,file = cell_seg_data_summary_file_path,sep = "\t",row.names = FALSE)
      cat("done\n")
      
      # format names for tissue seg data
      cat("\t\t\tcorrecting names : _tissue_seg_data ... ")
      tissue_seg_data_file_path = GetTissueSegDataFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      tissue_seg_data = read.csv(file = tissue_seg_data_file_path,sep = "\t")
      names(tissue_seg_data) = FormatVariableNames(data_names = names(tissue_seg_data))
      write.table(x = tissue_seg_data,file = tissue_seg_data_file_path,sep = "\t",row.names = FALSE)
      cat("done\n")
      
      # format names for tissue seg data summary
      cat("\t\t\tcorrecting names : _tissue_seg_data_summary ... ")
      tissue_seg_data_summary_file_path = GetTissueSegDataSummaryFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      tissue_seg_data_summary = read.csv(file = tissue_seg_data_summary_file_path,sep = "\t")
      names(tissue_seg_data_summary) = FormatVariableNames(data_names = names(tissue_seg_data_summary))
      write.table(x = tissue_seg_data_summary,file = tissue_seg_data_summary_file_path,sep = "\t",row.names = FALSE)
      cat("done\n")
      
      # format names for score data
      cat("\t\t\tcorrecting names : _score_data ... ")
      score_data_file_paths = GetScoreDataFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      for ( score_data_file_path in score_data_file_paths ) {
        score_data = read.csv(file = score_data_file_path,sep = "\t")
        names(score_data) = FormatVariableNames(data_names = names(score_data))
        write.table(x = score_data,file = score_data_file_path,sep = "\t",row.names = FALSE)  
      }
      cat("done\n")    
    }
  }
}


# ### START FormatVariableNamesinPatientFilePathTree check
# rm(list=ls())
# setwd(dir = "G:\\Script Rebuild\\Transform Script/")
# # sink dir
# sink_directory = "G:\\Script Rebuild/Training_Data_P/"
# # construct patient file tree
# file_path_tree = ConstructPatientFilePathTree(dir_name = sink_directory)
# FormatVariableNamesinPatientFilePathTree(file_path_tree)
# ### END FormatVariableNamesinPatientFilePathTree check





