#################
# Author: Thomas Hart
# Contact: thomas.hart@columbia.edu
#################

################# START updates list 
# 08/22/12016 @ 1848
#   1. file completed
################# END updatse list 

################# START source 

#source(file = "ConstructPatientFilePathTree.R")
#source(file = "ReformatFileNames.R")

################# END source 

################# START announcer 
cat(
  "TransformPhenotype.R contains:\n",
  "\tTransformPhenotype(patient_id,patient_image_id,file_path_tree,thresholds,phenotype_possibilities)\n",
  "\tTransformPhenotypesoinPatientFilePathTree(file_path_tree,file_path_tree,thresholds,phenotype_possibilities)\n",sep=""
)
################# END announcer 

CorrectTransformedPhenotype = function(phenotype_possibilities,transformed_phenotype) {
  possibilities = phenotype_possibilities[[transformed_phenotype[1]]]
  return( c(transformed_phenotype[1] , transformed_phenotype[ transformed_phenotype %in% possibilities ]) )
}

# ### START check CorrectTransformedPhenotype
# phenotype_possibilities = list(
#   "CD3+"=c("CD8+","CD8-","Ki67+","Ki67-","HLA-DR+","HLA-DR-"),
#   "CD8+"=c("Ki67+","Ki67-","HLA-DR+","HLA-DR-","CD68+","CD68-"),
#   "Sox10+"=c("HLA-DR+","HLA-DR-","Ki67+","Ki67-","CD68+","CD68-"),
#   "CD68+"=c("CD8+","CD8-","HLA-DR+","HLA-DR-"),
#   "Other" = c("HLA-DR+","HLA-DR-")
# )
# string = "Other CD8- HLA-DR- Ki67-"
# split = strsplit(x = string,split = " ")[[1]]
# CorrectTransformedPhenotype(phenotype_possibilities,transformed_phenotype = split)
# 
# string = "CD68+ Ki67- HLA-DR+"
# split = strsplit(x = string,split = " ")[[1]]
# CorrectTransformedPhenotype(phenotype_possibilities,transformed_phenotype = split)
# 
# string = "Sox10+ HLA-DR+ Ki67+ CD8+"
# split = strsplit(x = string,split = " ")[[1]]
# CorrectTransformedPhenotype(phenotype_possibilities,transformed_phenotype = split)
# ### END check CorrectTransformedPhenotype

TransformPhenotype = function(patient_id,patient_image_id,file_path_tree,thresholds,phenotype_possibilities) {
  # get cell seg data
  cell_seg_data_file_path = GetCellSegDataFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
  
  # read in cell seg data
  cell_seg_data = read.csv(file = cell_seg_data_file_path,sep = "\t")
  
  # get score data
  score_data_file_paths = GetScoreDataFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
  
  # read in score_data
  score_datas = list()
  for ( score_data_file_path in score_data_file_paths ) {
    temp = read.csv(score_data_file_path,sep = "\t")
    score_datas[[score_data_file_path]] = temp;
  }
  
  # extract thresholds from score_data
  threshold_values = array(data = 0,dim = length(thresholds[,1]))
  for ( score_data in score_datas ) {
    which_thresholds_in_score_data_names = which( thresholds[,1] %in% names(score_data) )
    threshold_values[ which_thresholds_in_score_data_names ] = as.numeric( score_data[ 1 , thresholds[,1][ which_thresholds_in_score_data_names ] ] )
  }
  cat("\t\t\tthreshold values:\n")
  for ( row in 1:length(thresholds[,1]) ) {
    cat(sprintf("\t\t\t%30s %7f",thresholds[row,1],threshold_values[row]),"\n")
  }
  
  # concatenate phenotypes
  cell_seg_data$Transformed_Phenotype = as.character(cell_seg_data$Phenotype)
  cell_seg_data$Transformed_Phenotype[ cell_seg_data$Phenotype=="" ] = "Other"
  for ( i in 1:length(threshold_values) ) {
    
    # concatenate
    positive_indices_to_cat = which( cell_seg_data[ thresholds[i,2] ] >= threshold_values[i]  )
    for ( index in positive_indices_to_cat ) {
      cell_seg_data$Transformed_Phenotype[index] = paste(c(cell_seg_data$Transformed_Phenotype[index],thresholds[i,3]),collapse = " ")
      string = cell_seg_data$Transformed_Phenotype[index];
      split = strsplit(x = string,split = " ")[[1]]
      cell_seg_data$Transformed_Phenotype[index] = paste(CorrectTransformedPhenotype(phenotype_possibilities,transformed_phenotype = split),collapse = " ")
    }
    negative_indices_to_cat = which( cell_seg_data[ thresholds[i,2] ] < threshold_values[i]  )
    for ( index in negative_indices_to_cat ) {
      cell_seg_data$Transformed_Phenotype[index] = paste(c(cell_seg_data$Transformed_Phenotype[index],thresholds[i,4]),collapse = " ")
      string = cell_seg_data$Transformed_Phenotype[index];
      split = strsplit(x = string,split = " ")[[1]]
      cell_seg_data$Transformed_Phenotype[index] = paste(CorrectTransformedPhenotype(phenotype_possibilities,transformed_phenotype = split),collapse = " ")
    }
  }
  
  # write data
  write.table(x = cell_seg_data,file = cell_seg_data_file_path,sep = "\t",row.names = F)
}

# ### START TransformPhenotype check
# rm(list=ls())
# setwd(dir = "G:\\Script Rebuild\\Transform Script/")
# # thresholds
# thresholds = rbind(
#   #c( "threshold value in _score_data", "_cell_seg_data mean data", "if >=", "if <" )
#   c("CD68_Opal_520_Threshold","Membrane_CD68_Opal_520_Mean_Normalized_Counts_Total_Weighting","CD68+","CD68-"),
#   c("CD8_Opal_620_Threshold","Membrane_CD8_Opal_620_Mean_Normalized_Counts_Total_Weighting","CD8+","CD8-"),
#   c("HLA_DR_Opal_570_Threshold","Membrane_HLA_DR_Opal_570_Mean_Normalized_Counts_Total_Weighting","HLA-DR+","HLA-DR-"),
#   c("Ki67_Opal_540_Threshold","Nucleus_Ki67_Opal_540_Mean_Normalized_Counts_Total_Weighting","Ki67+","Ki67-")
# )
# # phenotype possibilities
# phenotype_possibilities = list(
#   "CD3+"=c("CD8+","CD8-","Ki67+","Ki67-","HLA-DR+","HLA-DR-","CD68+","CD68-"),
#   "CD3-CD8+"=c("Ki67+","Ki67-","HLA-DR+","HLA-DR-","CD68+","CD68-"),
#   "Sox10+"=c("Ki67+","Ki67-","HLA-DR+","HLA-DR-","CD68+","CD68-"),
#   "CD68+"=c("CD8+","CD8-","HLA-DR+","HLA-DR-","Ki67+","Ki67-"),
#   "Other" = c("HLA-DR+","HLA-DR-")
# )
# # sink dir
# sink_directory = "G:\\Script Rebuild/Training_Data_P/"
# # construct patient file tree
# file_path_tree = ConstructPatientFilePathTree(dir_name = sink_directory)
# # patient id
# patient_id = "cu-04-06";
# # patient_image_id
# patient_image_id = "cu-04-06 D_1"
# # run
# TransformPhenotype(patient_id = patient_id,patient_image_id = patient_image_id,
#                    file_path_tree = file_path_tree,thresholds = thresholds,
#                    phenotype_possibilities = phenotype_possibilities)
# ### END TransformPhenotype check

TransformPhenotypesoinPatientFilePathTree = function(file_path_tree,thresholds,phenotype_possibilities) {
  # read in, correct and rewrite each file
  for ( patient_id in names(file_path_tree) ) {
    cat("\tpatient_id: ",patient_id,"\n")
    for ( patient_image_id in names(file_path_tree[[patient_id]]) ) {
      cat("\t\t image: ",patient_image_id,"\n")
      TransformPhenotype(patient_id = patient_id,patient_image_id = patient_image_id,
                         file_path_tree = file_path_tree,thresholds = thresholds,
                         phenotype_possibilities = phenotype_possibilities)
    }
  }
}


