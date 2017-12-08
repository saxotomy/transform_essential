################# 
# Author: Thomas Hart
# Contact: thomas.hart@columbia.edu
################# 

gc()
rm(list=ls())

setwd(dir = "F:\\Script Rebuild\\Transform Essential Scripts")
################# START source ################# 

cat("Loading required functions ... \n")
source(file = "ConstructPatientFilePathTree.R")
source(file = "ReformatFileNames.R")
source(file = "ReformatPhenotypes.R")
source(file = "TransformPhenotypes.R")
source(file = "ComputeImageDensities.R")
source(file = "ComputeNearestNeighbors.R")
source(file = "BindData.R")
cat("done\n\n\n")

################# END source ################# 

################# START main ################# 
cat('---------------------------- BEGIN MAIN -----------------------------\n')

# source directory
source_directory = "F:\\Script Rebuild\\Transform Essential Scripts\\Run/Datasets/"

# sink directory
sink_directory = "F:\\Script Rebuild\\Transform Essential Scripts\\Run/Datasets_P/"

# copy data from source directory to sink directory
cat("Copying: [",source_directory,"]->[",sink_directory,"] ...\n")
file_paths = list.files(path = source_directory,full.names = TRUE)
for ( file_path in file_paths ) {
  cat("\tCopying: ",sub(pattern = ".*[\\/]",replacement = "",x = file_path),"\n")
  file.copy(from = file_path,to = sink_directory,recursive = TRUE)
}
cat("done\n")

# build file tree from files in the sink directory
cat("Constructing file_path_tree ...\n")
file_path_tree = ConstructPatientFilePathTree(dir_name = sink_directory)
cat("done\n")

# reformat file names in all file tree files
cat("Correcting variable names ...\n")
FormatVariableNamesinPatientFilePathTree(file_path_tree = file_path_tree)
cat("done\n")

# reformat phenotype labels in _cell_seg_data.txts
label_corrections = list(
  c("CD3","CD3+"),
  c("CD68","CD68+"),
  c("Ki67","Ki67+"),
  c("HLA-DR","HLA-DR+"),
  c("CD8","CD3-CD8+"),
  c("Tumor","Sox10+"),
  c("other","Other")
)
cat("Correcting phenotypes labels ... ")
FormatPhenotypesinPatientFilePathTree(file_path_tree = file_path_tree,label_corrections = label_corrections)
cat("done\n")

# transform phenotypes in _cell_seg_data.txts using thresholds from _score_data.txt
thresholds = rbind(
  #c( "threshold value in _score_data", "_cell_seg_data mean data", "if >=", "if <" )
  c("CD68_Opal_520_Threshold","Membrane_CD68_Opal_520_Mean_Normalized_Counts_Total_Weighting","CD68+","CD68-"),
  c("CD8_Opal_620_Threshold","Membrane_CD8_Opal_620_Mean_Normalized_Counts_Total_Weighting","CD8+","CD8-"),
  c("HLA_DR_Opal_570_Threshold","Membrane_HLA_DR_Opal_570_Mean_Normalized_Counts_Total_Weighting","HLA-DR+","HLA-DR-"),
  c("Ki67_Opal_540_Threshold","Nucleus_Ki67_Opal_540_Mean_Normalized_Counts_Total_Weighting","Ki67+","Ki67-")
)
# phenotype possibilities
phenotype_possibilities = list(
  "CD3+"=c("CD8+","CD8-","Ki67+","Ki67-","HLA-DR+","HLA-DR-","CD68+","CD68-"),
  "CD3-CD8+"=c("Ki67+","Ki67-","HLA-DR+","HLA-DR-","CD68+","CD68-"),
  "Sox10+"=c("Ki67+","Ki67-","HLA-DR+","HLA-DR-","CD68+","CD68-"),
  "CD68+"=c("CD8+","CD8-","HLA-DR+","HLA-DR-","Ki67+","Ki67-"),
  "Other" = c("HLA-DR+","HLA-DR-")
)
cat("Transforming phenotypes ...\n")
TransformPhenotypesoinPatientFilePathTree(file_path_tree = file_path_tree,thresholds = thresholds,phenotype_possibilities = phenotype_possibilities)
cat("done\n")

# count phenotypes in _cell_seg_data.txt and compute densities using area from _tissue_seg_data_summary.txt
# output counts to _tumor_densities_by_image.txt or _stroma_densities_by_image.txt
phenotypes_to_count = list(
  c("CD3+"),
#  c("CD3+","CD68-"),
#  c("CD3+","CD8-"),
#  c("CD3+","HLA-DR-"),
  c("CD3+","Ki67-"),
#  c("CD3+","CD68+"),
  c("CD3+","CD8+"),
  c("CD3+","HLA-DR+"),
  c("CD3+","Ki67+"),
#  c("CD3+","CD68-","CD8-"),
#  c("CD3+","CD68-","HLA-DR-"),
#  c("CD3+","CD68-","Ki67-"),
#  c("CD3+","CD8-","HLA-DR-"),
#  c("CD3+","CD8-","Ki67-"),
#  c("CD3+","HLA-DR-","Ki67-"),
#  c("CD3+","CD68-","CD8+"),
#  c("CD3+","CD68-","HLA-DR+"),
#  c("CD3+","CD68-","Ki67+"),
#  c("CD3+","CD8-","HLA-DR+"),
#  c("CD3+","CD8-","Ki67+"),
#  c("CD3+","HLA-DR-","Ki67+"),
#  c("CD3+","CD68+","CD8-"),
#  c("CD3+","CD68+","HLA-DR-"),
#  c("CD3+","CD68+","Ki67-"),
  c("CD3+","CD8+","HLA-DR-"),
  c("CD3+","CD8+","Ki67-"),
#  c("CD3+","HLA-DR+","Ki67-"),
#  c("CD3+","CD68+","CD8+"),
#  c("CD3+","CD68+","HLA-DR+"),
#  c("CD3+","CD68+","Ki67+"),
  c("CD3+","CD8+","HLA-DR+"),
  c("CD3+","CD8+","Ki67+"),
  c("CD3+","HLA-DR+","Ki67+"),
#  c("CD3+","CD68-","CD8-","HLA-DR-"),
#  c("CD3+","CD68-","CD8-","Ki67-"),
#  c("CD3+","CD68-","HLA-DR-","Ki67-"),
#  c("CD3+","CD8-","HLA-DR-","Ki67-"),
#  c("CD3+","CD68-","CD8-","HLA-DR+"),
#  c("CD3+","CD68-","CD8-","Ki67+"),
#  c("CD3+","CD68-","HLA-DR-","Ki67+"),
#  c("CD3+","CD8-","HLA-DR-","Ki67+"),
#  c("CD3+","CD68-","CD8+","HLA-DR-"),
#  c("CD3+","CD68-","CD8+","Ki67-"),
#  c("CD3+","CD68-","HLA-DR+","Ki67-"),
#  c("CD3+","CD8-","HLA-DR+","Ki67-"),
#  c("CD3+","CD68-","CD8+","HLA-DR+"),
#  c("CD3+","CD68-","CD8+","Ki67+"),
#  c("CD3+","CD68-","HLA-DR+","Ki67+"),
#  c("CD3+","CD8-","HLA-DR+","Ki67+"),
#  c("CD3+","CD68+","CD8-","HLA-DR-"),
#  c("CD3+","CD68+","CD8-","Ki67-"),
#  c("CD3+","CD68+","HLA-DR-","Ki67-"),
#  c("CD3+","CD8+","HLA-DR-","Ki67-"),
#  c("CD3+","CD68+","CD8-","HLA-DR+"),
#  c("CD3+","CD68+","CD8-","Ki67+"),
#  c("CD3+","CD68+","HLA-DR-","Ki67+"),
#  c("CD3+","CD8+","HLA-DR-","Ki67+"),
#  c("CD3+","CD68+","CD8+","HLA-DR-"),
#  c("CD3+","CD68+","CD8+","Ki67-"),
#  c("CD3+","CD68+","HLA-DR+","Ki67-"),
#  c("CD3+","CD8+","HLA-DR+","Ki67-"),
#  c("CD3+","CD68+","CD8+","HLA-DR+"),
#  c("CD3+","CD68+","CD8+","Ki67+"),
#  c("CD3+","CD68+","HLA-DR+","Ki67+"),
#  c("CD3+","CD8+","HLA-DR+","Ki67+"),
#  c("CD3+","CD68-","CD8-","HLA-DR-","Ki67-"),
#  c("CD3+","CD68-","CD8-","HLA-DR-","Ki67+"),
#  c("CD3+","CD68-","CD8-","HLA-DR+","Ki67-"),
#  c("CD3+","CD68-","CD8-","HLA-DR+","Ki67+"),
#  c("CD3+","CD68-","CD8+","HLA-DR-","Ki67-"),
#  c("CD3+","CD68-","CD8+","HLA-DR-","Ki67+"),
#  c("CD3+","CD68-","CD8+","HLA-DR+","Ki67-"),
  # c("CD3+","CD68-","CD8+","HLA-DR+","Ki67+"),
  # c("CD3+","CD68+","CD8-","HLA-DR-","Ki67-"),
  # c("CD3+","CD68+","CD8-","HLA-DR-","Ki67+"),
  # c("CD3+","CD68+","CD8-","HLA-DR+","Ki67-"),
  # c("CD3+","CD68+","CD8-","HLA-DR+","Ki67+"),
  # c("CD3+","CD68+","CD8+","HLA-DR-","Ki67-"),
  # c("CD3+","CD68+","CD8+","HLA-DR-","Ki67+"),
  # c("CD3+","CD68+","CD8+","HLA-DR+","Ki67-"),
  # c("CD3+","CD68+","CD8+","HLA-DR+","Ki67+"),
  c("CD3-CD8+"),
# c("CD3-CD8+","CD68-"),
#  c("CD3-CD8+","HLA-DR-"),
# c("CD3-CD8+","Ki67-"),
  c("CD3-CD8+","CD68+"),
  c("CD3-CD8+","HLA-DR+"),
# c("CD3-CD8+","Ki67+"),
# c("CD3-CD8+","CD68-","HLA-DR-"),
# c("CD3-CD8+","CD68-","Ki67-"),
# c("CD3-CD8+","HLA-DR-","Ki67-"),
# c("CD3-CD8+","CD68-","HLA-DR+"),
# c("CD3-CD8+","CD68-","Ki67+"),
# c("CD3-CD8+","HLA-DR-","Ki67+"),
# c("CD3-CD8+","CD68+","HLA-DR-"),
# c("CD3-CD8+","CD68+","Ki67-"),
# c("CD3-CD8+","HLA-DR+","Ki67-"),
# c("CD3-CD8+","CD68+","HLA-DR+"),
# c("CD3-CD8+","CD68+","Ki67+"),
# c("CD3-CD8+","HLA-DR+","Ki67+"),
# c("CD3-CD8+","CD68-","HLA-DR-","Ki67-"),
# c("CD3-CD8+","CD68-","HLA-DR-","Ki67+"),
# c("CD3-CD8+","CD68-","HLA-DR+","Ki67-"),
# c("CD3-CD8+","CD68-","HLA-DR+","Ki67+"),
# c("CD3-CD8+","CD68+","HLA-DR-","Ki67-"),
# c("CD3-CD8+","CD68+","HLA-DR-","Ki67+"),
# c("CD3-CD8+","CD68+","HLA-DR+","Ki67-"),
# c("CD3-CD8+","CD68+","HLA-DR+","Ki67+"),
  c("Sox10+"),
  c("Sox10+","CD68-"),
  c("Sox10+","HLA-DR-"),
  c("Sox10+","Ki67-"),
  c("Sox10+","CD68+"),
  c("Sox10+","HLA-DR+"),
  c("Sox10+","Ki67+"),
  c("Sox10+","CD68-","HLA-DR-"),
  c("Sox10+","CD68-","Ki67-"),
  c("Sox10+","HLA-DR-","Ki67-"),
  c("Sox10+","CD68-","HLA-DR+"),
  c("Sox10+","CD68-","Ki67+"),
  c("Sox10+","HLA-DR-","Ki67+"),
  c("Sox10+","CD68+","HLA-DR-"),
  c("Sox10+","CD68+","Ki67-"),
  c("Sox10+","HLA-DR+","Ki67-"),
  c("Sox10+","CD68+","HLA-DR+"),
  c("Sox10+","CD68+","Ki67+"),
  c("Sox10+","HLA-DR+","Ki67+"),
  c("Sox10+","CD68-","HLA-DR-","Ki67-"),
  c("Sox10+","CD68-","HLA-DR-","Ki67+"),
  c("Sox10+","CD68-","HLA-DR+","Ki67-"),
  c("Sox10+","CD68-","HLA-DR+","Ki67+"),
  c("Sox10+","CD68+","HLA-DR-","Ki67-"),
  c("Sox10+","CD68+","HLA-DR-","Ki67+"),
  c("Sox10+","CD68+","HLA-DR+","Ki67-"),
  c("Sox10+","CD68+","HLA-DR+","Ki67+"),
  c("CD68+"),
  c("CD68+","CD8-"),
  c("CD68+","HLA-DR-"),
  c("CD68+","Ki67-"),
  c("CD68+","CD8+"),
  c("CD68+","HLA-DR+"),
  c("CD68+","Ki67+"),
  c("CD68+","CD8-","HLA-DR-"),
  c("CD68+","CD8-","Ki67-"),
  c("CD68+","HLA-DR-","Ki67-"),
  c("CD68+","CD8-","HLA-DR+"),
  c("CD68+","CD8-","Ki67+"),
  c("CD68+","HLA-DR-","Ki67+"),
  c("CD68+","CD8+","HLA-DR-"),
  c("CD68+","CD8+","Ki67-"),
  c("CD68+","HLA-DR+","Ki67-"),
  c("CD68+","CD8+","HLA-DR+"),
  c("CD68+","CD8+","Ki67+"),
  c("CD68+","HLA-DR+","Ki67+"),
  c("CD68+","CD8-","HLA-DR-","Ki67-"),
  c("CD68+","CD8-","HLA-DR-","Ki67+"),
  c("CD68+","CD8-","HLA-DR+","Ki67-"),
  c("CD68+","CD8-","HLA-DR+","Ki67+"),
  c("CD68+","CD8+","HLA-DR-","Ki67-"),
  c("CD68+","CD8+","HLA-DR-","Ki67+"),
  c("CD68+","CD8+","HLA-DR+","Ki67-"),
  c("CD68+","CD8+","HLA-DR+","Ki67+"),
  c("Other"),
  c("Other","HLA-DR-"),
  c("Other","HLA-DR+")
)

# cat("Counting phenotypes by image ... \n")
# stroma_output_data_by_image_file_path = "F:\\Script Rebuild\\Transform Essential Scripts\\Run/Datasets_P/CountsOutput/Stroma_Cell_Density_by_Patient_Image.txt"
# tumor_output_data_by_image_file_path = "F:\\Script Rebuild\\Transform Essential Scripts\\Run/Datasets_P/CountsOutput/Tumor_Cell_Density_by_Patient_Image.txt";
# cell_image_density_data = CalculateANDWriteImageCellDensitiesinPatientFilePathTree(file_path_tree = file_path_tree, 
#                                                                                    stroma_output_data_by_image_file_path = stroma_output_data_by_image_file_path,
#                                                                                    tumor_output_data_by_image_file_path = tumor_output_data_by_image_file_path)
# cat("done\n")
# 
# # output by patient counts to _tumor_densities_by_patient.txt or _stroma_densities_by_patient.txt
# cat("Counting phenotypes by patient ... \n")
# stroma_data_by_patient_output_file_path = "F:\\Script Rebuild\\Transform Essential Scripts\\Run/Datasets_P/CountsOutput/Stroma_Cell_Density_by_Patient.txt";
# tumor_data_by_patient_output_file_path = "F:\\Script Rebuild\\Transform Essential Scripts\\Run/Datasets_P/CountsOutput/Tumor_Cell_Density_by_Patient.txt";
# Output = CalculateANDWriteCellDensitiesbyPatientinPatientFilePathTree(file_path_tree = file_path_tree, 
#                                                                       cell_image_density_data = cell_image_density_data,
#                                                                       stroma_data_by_patient_output_file_path = stroma_data_by_patient_output_file_path,
#                                                                       tumor_data_by_patient_output_file_path = tumor_data_by_patient_output_file_path)
# cat("done\n")

# compute nearest neighbors
cat("Computing cell's nearest neighbors by image ... \n")
CalculateCellNearestNeighborbyImageinPatientFilePathTree(file_path_tree = file_path_tree)
cat("done")

# combine data by patient
cell_distance_by_patient_data_file_path = "F:\\Script Rebuild\\Transform Essential Scripts\\Run/Datasets_P\\DistanceOutput\\"
cat("Combining cell distance data by patient ... \n")
CombineNearestNeighborsbyPatient(file_path_tree = file_path_tree,cell_distance_by_patient_data_file_path = cell_distance_by_patient_data_file_path)
cat("done\n")

# # combine all data
cat("Combining all distance data ... \n")
bound_data_file_path = "F:\\Script Rebuild\\Transform Essential Scripts\\Run/Datasets_P/DistanceOutput/BoundData.txt"
bound_data = BindAllDataInDir(file_path = cell_distance_by_patient_data_file_path )
cat("done\n")
cat("Writing bound data to [ ",bound_data_file_path, "] \n")
write.table(x = bound_data,file = bound_data_file_path,sep = "\t",row.names = FALSE)
cat("done\n")
cat('---------------------------- END MAIN -----------------------------\n')
################# END main ################# 
