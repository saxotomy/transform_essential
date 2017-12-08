################# Watermark
# Author: Thomas Hart
# Contact: thomas.hart@columbia.edu
################# Watermark

# ########### START updates list
# 08/22/2016 @ 1424
#   1. file created
# 08/22/3016 @ 1457
#   1. file finished
# 09/16/2016 @ 13:42
#   1. modified GetPatientId to only check the tail end for [number]_cell_seg_data.txt
#   2. still requires a " " between the patient id and the image specifier
# ########### END updatse list

# ########### START announcer
cat(
  "ConstructPatientFilePathTree.R contains:\n",
  "\tGetPatientId(file_paths)\n",
  "\tGetPatientImageId(file_paths)\n",
  "\tExtractCellSegDataFilePath(file_names, patient_image_id)\n",
  "\tExtractCellSegDataSummaryFilePath(file_names, patient_image_id)\n",
  "\tExtractTissueSegDataFilePath(file_names, patient_image_id)\n",
  "\tExtractTissueSegDataSummaryFilePath(file_names, patient_image_id)\n",
  "\tExtractScoreDataFilePath(file_names, patient_image_id)\n",
  "\tBuildPatientFilePathTree(dir_name)\n",sep=""
)
# ########### END announcer

# rm(list=ls())

GetPatientId = function(file_paths) {
  temp = file_paths[ grepl(pattern = "_cell_seg_data.txt",x = file_paths) ]
  temp = sub(pattern = "[0-9]+_cell_seg_data.txt",replacement = "",x = temp)
  temp = sub(pattern = ".*[\\]",replacement = "",x = temp)
  temp = sub(pattern = ".*/",replacement = "",x = temp)
  ids = c(); 
  for ( t in temp ) { ids = c(ids,strsplit(x = t,split = " ")[[1]][1]) }
  return( ids )
}

# ### START GetPatientId Check
# file_paths = c(
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/Batch.log",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/batch_procedure.ifp",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_binary_seg_maps.tif",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_cell_seg_data_P.txt",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_cell_seg_data_summary.txt",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_component_data.tif",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_composite_image.tif ",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_image_with_phenotype_map.tif",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 D_1_cell_seg_data.txt",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_tissue_seg_data.txt "
# )
# patient_ids = GetPatientId(file_paths = file_paths)
# check = c("cu-04-06",'cu-05-06')
# ### END GetPatientId Check

GetPatientImageId = function(file_paths) {
  temp = file_paths[ grepl(pattern = "_cell_seg_data.txt",x = file_paths) ]
  temp = sub(pattern = ".*[\\]",replacement = "",x = temp)
  temp = sub(pattern = ".*/",replacement = "",x = temp)
  return( sub(pattern = "_cell_seg_data.txt",replacement = "",x = temp) )
}

# ### START GetPatientImageId Check
# file_paths = c(
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/Batch.log",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/batch_procedure.ifp",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_binary_seg_maps.tif",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_cell_seg_data.txt",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_cell_seg_data_summary.txt",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_component_data.tif",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_composite_image.tif ",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_image_with_phenotype_map.tif",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-08-06 D_1_cell_seg_data.txt",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-08-06C 2_3_cell_seg_data.txt",
#   "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_1_tissue_seg_data.txt"
# )
# patient_image_ids = GetPatientImageId(file_paths = file_paths)
# check = c("cu-04-06 D_1","cu-08-06 D_1","cu-08-06C 2_3")
# ### END GetPatientImageId Check

ExtractCellSegDataFilePath = function(file_names,patient_image_id) {
  return( file_names[ grepl(pattern = paste(patient_image_id,"_cell_seg_data.txt",sep="",collapse = ""),x = file_names,fixed = TRUE)] )
}

# ### START check ExtractCellSegDataFilePath
# file_names = c(
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_cell_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06D A7 S_3_cell_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_cell_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01D A3 D_5_cell_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_cell_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_cell_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_4_tissue_seg_data.txt"
# )
# path = ExtractCellSegDataFilePath(file_names = file_names,patient_image_id = "cu-05-06 A7 S_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_cell_seg_data.txt";
# 
# path = ExtractCellSegDataFilePath(file_names = file_names,patient_image_id = "cu-07-05 B D_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_cell_seg_data.txt";
# 
# path = ExtractCellSegDataFilePath(file_names = file_names,patient_image_id = "cu-12-01 A3 D_5")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_cell_seg_data.txt";
# 
# path = ExtractCellSegDataFilePath(file_names = file_names,patient_image_id = "cu-12-01 A3 D_5")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_cell_seg_data.txt";
# ### END check ExtractCellSegDataFilePath

ExtractCellSegDataSummaryFilePath = function(file_names,patient_image_id) {
  return( file_names[ grepl(pattern = paste(patient_image_id,"_cell_seg_data_summary.txt",sep = "",collapse = ""),x = file_names,fixed = TRUE) ] )
}

# ### START check ExtractCellSegDataFilePath
# file_names = c(
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_cell_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_cell_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_cell_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_cell_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_4_tissue_seg_data.txt"
# )
# path = ExtractCellSegDataSummaryFilePath(file_names = file_names,patient_image_id = "cu-05-06 A7 S_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_cell_seg_data_summary.txt";
# 
# path = ExtractCellSegDataSummaryFilePath(file_names = file_names,patient_image_id = "cu-07-05 B D_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_cell_seg_data_summary.txt";
# 
# path = ExtractCellSegDataSummaryFilePath(file_names = file_names,patient_image_id = "cu-12-01 A3 D_5")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_cell_seg_data_summary.txt";
# ### END check ExtractCellSegDataFilePath

ExtractTissueSegDataFilePath = function(file_names,patient_image_id) {
  return( file_names[ grepl(pattern = paste(patient_image_id,"_tissue_seg_data.txt",sep = "",collapse = ""),x = file_names,fixed = TRUE) ] )
}

# ### START check ExtractCellSegDataFilePath
# file_names = c(
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_tissue_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_tissue_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_tissue_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_tissue_seg_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_4_tissue_seg_data.txt"
# )
# path = ExtractTissueSegDataFilePath(file_names = file_names,patient_image_id = "cu-05-06 A7 S_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_tissue_seg_data.txt";
# 
# path = ExtractTissueSegDataFilePath(file_names = file_names,patient_image_id = "cu-07-05 B D_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_tissue_seg_data.txt";
# 
# path = ExtractTissueSegDataFilePath(file_names = file_names,patient_image_id = "cu-12-01 A3 D_5")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_tissue_seg_data.txt";
# ### END check ExtractCellSegDataFilePath

ExtractTissueSegDataSummaryFilePath = function(file_names,patient_image_id) {
  return( file_names[ grepl(pattern = paste(patient_image_id,"_tissue_seg_data_summary.txt",sep = "",collapse = ""),x = file_names,fixed = TRUE) ] )
}

# ### START check ExtractTissueSegDataSummaryFilePath
# file_names = c(
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_tissue_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_tissue_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_tissue_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_tissue_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_4_tissue_seg_data.txt"
# )
# path = ExtractTissueSegDataSummaryFilePath(file_names = file_names,patient_image_id = "cu-05-06 A7 S_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_tissue_seg_data_summary.txt";
# 
# path = ExtractTissueSegDataSummaryFilePath(file_names = file_names,patient_image_id = "cu-07-05 B D_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_tissue_seg_data_summary.txt";
# 
# path = ExtractTissueSegDataSummaryFilePath(file_names = file_names,patient_image_id = "cu-12-01 A3 D_5")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_tissue_seg_data_summary.txt";
# ### END check ExtractTissueSegDataSummaryFilePath

ExtractScoreDataFilePath = function(file_names,patient_image_id) {
  return( file_names[ grepl(pattern = paste(patient_image_id,"_score_data.txt",sep = "",collapse = ""),x = file_names,fixed = TRUE) ] )
}

# ### START check ExtractScoreDataFilePath
# file_names = c(
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_score_data.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_tissue_seg_data_summary.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-12-01 A3 D_5_score_data.txt.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_score_data.txt.txt",
# "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-04-06 D_4_tissue_seg_data.txt"
# )
# path = ExtractScoreDataFilePath(file_names = file_names,patient_image_id = "cu-05-06 A7 S_3")
# check = c("G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_score_data.txt",
#           "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-05-06 A7 S_3_score_data.txt.txt");
# 
# path = ExtractScoreDataFilePath(file_names = file_names,patient_image_id = "cu-07-05 B D_3")
# check = "G:\\Script Rebuild\\Data//Batch Original 3-22 - 4-11/cu-07-05 B D_3_tissue_seg_data_summary.txt";
# 
# path = ExtractScoreDataFilePath(file_names = file_names,patient_image_id = "cu-12-01 A3 D_5")
# check = character(0);
# ### END check ExtractScoreDataFilePath

ConstructPatientFilePathTree = function(dir_name) {
  # get all file paths
  file_names = list.files(
    path = dir_name,
    recursive = TRUE,
    full.names = TRUE
  )
  
  # cell seg data file paths
  cell_seg_data_file_paths = file_names[ grepl(pattern = "_cell_seg_data.txt",x = file_names,fixed = TRUE) ]
  
  # extract patient ids
  patient_ids = GetPatientId(file_paths = cell_seg_data_file_paths) 
  patient_image_ids = GetPatientImageId(file_paths = cell_seg_data_file_paths)
  
  # build file tree
  file_path_tree = list()
  
  # for each patient id
  # find all associated patient image id
  # for each patient image id
  # place _cell_seg_data under _cell_seg_data
  # place _cell_seg_data_summary under _cell_seg_data_summary
  # place _tissue_seg_data under _tissue_seg_data
  # place _tissue_seg_data_summary under _tissue_seg_data_summary
  # place _score_data under _score_data
  
  # ### START for loop check
  # patient_ids = patient_ids[1]
  # ### END for loop check
  
  for ( patient_id in unique(patient_ids) ) {
    # cat(patient_id,'\n')
    associated_patient_image_ids = patient_image_ids[ grepl(pattern = patient_id,x = patient_image_ids,fixed = TRUE) ]
    for ( patient_image_id in associated_patient_image_ids ) {
      # cat(patient_id,"$",patient_image_id,'\n')
      file_path_tree[[patient_id]][[patient_image_id]]$cell_seg_data = ExtractCellSegDataFilePath(file_names = file_names,patient_image_id = patient_image_id)
      file_path_tree[[patient_id]][[patient_image_id]]$cell_seg_data_summary = ExtractCellSegDataSummaryFilePath(file_names = file_names,patient_image_id = patient_image_id)
      file_path_tree[[patient_id]][[patient_image_id]]$tissue_seg_data = ExtractTissueSegDataFilePath(file_names = file_names,patient_image_id = patient_image_id)
      file_path_tree[[patient_id]][[patient_image_id]]$tissue_seg_data_summary = ExtractTissueSegDataSummaryFilePath(file_names = file_names,patient_image_id = patient_image_id)
      file_path_tree[[patient_id]][[patient_image_id]]$score_data = ExtractScoreDataFilePath(file_names = file_names,patient_image_id = patient_image_id)
    }
  }
  
  # return tree
  return(file_path_tree)
}

# ### START BuildPatientFileTree check
# source_dir = "G:\\Script Rebuild\\Data/"
# file_path_tree = ConstructPatientFilePathTree(dir_name = source_dir)
# output_names = names(file_path_tree)
# check_names = unique( GetPatientId(list.files(path = source_dir,full.names = TRUE,rec=TRUE)) )
# ### END BuildPatientFileTree Check

GetCellSegDataFilePath = function(patient_id,patient_image_id,patient_file_path_tree) {
  return( patient_file_path_tree[[patient_id]][[patient_image_id]]$cell_seg_data )
}

# ### START GetCellSegDataFilePath check
# source_dir = "G:\\Script Rebuild\\Data/"
# file_path_tree = ConstructPatientFilePathTree(dir_name = source_dir)
# GetCellSegDataFilePath(patient_id = 'cu-04-06',patient_image_id = 'cu-04-06 D_1',patient_file_path_tree = file_path_tree)
# GetCellSegDataFilePath(patient_id = 'cu-05-06',patient_image_id = 'cu-05-06 A7 S_1',patient_file_path_tree = file_path_tree)
# ### END GetCellSegDataFilePath check

GetCellSegDataSummaryFilePath = function(patient_id,patient_image_id,patient_file_path_tree) {
  return( patient_file_path_tree[[patient_id]][[patient_image_id]]$cell_seg_data_summary )
}

# ### START GetCellSegDataSummaryFilePath check
# source_dir = "G:\\Script Rebuild\\Data/"
# file_path_tree = ConstructPatientFilePathTree(dir_name = source_dir)
# GetCellSegDataSummaryFilePath(patient_id = 'cu-04-06',patient_image_id = 'cu-04-06 D_1',patient_file_path_tree = file_path_tree)
# GetCellSegDataSummaryFilePath(patient_id = 'cu-05-06',patient_image_id = 'cu-05-06 A7 S_1',patient_file_path_tree = file_path_tree)
# ### END GetCellSegDataSummaryFilePath check


GetTissueSegDataFilePath = function(patient_id,patient_image_id,patient_file_path_tree) {
  return( patient_file_path_tree[[patient_id]][[patient_image_id]]$tissue_seg_data )
}

# ### START GetTissueSegDataFilePath check
# source_dir = "G:\\Script Rebuild\\Data/"
# file_path_tree = ConstructPatientFilePathTree(dir_name = source_dir)
# GetTissueSegDataFilePath(patient_id = 'cu-04-06',patient_image_id = 'cu-04-06 D_1',patient_file_path_tree = file_path_tree)
# GetTissueSegDataFilePath(patient_id = 'cu-05-06',patient_image_id = 'cu-05-06 A7 S_1',patient_file_path_tree = file_path_tree)
# ### END GetTissueSegDataFilePath check

GetTissueSegDataSummaryFilePath = function(patient_id,patient_image_id,patient_file_path_tree) {
  return( patient_file_path_tree[[patient_id]][[patient_image_id]]$tissue_seg_data_summary )
}

# ### START GetTissueSegDataSummaryFilePath check
# source_dir = "G:\\Script Rebuild\\Data/"
# file_path_tree = ConstructPatientFilePathTree(dir_name = source_dir)
# GetTissueSegDataSummaryFilePath(patient_id = 'cu-04-06',patient_image_id = 'cu-04-06 D_1',patient_file_path_tree = file_path_tree)
# GetTissueSegDataSummaryFilePath(patient_id = 'cu-05-06',patient_image_id = 'cu-05-06 A7 S_1',patient_file_path_tree = file_path_tree)
# ### END GetTissueSegDataSummaryFilePath check

GetScoreDataFilePath = function(patient_id,patient_image_id,patient_file_path_tree) {
  return( patient_file_path_tree[[patient_id]][[patient_image_id]]$score_data )
}

# ### START GetScoreDataFilePath check
# source_dir = "G:\\Script Rebuild\\Data/"
# file_path_tree = ConstructPatientFilePathTree(dir_name = source_dir)
# GetScoreDataFilePath(patient_id = 'cu-04-06',patient_image_id = 'cu-04-06 D_1',patient_file_path_tree = file_path_tree)
# GetScoreDataFilePath(patient_id = 'cu-05-06',patient_image_id = 'cu-05-06 A7 S_1',patient_file_path_tree = file_path_tree)
# ### END GetScoreDataFilePath check


