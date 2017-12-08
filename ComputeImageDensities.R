#################
# Author: Thomas Hart
# Contact: thomas.hart@columbia.edu
#################

################# START updates list 
# 08/24/2016 @ 0935
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
  "ComputeImageDensities.R contains:\n",
  "\tCountTransformedPhenotypes( cell_seg_data, tissue_category, phenotypes_to_count  )\n",
  "\tCovertPixelSqtoMMSq(pixels_sq)\n",
  "\tCalculateTissueAreainMM(tissue_seg_data_summary,tissue_category)\n",
  "\tCalculateCellDensitiesbyImage(cell_seg_data, tissue_seg_data_summary, tissue_category, phenotypes_to_count)\n",
  "\tCalculateCellDensitiesbyImageinPatientFilePathTree = function(file_path_tree)\n",
  "\tCalculateANDWriteImageCellDensitiesinPatientFilePathTree = function(file_path_tree,stroma_output_data_file_path,tumor_output_data_file_path)\n"
)
################# END announcer 

# ################# START USER INPUT ################# 
# 
# # totals
# phenotypes_to_count = list(
#   c("CD3+","CD8+"),
#   c("CD3+"),
#   c("CD3-CD8+"),
#   c("Sox10+"),
#   c("Sox10+","Ki67+"),
#   c("Sox10+","Ki67-"),
#   c("CD68+"),
#   c("Sox10+","CD68+"),
#   c("CD3+","CD8+","Ki67+"),
#   c("CD3+","CD8-","HLA-DR+"),
#   c("CD3+","CD8-","CD68+","HLA-DR+"),
#   # HLA-DR+
#   c("CD3+","CD8+","HLA-DR+"),
#   c("CD3+","HLA-DR+"),
#   c("CD3-CD8+","HLA-DR+"),
#   c("Sox10+","HLA-DR+"),
#   c("Sox10+","Ki67+","HLA-DR+"),
#   c("Sox10+","Ki67-","HLA-DR+"),
#   c("CD68+","HLA-DR+"),
#   c("Sox10+","CD68+","HLA-DR+"),
#   # HLA-DR-
#   c("CD3+","CD8+","HLA-DR-"),
#   c("CD3+","HLA-DR-"),
#   c("CD3-CD8+","HLA-DR-"),
#   c("Sox10+","HLA-DR-"),
#   c("Sox10+","Ki67+","HLA-DR-"),
#   c("Sox10+","Ki67-","HLA-DR-"),
#   c("CD68+","HLA-DR-"),
#   c("Sox10+","CD68+","HLA-DR-")
# )
# 
# ################# END USER INPUT ################# 

CountTransformedPhenotypes = function( cell_seg_data, tissue_category, phenotypes_to_count  ) {
  Phenotypes = as.character(cell_seg_data$Transformed_Phenotype)
  Phenotypes = Phenotypes[ as.character(cell_seg_data$Tissue_Category) == tissue_category ]
  counts = list()
  for ( i in 1:length(phenotypes_to_count) ) {
    logical = matrix(data = logical(),nrow = length(Phenotypes),ncol = length(phenotypes_to_count[[i]]))
    for ( j in 1:length(phenotypes_to_count[[i]]) ) {
      if ( j==1 ) {
        temp = sub(pattern = "[+]",replacement = "[+]",x = phenotypes_to_count[[i]][j]);
        temp = sub(pattern = "[-]",replacement = "[-]",x = temp);
        logical[,j] = grepl(pattern = paste("^",temp,sep = ""),x = Phenotypes,fixed = FALSE)
      } else {
        logical[,j] = grepl(pattern = phenotypes_to_count[[i]][j],x = Phenotypes,fixed = TRUE) 
      }
    }
    counts[[ paste(phenotypes_to_count[[i]],collapse = "") ]] = sum( apply(X = logical,MARGIN = 1,FUN = all) )
  }
  return( counts )
}

# ### START CountTransformedPhenotypes check
# # read cell_seg_data
# cell_seg_data = read.csv(file = "F:\\Script Rebuild\\Training_Data_P/Batch Original 3-22 - 4-11/cu-04-06 D_1_cell_seg_data.txt",sep = "\t")
# # calculate counts
# tumor_counts = CountTransformedPhenotypes(cell_seg_data = cell_seg_data,tissue_category = "Tumor",phenotypes_to_count = phenotypes_to_count)
# tumor_counts_check = list(
#   "CD3+_CD8+" = 7,
#   "CD3+" = 75,
#   "CD3-CD8+" = 21,
#   "Sox10+" = 99,
#   "Sox10+_Ki67+" = 1,
#   "Sox10+_Ki67-" = 98
# )
# 
# # calculate counts
# stroma_counts = CountTransformedPhenotypes(cell_seg_data = cell_seg_data,tissue_category = "Stroma",phenotypes_to_count = phenotypes_to_count)
# stroma_counts_check = list(
#   "CD3+_CD8+" = 52,
#   "CD3+" = 376,
#   "CD3-CD8+" = 80,
#   "Sox10+" = 13,
#   "Sox10+_Ki67+" = 0,
#   "Sox10+_Ki67-" = 13
# )
# ### END CountTransformedPhenotypes check

CovertPixelSqtoMMSq = function(pixels_sq) { return( pixels_sq * 0.00246 ) } # approved by Ed 2016/09/06 17:18

CalculateTissueAreainMM = function(tissue_seg_data_summary,tissue_category) {
  tissue_area_pixels = tissue_seg_data_summary$Region_Area_pixels[ which( as.character(tissue_seg_data_summary$Tissue_Category) == tissue_category) ]
  tissue_area_mm = CovertPixelSqtoMMSq(pixels_sq = tissue_area_pixels) 
  return( list("tissue_area_mm"=tissue_area_mm,
               "tissue_area_pixels"=tissue_area_pixels))
}

# ### START CalculateTissueAreainMM check
# # read cell_seg_data
# tissue_seg_data_summary = read.csv(file = "F:\\Script Rebuild\\Training_Data_P/Batch Original 3-22 - 4-11/cu-04-06 D_1_tissue_seg_data_summary.txt",sep = "\t")
# # tissue category
# tissue_category = "Tumor"
# # calculate tissue arrea
# CalculateTissueAreainMM(tissue_seg_data_summary = tissue_seg_data_summary,tissue_category = tissue_category)
# Check = 7951.373;
# 
# # read cell_seg_data
# tissue_seg_data_summary = read.csv(file = "F:\\Script Rebuild\\Training_Data_P/Batch Original 3-22 - 4-11/cu-04-06 D_1_tissue_seg_data_summary.txt",sep = "\t")
# # tissue category
# tissue_category = "Stroma"
# # calculate tissue arrea
# CalculateTissueAreainMM(tissue_seg_data_summary = tissue_seg_data_summary,tissue_category = tissue_category)
# Check = 93392.51;
# ### END CalculateTissueAreainMM check

# for each patietn
  # for each image 
    # read in data
    # for tissue in {"Tumor", "Stroma"}
      # CountTransformedPhenotypes( cell_seg_data, tissue_category, transformed_phenotypes  )

CalculateCellDensitiesbyImage = function(cell_seg_data, tissue_seg_data_summary, tissue_category, phenotypes_to_count) {
  counts = CountTransformedPhenotypes(cell_seg_data = cell_seg_data,tissue_category = tissue_category,phenotypes_to_count = phenotypes_to_count)
  count_names = names(counts)
  area = CalculateTissueAreainMM(tissue_seg_data_summary = tissue_seg_data_summary,tissue_category = tissue_category)
  data_list = list(); 
  data_list[["Area_mm"]] = area$tissue_area_mm
  data_list[["Area_px"]] = area$tissue_area_pixels
  for ( i in 1:length(counts) ) {
    data_list[[ paste("Count",paste(count_names[i],collapse=""),sep="_") ]] = counts[[i]]
    data_list[[ paste("Density",paste(count_names[i],collapse=""),sep = "_")]] = counts[[i]]/area$tissue_area_mm
  } 
  data_frame = data.frame(data_list); names(data_frame) = names(data_list);
  return(data_frame)
}

# ### START CalculateCellDensitiesbyImage check
# # read cell_seg_data
# cell_seg_data = read.csv(file = "G:\\Script Rebuild\\Training_Data_P/Batch Original 3-22 - 4-11/cu-04-06 D_1_cell_seg_data.txt",sep = "\t")
# # read cell_seg_data
# tissue_seg_data_summary = read.csv(file = "G:\\Script Rebuild\\Training_Data_P/Batch Original 3-22 - 4-11/cu-04-06 D_1_tissue_seg_data_summary.txt",sep = "\t")
# # tissue category
# tissue_category = "Tumor"
# # phenotypes to count
# # see above
# Output = CalculateCellDensitiesbyImage(cell_seg_data = cell_seg_data,tissue_seg_data_summary = tissue_seg_data_summary,
#                             tissue_category = tissue_category, phenotypes_to_count = phenotypes_to_count)
# ### END CalculateCellDensitiesbyImage check

CalculateCellDensitiesbyImageinPatientFilePathTree = function(file_path_tree) {
  # read in, correct and rewrite each file
  tumor_df_list = list()
  stroma_df_list = list()
  for ( patient_id in names(file_path_tree) ) {
    cat("\tpatient_id: ",patient_id,"\n")
    for ( patient_image_id in names(file_path_tree[[patient_id]]) ) {
      cat("\t\timage_id: ",patient_image_id,"\n")
      cat("\t\t\tcomputing counts ... ")
      # read cell seg data
      cell_seg_data_file_path = GetCellSegDataFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      cell_seg_data = read.csv(file = cell_seg_data_file_path,sep = "\t")
      
      # read tissue seg data summary
      tissue_seg_data_summary_file_path = GetTissueSegDataSummaryFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      tissue_seg_data_summary = read.csv(file = tissue_seg_data_summary_file_path,sep = "\t")
      
      # compute
      patient_frame = data.frame(Patient_Id = patient_id,Patient_Image_Id = patient_image_id)
      tumor_density_frame = CalculateCellDensitiesbyImage(cell_seg_data = cell_seg_data,tissue_seg_data_summary = tissue_seg_data_summary,tissue_category = "Tumor",phenotypes_to_count = phenotypes_to_count)
      stroma_density_frame = CalculateCellDensitiesbyImage(cell_seg_data = cell_seg_data,tissue_seg_data_summary = tissue_seg_data_summary,tissue_category = "Stroma",phenotypes_to_count = phenotypes_to_count)
    
      # add to df list
      tumor_df = cbind(patient_frame,tumor_density_frame); names(tumor_df) = c(names(patient_frame),names(tumor_density_frame))
      tumor_df_list[[ length(tumor_df_list)+1 ]] = tumor_df
      stroma_df = cbind(patient_frame,stroma_density_frame); names(stroma_df) = c(names(patient_frame),names(stroma_density_frame))
      stroma_df_list[[ length(stroma_df_list)+1 ]] = stroma_df
      cat("done\n")
    }
  }
  tumor_df = tumor_df_list[[1]]
  for ( i in 2:length(tumor_df_list) ) { tumor_df = rbind(tumor_df,tumor_df_list[[i]]) }
  stroma_df = stroma_df_list[[1]]
  for ( i in 2:length(stroma_df_list) ) { stroma_df = rbind(stroma_df,stroma_df_list[[i]]) }
  return( list( "Tumor"=tumor_df, "Stroma"=stroma_df ) )
}

# ### START CalculateCellDensitiesbyImageinPatientFilePathTree check
# # construct file path tree
# file_path_tree = ConstructPatientFilePathTree(dir_name = "F:\\Script Rebuild\\Training_Data_P\\")
# # calculate densities
# Output = CalculateCellDensitiesbyImageinPatientFilePathTree(file_path_tree = file_path_tree)
# ### END CalculateCellDensitiesbyImageinPatientFilePathTree check

CalculateANDWriteImageCellDensitiesinPatientFilePathTree = function(file_path_tree,
                                                                    stroma_output_data_by_image_file_path,
                                                                    tumor_output_data_by_image_file_path) {
  # calculate densities
  Output = CalculateCellDensitiesbyImageinPatientFilePathTree(file_path_tree = file_path_tree)
  # write outputs
  cat("writing : stroma count data -> [ ",stroma_output_data_by_image_file_path,"] ... ")
  write.table(x = Output$Stroma,file = stroma_output_data_by_image_file_path,sep = "\t",row.names = FALSE)
  cat("done\n")
  cat("writing : tumor count data  -> [ ",tumor_output_data_by_image_file_path,"] ... ")
  write.table(x = Output$Tumor,file = tumor_output_data_by_image_file_path,sep = "\t",row.names = FALSE)
  cat("done\n")
  return(Output)
}

# ### START CalculateANDWriteImageCellDensitiesinPatientFilePathTree check
# stroma_output_data_by_image_file_path = "F:\\Script Rebuild\\Training_Data_P\\/CountsOutput/Stroma_Cell_Density_by_Patient_Image.txt";
# tumor_output_data_by_image_file_path = "F:\\Script Rebuild\\Training_Data_P\\/CountsOutput/Tumor_Cell_Density_by_Patient_Image.txt";
# # file path tree
# file_path_tree = ConstructPatientFilePathTree(dir_name = "F:\\Script Rebuild\\Training_Data_P\\")
# #
# CalculateANDWriteImageCellDensitiesinPatientFilePathTree(file_path_tree = file_path_tree,
#                                                          stroma_output_data_by_image_file_path = stroma_output_data_by_image_file_path,
#                                                          tumor_output_data_by_image_file_path = tumor_output_data_by_image_file_path)
# ### END CalculateANDWriteImageCellDensitiesinPatientFilePathTree check

CalculateCellDensitiesbyPatient = function(cell_image_density_data, patient_id, tissue_category) {
  tissue_data = cell_image_density_data[[tissue_category]]
  patient_tissue_data = tissue_data[ as.character(tissue_data$Patient_Id) == patient_id , ]
  
  # get count
  tissue_total_area_px = sum(patient_tissue_data$Area_px)
  tissue_total_area_mm = CovertPixelSqtoMMSq(pixels_sq = tissue_total_area_px)
  tissue_total_counts_list = lapply(X = patient_tissue_data[ , grepl(pattern = "Count_",x = names(patient_tissue_data),fixed = TRUE) ],FUN = sum)
  total_counts_names = names(tissue_total_counts_list);
  
  # create data list
  data_list = list()
  data_list[["Patient_Id"]] = patient_id;
  data_list[["Total_Area_mm"]] = tissue_total_area_mm;
  data_list[["Total_Area_px"]] = tissue_total_area_px;
  
  for ( i in 1:length(tissue_total_counts_list) ) {
    data_list[[ paste("Total",total_counts_names[i],sep="_") ]] = tissue_total_counts_list[[i]]
    temp = sub(pattern = "Count",replacement = "Density",x = total_counts_names[i],fixed = TRUE)
    data_list[[ temp ]] = tissue_total_counts_list[[i]]/tissue_total_area_mm
  }
  data_frame = data.frame(data_list); names(data_frame) = names(data_list);
  return(data_frame)
}

# ### START CalculateCellDensitiesbyPatient check
# file_path_tree = ConstructPatientFilePathTree(dir_name = "F:\\Script Rebuild\\Training_Data_P\\");
# cell_image_density_data = CalculateCellDensitiesbyImageinPatientFilePathTree(file_path_tree = file_path_tree);
# patient_tumor_density_data = CalculateCellDensitiesbyPatient(cell_image_density_data = cell_image_density_data,
#                                                        patient_id = "cu-04-06",tissue_category = "Tumor")
# patient_stroma_density_data = CalculateCellDensitiesbyPatient(cell_image_density_data = cell_image_density_data,
#                                                              patient_id = "cu-04-06",tissue_category = "Stroma")
# ### END CalculateCellDensitiesbyPatient check

CalculateCellDensitiesbyPatientinPatientFilePathTree = function(file_path_tree,
                                                                cell_image_density_data) {
  total_tumor_df_list = list()
  total_stroma_df_list = list()
  for ( patient_id in names(file_path_tree) ) {
    cat("\tpatient_id: ",patient_id,"\n")
    cat("\t\tCalculating ... ")
    total_tumor_density_df = CalculateCellDensitiesbyPatient(cell_image_density_data = cell_image_density_data,
                                                             patient_id = patient_id,tissue_category = "Tumor")
    total_tumor_df_list[[ length(total_tumor_df_list)+1 ]] = total_tumor_density_df

    total_stroma_density_df = CalculateCellDensitiesbyPatient(cell_image_density_data = cell_image_density_data,
                                                              patient_id = patient_id,tissue_category = "Stroma")
    total_stroma_df_list[[ length(total_stroma_df_list)+1 ]] = total_stroma_density_df
    cat("done\n")
  }
  tumor_df = total_tumor_df_list[[1]]
  for ( i in 2:length(total_tumor_df_list) ) { tumor_df = rbind(tumor_df,total_tumor_df_list[[i]]) }
  stroma_df = total_stroma_df_list[[1]]
  for ( i in 2:length(total_stroma_df_list) ) { stroma_df = rbind(stroma_df,total_stroma_df_list[[i]]) }
  return( list( "Tumor"=tumor_df, "Stroma"=stroma_df ) )
}

# ### START CalculateCellDensitiesbyImage check
# file_path_tree = ConstructPatientFilePathTree(dir_name = "F:\\Script Rebuild\\Training_Data_P\\");
# cell_image_density_data = CalculateCellDensitiesbyImageinPatientFilePathTree(file_path_tree = file_path_tree);
# cell_patient_density_data = CalculateCellDensitiesbyPatientinPatientFilePathTree(file_path_tree = file_path_tree,
#                                                                                  cell_image_density_data = cell_image_density_data)
# ### END CalculateCellDensitiesbyImage check

CalculateANDWriteCellDensitiesbyPatientinPatientFilePathTree = function(file_path_tree,
                                                                cell_image_density_data,
                                                                stroma_data_by_patient_output_file_path,
                                                                tumor_data_by_patient_output_file_path) {
  Output = CalculateCellDensitiesbyPatientinPatientFilePathTree(file_path_tree = file_path_tree,
                                                                                   cell_image_density_data = cell_image_density_data);
  # write outputs
  cat("writing : stroma count data -> [ ",stroma_data_by_patient_output_file_path,"] ... ")
  write.table(x = Output$Stroma,file = stroma_data_by_patient_output_file_path,sep = "\t",row.names = FALSE)
  cat("done\n")
  cat("writing : tumor count data  -> [ ",tumor_data_by_patient_output_file_path,"] ... ")
  write.table(x = Output$Tumor,file = tumor_data_by_patient_output_file_path,sep = "\t",row.names = FALSE)
  cat("done\n")
  return( Output )
}

# ### START CalculateANDWriteCellDensitiesbyPatientinPatientFilePathTree check
# file_path_tree = ConstructPatientFilePathTree(dir_name = "F:\\Script Rebuild\\Training_Data_P\\");
# cell_image_density_data = CalculateCellDensitiesbyImageinPatientFilePathTree(file_path_tree = file_path_tree);
# stroma_data_by_patient_output_file_path = "F:\\Script Rebuild\\Training_Data_P\\/CountsOutput/Stroma_by_Patient_Cell_Density.txt";
# tumor_data_by_patient_output_file_path = "F:\\Script Rebuild\\Training_Data_P\\/CountsOutput/Tumor_by_Patient_Cell_Density.txt";
# cell_patient_density_data = CalculateANDWriteCellDensitiesbyPatientinPatientFilePathTree(file_path_tree = file_path_tree,
#                                                                                  cell_image_density_data = cell_image_density_data,
#                                                                                  stroma_data_by_patient_output_file_path = stroma_data_by_patient_output_file_path,
#                                                                                  tumor_data_by_patient_output_file_path = tumor_data_by_patient_output_file_path);
# ### END CalculateANDWriteCellDensitiesbyPatientinPatientFilePathTree check
