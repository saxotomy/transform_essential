#################
# Author: Thomas Hart
# Contact: thomas.hart@columbia.edu
#################

################# START updates list 
# 08/30/2016 @ 0956
#   1. file created
# 09/19/2016 @ 1230
#   1. CombineNearestNeighborsbyPatient modified to process patients with only one image
# 10/13/2016 @ 15:49
#   1. CombineNearestNeighborsbyPatient modified output x and y coordinates of nearest neighbors as well
################# END updatse list 

# rm(list=ls())

# ################# START source 
# 
# setwd("G:\\Script Rebuild\\Transform Script\\")
# source(file = "ComputeImageDensities.R")
# source(file = "ConstructPatientFilePathTree.R")
# 
# ################# END source 

################# START announcer 
cat(
  "ComputeNearestNeighbors.R contains:\n",
  "Distance(x1,y1,x2,y2)\n",
  "GetCellCoordinatesByPhenotype( cell_seg_data, phenotypes_to_count  )\n",
  "GetNearestNeighborByPhenotype( x, y, cell_coordinates_by_phenotype, exclude_row = 0 )\n",
  "CalculateCellNearestNeighborbyImageinPatientFilePathTree( file_path_tree )\n",
  "CombineNearestNeighborsbyPatient( file_path_tree , cell_distance_by_patient_data_file_path)\n"
)
################# END announcer 

# ################# START user input
# phenotypes_to_count = list(
#   c("CD3+"),
#   c("CD3+","CD68-"),
#   c("CD3+","CD8-"),
#   c("CD3+","HLA-DR-"),
#   c("CD3+","Ki67-"),
#   c("CD3+","CD68+"),
#   c("CD3+","CD8+"),
#   c("CD3+","HLA-DR+"),
#   c("CD3+","Ki67+"),
#   c("CD3+","CD68-","CD8-"),
#   c("CD3+","CD68-","HLA-DR-"),
#   c("CD3+","CD68-","Ki67-"),
#   c("CD3+","CD8-","HLA-DR-"),
#   c("CD3+","CD8-","Ki67-"),
#   c("CD3+","HLA-DR-","Ki67-"),
#   c("CD3+","CD68-","CD8+"),
#   c("CD3+","CD68-","HLA-DR+"),
#   c("CD3+","CD68-","Ki67+"),
#   c("CD3+","CD8-","HLA-DR+"),
#   c("CD3+","CD8-","Ki67+"),
#   c("CD3+","HLA-DR-","Ki67+"),
#   c("CD3+","CD68+","CD8-"),
#   c("CD3+","CD68+","HLA-DR-"),
#   c("CD3+","CD68+","Ki67-"),
#   c("CD3+","CD8+","HLA-DR-"),
#   c("CD3+","CD8+","Ki67-"),
#   c("CD3+","HLA-DR+","Ki67-"),
#   c("CD3+","CD68+","CD8+"),
#   c("CD3+","CD68+","HLA-DR+"),
#   c("CD3+","CD68+","Ki67+"),
#   c("CD3+","CD8+","HLA-DR+"),
#   c("CD3+","CD8+","Ki67+"),
#   c("CD3+","HLA-DR+","Ki67+"),
#   c("CD3+","CD68-","CD8-","HLA-DR-"),
#   c("CD3+","CD68-","CD8-","Ki67-"),
#   c("CD3+","CD68-","HLA-DR-","Ki67-"),
#   c("CD3+","CD8-","HLA-DR-","Ki67-"),
#   c("CD3+","CD68-","CD8-","HLA-DR+"),
#   c("CD3+","CD68-","CD8-","Ki67+"),
#   c("CD3+","CD68-","HLA-DR-","Ki67+"),
#   c("CD3+","CD8-","HLA-DR-","Ki67+"),
#   c("CD3+","CD68-","CD8+","HLA-DR-"),
#   c("CD3+","CD68-","CD8+","Ki67-"),
#   c("CD3+","CD68-","HLA-DR+","Ki67-"),
#   c("CD3+","CD8-","HLA-DR+","Ki67-"),
#   c("CD3+","CD68-","CD8+","HLA-DR+"),
#   c("CD3+","CD68-","CD8+","Ki67+"),
#   c("CD3+","CD68-","HLA-DR+","Ki67+"),
#   c("CD3+","CD8-","HLA-DR+","Ki67+"),
#   c("CD3+","CD68+","CD8-","HLA-DR-"),
#   c("CD3+","CD68+","CD8-","Ki67-"),
#   c("CD3+","CD68+","HLA-DR-","Ki67-"),
#   c("CD3+","CD8+","HLA-DR-","Ki67-"),
#   c("CD3+","CD68+","CD8-","HLA-DR+"),
#   c("CD3+","CD68+","CD8-","Ki67+"),
#   c("CD3+","CD68+","HLA-DR-","Ki67+"),
#   c("CD3+","CD8+","HLA-DR-","Ki67+"),
#   c("CD3+","CD68+","CD8+","HLA-DR-"),
#   c("CD3+","CD68+","CD8+","Ki67-"),
#   c("CD3+","CD68+","HLA-DR+","Ki67-"),
#   c("CD3+","CD8+","HLA-DR+","Ki67-"),
#   c("CD3+","CD68+","CD8+","HLA-DR+"),
#   c("CD3+","CD68+","CD8+","Ki67+"),
#   c("CD3+","CD68+","HLA-DR+","Ki67+"),
#   c("CD3+","CD8+","HLA-DR+","Ki67+"),
#   c("CD3+","CD68-","CD8-","HLA-DR-","Ki67-"),
#   c("CD3+","CD68-","CD8-","HLA-DR-","Ki67+"),
#   c("CD3+","CD68-","CD8-","HLA-DR+","Ki67-"),
#   c("CD3+","CD68-","CD8-","HLA-DR+","Ki67+"),
#   c("CD3+","CD68-","CD8+","HLA-DR-","Ki67-"),
#   c("CD3+","CD68-","CD8+","HLA-DR-","Ki67+"),
#   c("CD3+","CD68-","CD8+","HLA-DR+","Ki67-"),
#   c("CD3+","CD68-","CD8+","HLA-DR+","Ki67+"),
#   c("CD3+","CD68+","CD8-","HLA-DR-","Ki67-"),
#   c("CD3+","CD68+","CD8-","HLA-DR-","Ki67+"),
#   c("CD3+","CD68+","CD8-","HLA-DR+","Ki67-"),
#   c("CD3+","CD68+","CD8-","HLA-DR+","Ki67+"),
#   c("CD3+","CD68+","CD8+","HLA-DR-","Ki67-"),
#   c("CD3+","CD68+","CD8+","HLA-DR-","Ki67+"),
#   c("CD3+","CD68+","CD8+","HLA-DR+","Ki67-"),
#   c("CD3+","CD68+","CD8+","HLA-DR+","Ki67+"),
#   c("CD3-CD8+"),
#   c("CD3-CD8+","CD68-"),
#   c("CD3-CD8+","HLA-DR-"),
#   c("CD3-CD8+","Ki67-"),
#   c("CD3-CD8+","CD68+"),
#   c("CD3-CD8+","HLA-DR+"),
#   c("CD3-CD8+","Ki67+"),
#   c("CD3-CD8+","CD68-","HLA-DR-"),
#   c("CD3-CD8+","CD68-","Ki67-"),
#   c("CD3-CD8+","HLA-DR-","Ki67-"),
#   c("CD3-CD8+","CD68-","HLA-DR+"),
#   c("CD3-CD8+","CD68-","Ki67+"),
#   c("CD3-CD8+","HLA-DR-","Ki67+"),
#   c("CD3-CD8+","CD68+","HLA-DR-"),
#   c("CD3-CD8+","CD68+","Ki67-"),
#   c("CD3-CD8+","HLA-DR+","Ki67-"),
#   c("CD3-CD8+","CD68+","HLA-DR+"),
#   c("CD3-CD8+","CD68+","Ki67+"),
#   c("CD3-CD8+","HLA-DR+","Ki67+"),
#   c("CD3-CD8+","CD68-","HLA-DR-","Ki67-"),
#   c("CD3-CD8+","CD68-","HLA-DR-","Ki67+"),
#   c("CD3-CD8+","CD68-","HLA-DR+","Ki67-"),
#   c("CD3-CD8+","CD68-","HLA-DR+","Ki67+"),
#   c("CD3-CD8+","CD68+","HLA-DR-","Ki67-"),
#   c("CD3-CD8+","CD68+","HLA-DR-","Ki67+"),
#   c("CD3-CD8+","CD68+","HLA-DR+","Ki67-"),
#   c("CD3-CD8+","CD68+","HLA-DR+","Ki67+"),
#   c("Sox10+"),
#   c("Sox10+","CD68-"),
#   c("Sox10+","HLA-DR-"),
#   c("Sox10+","Ki67-"),
#   c("Sox10+","CD68+"),
#   c("Sox10+","HLA-DR+"),
#   c("Sox10+","Ki67+"),
#   c("Sox10+","CD68-","HLA-DR-"),
#   c("Sox10+","CD68-","Ki67-"),
#   c("Sox10+","HLA-DR-","Ki67-"),
#   c("Sox10+","CD68-","HLA-DR+"),
#   c("Sox10+","CD68-","Ki67+"),
#   c("Sox10+","HLA-DR-","Ki67+"),
#   c("Sox10+","CD68+","HLA-DR-"),
#   c("Sox10+","CD68+","Ki67-"),
#   c("Sox10+","HLA-DR+","Ki67-"),
#   c("Sox10+","CD68+","HLA-DR+"),
#   c("Sox10+","CD68+","Ki67+"),
#   c("Sox10+","HLA-DR+","Ki67+"),
#   c("Sox10+","CD68-","HLA-DR-","Ki67-"),
#   c("Sox10+","CD68-","HLA-DR-","Ki67+"),
#   c("Sox10+","CD68-","HLA-DR+","Ki67-"),
#   c("Sox10+","CD68-","HLA-DR+","Ki67+"),
#   c("Sox10+","CD68+","HLA-DR-","Ki67-"),
#   c("Sox10+","CD68+","HLA-DR-","Ki67+"),
#   c("Sox10+","CD68+","HLA-DR+","Ki67-"),
#   c("Sox10+","CD68+","HLA-DR+","Ki67+"),
#   c("CD68+"),
#   c("CD68+","CD8-"),
#   c("CD68+","HLA-DR-"),
#   c("CD68+","Ki67-"),
#   c("CD68+","CD8+"),
#   c("CD68+","HLA-DR+"),
#   c("CD68+","Ki67+"),
#   c("CD68+","CD8-","HLA-DR-"),
#   c("CD68+","CD8-","Ki67-"),
#   c("CD68+","HLA-DR-","Ki67-"),
#   c("CD68+","CD8-","HLA-DR+"),
#   c("CD68+","CD8-","Ki67+"),
#   c("CD68+","HLA-DR-","Ki67+"),
#   c("CD68+","CD8+","HLA-DR-"),
#   c("CD68+","CD8+","Ki67-"),
#   c("CD68+","HLA-DR+","Ki67-"),
#   c("CD68+","CD8+","HLA-DR+"),
#   c("CD68+","CD8+","Ki67+"),
#   c("CD68+","HLA-DR+","Ki67+"),
#   c("CD68+","CD8-","HLA-DR-","Ki67-"),
#   c("CD68+","CD8-","HLA-DR-","Ki67+"),
#   c("CD68+","CD8-","HLA-DR+","Ki67-"),
#   c("CD68+","CD8-","HLA-DR+","Ki67+"),
#   c("CD68+","CD8+","HLA-DR-","Ki67-"),
#   c("CD68+","CD8+","HLA-DR-","Ki67+"),
#   c("CD68+","CD8+","HLA-DR+","Ki67-"),
#   c("CD68+","CD8+","HLA-DR+","Ki67+"),
#   c("Other"),
#   c("Other","HLA-DR-"),
#   c("Other","HLA-DR+")
# )
# ################# END user input

Distance = function(x1,y1,x2,y2) { sqrt( (y2-y1)^2 + (x2-x1)^2 )  }

# ### START Distance Check
# set.seed(123)
# X <- matrix(rnorm(10),nrow=5)
# check_mat = dist(x = X,diag = T,upper = T)
# test_mat = matrix(data = NA,nrow = 5,ncol = 5)
# for ( i in 1:nrow(X) ) {
#   for ( j in 1:nrow(X) ) {
#     test_mat[i,j] = Distance(x1 = X[i,1],y1 = X[i,2],x2 = X[j,1],y2 = X[j,2])
#   }
# }
# round(x = test_mat,digits = 7)==round(x = as.matrix(check_mat),digits = 7)
# ### END Distance Check

# # read cell seg data
# cell_seg_data = read.csv(file = "G:\\Script Rebuild\\2016_08_25_12_12_Run\\Datasets_P\\Data\\Batch 8-18\\cu-04-06 D_1_cell_seg_data.txt",sep = "\t")

GetCellCoordinatesByPhenotype = function( cell_seg_data, phenotypes_to_count  ) {
  Phenotypes = as.character(cell_seg_data$Transformed_Phenotype)
  coordinates = list()
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
    selection = apply(X = logical,MARGIN = 1,FUN = all);
    coordinates[[ paste(phenotypes_to_count[[i]],collapse = "") ]] = cell_seg_data[ selection , c("Cell_X_Position","Cell_Y_Position")]
  }
  return( coordinates )
}

# ### START GetCellCoordinatesWithPhenotype check
# cell_seg_data = read.csv(file = "F:\\Script Rebuild\\2016_09_01_12_12_Run\\Datasets_P\\Data\\Batch 8-18\\cu-04-06 D_1_cell_seg_data.txt",sep = "\t")
# CellCoordinateByPhenotypes = GetCellCoordinatesByPhenotype(cell_seg_data = cell_seg_data,phenotypes_to_count = phenotypes_to_count)
# # check
# CountsofCoordinatesByPhenotypes = lapply(X = CellCoordinateByPhenotypes,FUN = nrow)
# CountsofTransformedPhenotypes = list(
#   "Stroma" = CountTransformedPhenotypes(cell_seg_data = cell_seg_data,tissue_category = "Stroma",phenotypes_to_count = phenotypes_to_count),
#   "Tumor" = CountTransformedPhenotypes(cell_seg_data = cell_seg_data,tissue_category = "Tumor",phenotypes_to_count = phenotypes_to_count)
#               )
# Names = names(CountsofTransformedPhenotypes$Stroma);
# CountsCheck = list()
# for ( name in Names ) {
#   CountsCheck[[name]] = CountsofTransformedPhenotypes$Stroma[[name]] + CountsofTransformedPhenotypes$Tumor[[name]]
# }
# ### END GetCellCoordinatesWithPhenotype check

GetNearestNeighborByPhenotype = function( x, y, cell_coordinates_by_phenotype, exclude_row = 0 ) {
  min_dist_by_phenotype = list()
  for ( name in names(cell_coordinates_by_phenotype) ) {
    if ( nrow(cell_coordinates_by_phenotype[[name]]) == 0 ) {
      min_dist_by_phenotype[[name]] = c(NA);
      # min_dist_by_phenotype[[paste(name,"_X_coord",sep="")]] = c(NA);
      # min_dist_by_phenotype[[paste(name,"_Y_coord",sep="")]] = c(NA);
    } else if ( exclude_row != 0 ) {
      row_vals = as.integer(as.character(row.names(cell_coordinates_by_phenotype[[name]])))
      x2 = cell_coordinates_by_phenotype[[name]]$Cell_X_Position[row_vals!=exclude_row]
      y2 = cell_coordinates_by_phenotype[[name]]$Cell_Y_Position[row_vals!=exclude_row]
      if ( length(x2) == 0 || length(y2) == 0 ) {
        min_dist_by_phenotype[[name]] = c(NA);
        # min_dist_by_phenotype[[paste(name,"_X_coord",sep="")]] = c(NA);
        # min_dist_by_phenotype[[paste(name,"_Y_coord",sep="")]] = c(NA);
      } else {
        distances = Distance( x1 = x,y1 = y, x2 = x2, y2 = y2 )
        min_distance = min(distances)
        min_dist_by_phenotype[[name]] = min_distance
        # coordinates = cell_coordinates_by_phenotype[[name]][which(distances==min_distance),]
        # min_dist_by_phenotype[[paste(name,"_X_coord",sep="")]] = as.numeric(coordinates[1]);
        # min_dist_by_phenotype[[paste(name,"_Y_coord",sep="")]] = as.numeric(coordinates[2]);
      }
    } else {
      distances = Distance( x1 = x,y1 = y, x2 = cell_coordinates_by_phenotype[[name]]$Cell_X_Position, y2 = cell_coordinates_by_phenotype[[name]]$Cell_Y_Position) 
      min_distance = min(distances)
      min_dist_by_phenotype[[name]] = min(distances)
      # coordinates = cell_coordinates_by_phenotype[[name]][which(distances==min_distance),]
      # min_dist_by_phenotype[[paste(name,"_X_coord",sep="")]] = as.numeric(coordinates[1]);
      # min_dist_by_phenotype[[paste(name,"_Y_coord",sep="")]] = as.numeric(coordinates[2]);
    }
  }
  return( min_dist_by_phenotype )
}

# ### START GetNearestNeighborByPhenotype check
# cell_seg_data = read.csv(file = "F:\\Script Rebuild\\2016_09_01_12_12_Run\\Datasets_P\\Data\\Batch 8-18\\cu-04-06 D_1_cell_seg_data.txt",sep = "\t")
# cell_coordinates_by_phenotype = GetCellCoordinatesByPhenotype(cell_seg_data = cell_seg_data,phenotypes_to_count = phenotypes_to_count)
# row = 6;
# nearest_neighbor_by_phenotype = GetNearestNeighborByPhenotype(x = cell_seg_data[row,"Cell_X_Position"],y = cell_seg_data[row,"Cell_Y_Position"],cell_coordinates_by_phenotype = cell_coordinates_by_phenotype,exclude_row = row)
# # check
# phenotype = "CD3+HLA-DR-"
# row_vals = as.integer(as.character(row.names(cell_coordinates_by_phenotype[[phenotype]])))
# x1 = cell_seg_data[row,"Cell_X_Position"]; y1 = cell_seg_data[row,"Cell_Y_Position"];
# min( Distance(x1 = x1,y1 = y1,x2 = cell_coordinates_by_phenotype[[phenotype]][row_vals!=row,"Cell_X_Position"],y2 = cell_coordinates_by_phenotype[[phenotype]][row_vals!=row,"Cell_Y_Position"]) )
# ### ENd GetNearestNeighborByPhenotype check

CalculateCellNearestNeighborbyImageinPatientFilePathTree = function( file_path_tree ) {
  for ( patient_id in names(file_path_tree) ) {
    cat("\tpatient_id: ",patient_id,"\n")
    patient_df_list = list()
    for ( patient_image_id in names(file_path_tree[[patient_id]]) ) {
      cat("\t\timage_id: ",patient_image_id,"\n")
      cat("\t\t\tcomputing nearest neighbors ... \n")
      
      # read cell seg data
      cell_seg_data_file_path = GetCellSegDataFilePath(patient_id = patient_id,patient_image_id = patient_image_id,patient_file_path_tree = file_path_tree)
      cell_seg_data = read.csv(file = cell_seg_data_file_path,sep = "\t")
      
      # compute
      cell_coordinates_by_phenotype = GetCellCoordinatesByPhenotype(cell_seg_data = cell_seg_data,phenotypes_to_count = phenotypes_to_count)
      cell_distance_matrix = data.frame(matrix(data = numeric(),nrow = nrow(cell_seg_data),ncol = 3*length(phenotypes_to_count)))
      for ( i in 1:nrow(cell_seg_data) ) {
        if ( i%%250==0 ) { cat(sprintf(fmt = "\t\t\t\t @row: %5i / %5i\n",i,nrow(cell_seg_data))) }
        nearest_neighbor_of_cell_i = GetNearestNeighborByPhenotype(x = cell_seg_data[i,"Cell_X_Position"],y = cell_seg_data[i,"Cell_Y_Position"],cell_coordinates_by_phenotype = cell_coordinates_by_phenotype,exclude_row = i)
        cell_distance_matrix[i,] = as.vector(unlist(nearest_neighbor_of_cell_i))
      }
      names(cell_distance_matrix) = paste("Dist_",names(nearest_neighbor_of_cell_i),sep = "")
      # write
      write.table(x = cbind(cell_seg_data,cell_distance_matrix),file = cell_seg_data_file_path,sep = "\t",row.names = FALSE)
      cat("\t\t\tdone\n")
    }
  }
}

# ### START CalculateCellNearestNeighborbyImageinPatientFilePathTree end
# # file path tree
# file_path_tree = ConstructPatientFilePathTree(dir_name = "G:\\Script Rebuild\\Training Datasets\\Training_Data_P")
# ### END CalculateCellNearestNeighborbyImageinPatientFilePathTree end

CombineNearestNeighborsbyPatient = function( file_path_tree , cell_distance_by_patient_data_file_path) {
  for ( patient_id in names(file_path_tree) ) {
    cat("\tpatient_id: ",patient_id,"\n")
    data_by_patient = list()
    # read in associated image files
    for ( patient_image_id in names(file_path_tree[[patient_id]]) ) {
      cat("\t\timage_id: ",patient_image_id,"\n")
      cat("\t\t\treading in data ... ")
      data_by_patient[[patient_image_id]] = read.csv(
        file = GetCellSegDataFilePath(patient_id = patient_id,
                                      patient_image_id = patient_image_id,
                                      patient_file_path_tree = file_path_tree),
                                      sep= "\t",check.names = F)
      cat("done\n")
    }
    # row bind 
    cat("\t\t\tcombining data ... \n")
    combined_data = data_by_patient[[1]]
    cat(sprintf(fmt = "\t\t\t\t%1i / %1i",1,length(data_by_patient)),"\n")
    if ( length(data_by_patient) > 1 ) {
      for ( index in 2:length(data_by_patient) ) {
        cat(sprintf(fmt = "\t\t\t\t%1i / %1i",index,length(data_by_patient)),"\n")
        combined_data = rbind(combined_data,data_by_patient[[index]])
      }
    }
    cat("\t\t\tdone\n")
    # write 
    cat("\t\t\twriting data to [ ",paste(cell_distance_by_patient_data_file_path,patient_id,"_Combined_Data.txt",sep = ""), "] ... ")
    write.table(
      x = combined_data,
      file = paste(cell_distance_by_patient_data_file_path,patient_id,"_Combined_Data.txt",sep = ""),
      sep = "\t",row.names = FALSE)
    cat("done\n")
  }
}

# ### START CombineNearestNeighborsbyPatient start
# # file path tree
# file_path_tree = ConstructPatientFilePathTree(dir_name = "G:\\Script Rebuild\\2016_09_01_12_12_Run\\Datasets_P")
# # dir path
# cell_distance_by_patient_data_file_path = "G:\\Script Rebuild\\2016_09_01_12_12_Run\\Datasets_P\\DistanceOutput\\"
# # combine
# x = CombineNearestNeighborsbyPatient(file_path_tree = file_path_tree,cell_distance_by_patient_data_file_path = cell_distance_by_patient_data_file_path)
# ### END CombineNearestNeighborsbyPatient end

GetNearestNeighborCoordinatesByPhenotype = function( x, y, cell_coordinates_by_phenotype, exclude_row = 0 ) {
  min_dist_by_phenotype = list()
  for ( name in names(cell_coordinates_by_phenotype) ) {
    if ( nrow(cell_coordinates_by_phenotype[[name]]) == 0 ) {
      min_dist_by_phenotype[[name]]$Distance = c(NA);
      M = matrix(data = c(NA,NA),nrow=1,ncol=2); colnames(M) = c("Cell_X_Position","Cell_Y_Position");
      min_dist_by_phenotype[[name]]$Coordinates = M
    } else if ( exclude_row != 0 ) {
      row_vals = as.integer(as.character(row.names(cell_coordinates_by_phenotype[[name]])))
      x2 = cell_coordinates_by_phenotype[[name]]$Cell_X_Position[row_vals!=exclude_row]
      y2 = cell_coordinates_by_phenotype[[name]]$Cell_Y_Position[row_vals!=exclude_row]
      if ( length(x2) == 0 || length(y2) == 0 ) {
        min_dist_by_phenotype[[name]]$Distance = c(NA);
        M = matrix(data = c(NA,NA),nrow=1,ncol=2); colnames(M) = c("Cell_X_Position","Cell_Y_Position");
        min_dist_by_phenotype[[name]]$Coordinates = M
      } else {
        distances = Distance( x1 = x,y1 = y, x2 = x2, y2 = y2 )
        min_distance = min(distances)
        min_dist_by_phenotype[[name]]$Distance = min(distances)
        coordinates = cell_coordinates_by_phenotype[[name]][which(distances==min_distance),]
        min_dist_by_phenotype[[name]]$Coordinates = coordinates 
      }
    } else {
      distances = Distance( x1 = x,y1 = y, x2 = cell_coordinates_by_phenotype[[name]]$Cell_X_Position, y2 = cell_coordinates_by_phenotype[[name]]$Cell_Y_Position) 
      min_distance = min(distances)
      min_dist_by_phenotype[[name]]$Distance = min(distances)
      coordinates = cell_coordinates_by_phenotype[[name]][which(distances==min_distance),]
      min_dist_by_phenotype[[name]]$Coordinates = coordinates 
    }
  }
  return( min_dist_by_phenotype )
}

# ### START GetNearestNeighborPositionByPhenotype check
# rm(list=ls())
# cell_seg_data = read.csv(file = "F:\\Script Rebuild\\2016_09_01_12_12_Run\\Datasets_P\\Data\\Batch 8-18\\cu-04-06 D_1_cell_seg_data.txt",sep = "\t")
# cell_coordinates_by_phenotype = GetCellCoordinatesByPhenotype(cell_seg_data = cell_seg_data,phenotypes_to_count = phenotypes_to_count)
# row = 6;
# nearest_neighbor_coordinate_by_phenotype = GetNearestNeighborCoordinatesByPhenotype(x = cell_seg_data[row,"Cell_X_Position"],y = cell_seg_data[row,"Cell_Y_Position"],cell_coordinates_by_phenotype = cell_coordinates_by_phenotype,exclude_row = row)
# # check
# phenotype = "CD3+HLA-DR-"
# row_vals = as.integer(as.character(row.names(cell_coordinates_by_phenotype[[phenotype]])))
# x1 = cell_seg_data[row,"Cell_X_Position"]; y1 = cell_seg_data[row,"Cell_Y_Position"];
# min( Distance(x1 = x1,y1 = y1,x2 = cell_coordinates_by_phenotype[[phenotype]][row_vals!=row,"Cell_X_Position"],y2 = cell_coordinates_by_phenotype[[phenotype]][row_vals!=row,"Cell_Y_Position"]) )
# ### ENd GetNearestNeighborPositionByPhenotype check

# FlourescentBaseMarkerColorCoding = list(
#   "CD3+" = '#E41A1C',
#   "CD3-CD8+" = '#377EB8',
#   "Sox10+" = '#4DAF4A',
#   "CD68+" = '#984EA3',
#   "Other" = '#FF7F00'
# )
# 
# library(tiff)
# readTIFF()



