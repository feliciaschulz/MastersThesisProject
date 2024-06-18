##############################################
# CELL TYPE INFERENCE - create_data.R
# 
# This script takes cell intensity data frames as input csv files and predicts
# cell types for each cell based on its intensity values for different markers.
# It uses an input csv file with thresholds for each marker that decide at which 
# intensity value the signal is seen as true and positive.
# Based on that, it creates a binary data frame with positive / negative expressions
# for each cell for each marker. 
# Then, a binary marker expression profiles file is used as input to determine which cell 
# type profile fits each cell.
# The final inferred cell type is then written to a data frame along with area, 
# centroids and intensity information. It is then saved.
#
# Usage: 
# Run the script either line-by-line in RStudio or source it (top right corner button 
# in RStudio). If you want to change the input image, change the string assigned to
# the variable "my_img" below.
# 
# Input: 
#   - Data frames with cells in rows, biomarker intensity data in columns,
#     two columns for x and y coordinates (X_coord, Y_coord), area of cell (area)
#     from folder "../00_Data/IntensityDataFrames/"
#   - Thresholds csv file with biomarkers as rows and threshold minimum values as 
#     columns from "../00_Data/Thresholds/"
#   - Marker Expression Profiles csv file with cell types as rows, biomarkers as columns
#     and binary values indicating positive or negative expression from "../00_Data/"
#
# Output: Two csv files in the folder "../00_Data/IntensitiesCelltypes/"
#   - ..._intensitiesdf.csv: intensity data as well as centroids, area, inferred cell type
#   - ..._celltypesdf.csv: binary marker expression data as well as inferred cell type
#
# Functions: Function create_data() that returns the two data frames as a list
###############################################


create_data <- function(markers = c("actin", "betatubulin3", "cd138", "cd248", "cd279", "cd3", "cd34", "cd38", "cd4",
                                "cd44", "cd45", "cd45ro", "cd79a", "collageni", "cytokeratin",
                                "fibulin2", "hladrdpdq", "ki67", "lyve1", "pdgfra", "podoplanin"),
                    thresholds_data = " ",
                    marker_celltypes_data = "",
                    rel_path) {
  
  
  print(paste("You are computing dimensionality reduction methods for the markers ", list(markers)))
  
  #####################
  #### Intensities ####
  #####################
  
  
  # Create df "intensities": all intensities for all segmented cells
  # length(intensities)/num of cols = amount of biomarkers
  # num of rows = amount of segmented cells for this image
  for (marker in markers) {   # Loop through selected markers
    
    file_name <- paste0("../00_Data/IntensityDataFrames/", rel_path, "/", marker, "df.csv")  # Construct the file name
    df <- read.csv(file_name) # Read the CSV file and save to df
    if (marker == markers[1]) {
      intensities = data.frame(Cells = df$X) # Initialise intensities df with the first marker
    }
    intensities[[marker]] <- df$intensity_max # Save to data frame intensitites
    
  }
  intensities$Cells =  NULL # Remove cells column
  

  ####################
  #### Cell types ####
  ####################
  
  ## Make cell types data frame
  
  # Create new df celltypes with each cell as row and column names as "intensities"
  celltypes <- data.frame(matrix(nrow = length(rownames(intensities)), ncol = length(colnames(intensities))))
  rownames(celltypes) <- rownames(intensities)
  colnames(celltypes) <- colnames(intensities)
  celltypes["InferredCellType"] = NA
  
  # Get thresholds as created in MACSiQ View with intensity histogram settings
  thresholds = read.csv(paste0("../00_Data/Thresholds/", rel_path, "_Thresholds.csv"))
  thresholds_clean <- subset(thresholds, !grepl("DAPI", Name) & !grepl("background", Name) & !grepl("AF ", Name)) # remove DAPIs, bleached images and autofluorescence
  
  # Change threshold marker names so that I can loop through the thresholds easily
  for (i in seq_len(nrow(thresholds_clean))) {
    words <- strsplit(thresholds_clean$Name[i], " ")[[1]]   # Split the name into words as a list
    if (length(words)>1) {
      words <- words[-length(words)]   # Remove the last word
    }
    words <- tolower(words)   # Convert all words to lowercase
    new_name <- paste(words, collapse = "")   # Concatenate remaining words
    
    # Assign the new name to the corresponding row in the dataframe
    thresholds_clean$Name[i] <- new_name
  }
  rownames(thresholds_clean) <- thresholds_clean[,1]
  thresholds_clean[,1] <- NULL
  
  # Import ground truth phenotypes depending on marker combinations, binary df
  # OBS: this df has the celltypes in [,1] because there are duplicate names so they cannot be made row names
  edited_marker_celltypes <- read.csv(paste0("../00_Data/", marker_celltypes_data), sep=";", header=TRUE)  
  row_names <- edited_marker_celltypes$X
  edited_marker_celltypes$X <- NULL
  missing_columns <- setdiff(names(intensities), names(edited_marker_celltypes)) # In case any are missing
  edited_marker_celltypes[missing_columns] <- 0 # Add missing columns to edited_marker_celltypes with all values set to zero
  edited_marker_celltypes <- edited_marker_celltypes[, names(intensities)] # Reorder columns to match intensities dataframe
  edited_marker_celltypes <- cbind(X = row_names, edited_marker_celltypes) # Make into one data frame

  

  # Assign binary marker expression values to celltypes data
  # OBS: intensities df marker names must match thresholds marker names, otherwise error in TRUE/FALSE comparison
  assign_binary_marker_expressions <- function(celltypes, intensities, thresholds, edited_marker_celltypes) {
    for (i in 1:nrow(intensities)) {
      
      # Loop through each intensities matrix cell to assign binary marker expressions (1 if above threshold, 0 if below threshold)
      for (j in 1:ncol(intensities)) {
  
        marker = colnames(intensities[j])
        
        if (intensities[i,j] >= thresholds_clean[marker,"RangeMin"] & intensities[i,j] <= thresholds_clean[marker, "RangeMax"]) {
          celltypes[i,j] <- 1
        } else {
          celltypes[i,j] <- 0
        }
      }
      
      cellrow = celltypes[i,][1:length(celltypes[i,])-1] # All columns except Inf.CellType
      
      # Loop through celltypes from binary edited_marker_celltypes decision data frames
      for (p in 1:nrow(edited_marker_celltypes)) {
        phenotyperow = edited_marker_celltypes[p,][2:length(edited_marker_celltypes)] # Row from edited marker celltypes without column X (=celltype name)
        
        if (all(cellrow == phenotyperow)) { # Check if celltype row is the same as one of these
          celltypes[i,"InferredCellType"] = edited_marker_celltypes[p,][1] # add cell type to appropriate column in celltypes df
        }
        
      }
      
      if (is.na(celltypes[i,"InferredCellType"])) { # This means that nothing has changed and no column has been added 
        celltypes[i,"InferredCellType"] = "Mixed" # Means that there were conflicting cell types (no matching row)
      }

    }
    return(celltypes)
  }
  celltypes <- assign_binary_marker_expressions(celltypes, intensities, thresholds_clean, edited_marker_celltypes)
  
  # Add in other data
  intensities$area = df$area
  intensities$X_coord = df$centroid.1
  intensities$Y_coord = df$centroid.0
  intensities$InferredCellType = celltypes$InferredCellType

  print("Done.")
  return(list(intensities, celltypes))
}


###############################
######## RUN CODE HERE ########
###############################

# Define image name
my_img = "R1B1ROI1"

# Run create_data() function
datalist <- create_data(markers = c("actin", "cd3", "cd4",
                                "cd45", "cd45ro", "collageni", "cytokeratin",
                                "fibulin2", "lyve1", "podoplanin", "cd38", "cd138"),
                    marker_celltypes_data = "MarkerExpressionProfiles.csv",
                    rel_path = my_img)

# Save output to variables (function returns a list of two objects)
intensities <- datalist[[1]]
celltypes <- datalist[[2]]

r1d1roi2_ct = celltypes


# Save to folder
write.csv(intensities, paste0("../00_Data/IntensitiesCelltypes/", my_img, "_intensitiesdf.csv"))
write.csv(celltypes, paste0("../00_Data/IntensitiesCelltypes/", my_img, "_celltypesdf.csv"))











