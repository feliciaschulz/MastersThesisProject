
create_data <- function(markers = c("actin", "betatubulin3", "cd138", "cd248", "cd279", "cd3", "cd34", "cd38", "cd4",
                                "cd44", "cd45", "cd45ro", "cd79a", "collageni", "cytokeratin",
                                "fibulin2", "hladrdpdq", "ki67", "lyve1", "pdgfra", "podoplanin"),
                    thresholds_data = " ",
                    marker_celltypes_data = "",
                    rel_path) {
  
  
  all_markers = c("actin", "betatubulin3", "cd138", "cd248", "cd279", "cd3", "cd34", "cd38", "cd4",
                  "cd44", "cd45", "cd45ro", "cd79a", "collageni", "cytokeratin", "dapi",
                  "fibulin2", "hladrdpdq", "ki67", "lyve1", "pdgfra", "podoplanin")
  
  print(paste("You are computing dimensionality reduction methods for the markers",
              list(markers), "and NOT for the markers", setdiff(markers, all_markers)))
  
  ############################
  #### Intensity averages ####
  ############################
  
  
  # Create df "averages": all intensity averages for all segmented cells
  # length(averages)/num of cols = amount of biomarkers
  # num of rows = amount of segmented cells for this image
  for (marker in markers) {   # Loop through selected markers
    
    file_name <- paste0("../00_Data/IntensityDataFrames/", rel_path, "/", marker, "df.csv")  # Construct the file name
    df <- read.csv(file_name) # Read the CSV file and save to df
    print("df for marker:")
    print(marker)
    print(head(df))
    if (marker == markers[1]) {
      averages = data.frame(Cells = df$X) # Initialise averages df with the first marker
    }
    averages[[marker]] <- df$intensity_max # Save to data frame averages
    
  }
  averages$Cells =  NULL # Remove cells column

  print("First step done")
  ####################
  #### Cell types ####
  ####################
  
  ## Make cell types data frame
  
  # Create new df celltypes with each cell as row and column names as "averages"
  celltypes <- data.frame(matrix(nrow = length(rownames(averages)), ncol = length(colnames(averages))))
  rownames(celltypes) <- rownames(averages)
  colnames(celltypes) <- colnames(averages)
  celltypes["InferredCellType"] = NA
  
  print("trying to find thresholds data")
  # Get thresholds as created in MACSiQ View with intensity histogram settings
  thresholds = read.csv(paste0("../00_Data/Thresholds/", rel_path, "_Thresholds.csv"))
  thresholds_clean <- subset(thresholds, !grepl("DAPI", Name) & !grepl("background", Name) & !grepl("AF ", Name)) # remove DAPIs, bleached images and autofluorescence
  print("found thresholds data")
  
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
  
  print("thresholds have been fixed")
  print("triyng to find marker celltypes")
  
  # Import ground truth phenotypes depending on marker combinations, binary df
  # OBS: this df has the celltypes in [,1] because there are duplicate names so they cannot be made row names
  edited_marker_celltypes <- read.csv(marker_celltypes_data, sep=";", header=TRUE)  
  #### OBS: Doing this code because the edited_marker_celltypes have not been properly updated yet. TEMPORARY CODE
  row_names <- edited_marker_celltypes$X
  edited_marker_celltypes$X <- NULL
  missing_columns <- setdiff(names(averages), names(edited_marker_celltypes))
  edited_marker_celltypes[missing_columns] <- 0 # Add missing columns to edited_marker_celltypes with all values set to zero
  edited_marker_celltypes <- edited_marker_celltypes[, names(averages)] # Reorder columns to match averages dataframe
  edited_marker_celltypes <- cbind(X = row_names, edited_marker_celltypes)
  
  print("marker celtypes found")
  print("assigning celltypes")
  ### FIRST LOOP THAT TAKES LONG
  # Assign binary marker expression values to celltypes data
  # OBS: averages df marker names must match thresholds marker names, otherwise error in TRUE/FALSE comparison
  assign_binary_marker_expressions <- function(celltypes, averages, thresholds, edited_marker_celltypes) {
    for (i in 1:nrow(averages)) {
      
      # Loop through each averages matrix cell to assign binary marker expressions (1 if above threshold, 0 if below threshold)
      for (j in 1:ncol(averages)) {
  
        marker = colnames(averages[j])
        
        if (averages[i,j] >= thresholds_clean[marker,"RangeMin"] & averages[i,j] <= thresholds_clean[marker, "RangeMax"]) {
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
  celltypes <- assign_binary_marker_expressions(celltypes, averages, thresholds_clean, edited_marker_celltypes)
  
  return(list(averages, celltypes, thresholds_clean))
  
}


###############################
######## RUN CODE HERE ########
###############################
# rm(list=ls())

# You can choose from these markers: (default = all except dapi)
# "actin", "betatubulin3", "cd138", "cd248", "cd279", "cd3", "cd34", "cd38", "cd4",
# "cd44", "cd45", "cd45ro", "cd79a", "collageni", "cytokeratin", "dapi",
# "fibulin2", "hladrdpdq", "ki67", "lyve1", "pdgfra", "podoplanin"

datalist <- create_data(markers = c("actin", "cd3", "cd4",
                                "cd45", "cd45ro", "collageni", "cytokeratin",
                                "fibulin2", "lyve1", "podoplanin", "cd38", "cd138"),
                    thresholds_data = "export-image-thresholds.csv",
                    marker_celltypes_data = "New_Marker_Celltypes.csv",
                    rel_path = "R1D1ROI1")
averages <- datalist[[1]]
celltypes <- datalist[[2]]
thresholds_clean <- datalist[[3]]





# "Madelenes selection": "actin", "cd3", "cd4",
#   "cd45", "cd45ro", "collageni", "cytokeratin",
#   "fibulin2", "lyve1", "podoplanin"














  