
create_data <- function(markers = c("actin", "betatubulin3", "cd138", "cd248", "cd279", "cd3", "cd34", "cd38", "cd4",
                                "cd44", "cd45", "cd45ro", "cd79a", "collageni", "cytokeratin",
                                "fibulin2", "hladrdpdq", "ki67", "lyve1", "pdgfra", "podoplanin"),
                    thresholds_data = "export-image-thresholds.csv",
                    marker_celltypes_data = "Marker_Celltypes_ToEdit.csv",
                    rel_path,
                    cyto) {
  
  
  all_markers = c("actin", "betatubulin3", "cd138", "cd248", "cd279", "cd3", "cd34", "cd38", "cd4",
                  "cd44", "cd45", "cd45ro", "cd79a", "collageni", "cytokeratin", "dapi",
                  "fibulin2", "hladrdpdq", "ki67", "lyve1", "pdgfra", "podoplanin")
  
  print(paste("You are computing dimensionality reduction methods for the markers",
              list(markers), "and NOT for the markers", setdiff(markers, all_markers)))
  
  ############################
  #### Intensity averages ####
  ############################
  
  
  for (marker in markers) { 
    
    # if (cyto==TRUE) {
    #   extension="_cytodf.csv"
    # } else (
    #   extension="df.csv"
    # )
    
    extension="df.csv"
    
    file_name <- paste0(rel_path, "/", marker, extension)  # Construct the file name
    #assign(marker, read.csv(file_name))    # Read the CSV file
    # Read the CSV file
    df <- read.csv(file_name)
    
    if (marker == markers[1]) {
      averages = data.frame(Cells = df$X)
      area_data <- df$Area
    }
    
    # Assign the mean to the averages data frame
    averages[[marker]] <- df$Mean
  }
  averages$Cells =  NULL # remove cells column
  # Create df "averages": all intensity averages for all segmented cells
  # length(averages)/num of cols = amount of biomarkers
  # num of rows = amount of segmented cells for this image
  
  
  
  
  
  ####################
  #### Cell types ####
  ####################
  
  ## Make cell types data
  
  # Create new df celltypes with each cell as row and column names as "averages"
  celltypes <- data.frame(matrix(nrow = length(rownames(averages)), ncol = length(colnames(averages))))
  rownames(celltypes) <- rownames(averages)
  colnames(celltypes) <- colnames(averages)
  
  
  # Get thresholds as created in MACSiQ View with intensity histogram settings
  thresholds = read.csv(paste0(rel_path, "/", rel_path, "_Thresholds.csv"))
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
  
  
  
  # Assign binary marker expression values to celltypes data
  # OBS: averages df marker names must match thresholds marker names, otherwise error in TRUE/FALSE comparison
  assign_binary_marker_expressions <- function(celltypes, averages, thresholds) {
    for (i in 1:nrow(averages)) {
      
      for (j in 1:ncol(averages)) {
        
        marker = colnames(averages[j])
        
        if (averages[i,j] >= thresholds_clean[marker,"DispRangeMin"]) {
          celltypes[i,j] <- 1
        } else {
          celltypes[i,j] <- 0
        }
        
      }
    }
    return(celltypes)
  }
  celltypes <- assign_binary_marker_expressions(celltypes, averages, thresholds_clean)
  
  # Import ground truth phenotypes depending on marker combinations, binary df
  # Should have been created in create_celltype_GT.R, and manually edited in the Mac app Numbers
  # then exported as csv and imported here
  edited_marker_celltypes <- read.csv(marker_celltypes_data, sep=";", header=TRUE)
  
  # OBS: this df has the celltypes in [,1] because there are duplicate names so they cannot be made row names
  
  
  #### OBS: Doing this code because the edited_marker_celltypes have not been properly updated yet. TEMPORARY CODE
  row_names <- edited_marker_celltypes$X
  edited_marker_celltypes$X <- NULL
  missing_columns <- setdiff(names(averages), names(edited_marker_celltypes))
  edited_marker_celltypes[missing_columns] <- 0 # Add missing columns to edited_marker_celltypes with all values set to zero
  edited_marker_celltypes <- edited_marker_celltypes[, names(averages)] # Reorder columns to match averages dataframe
  edited_marker_celltypes <- cbind(X = row_names, edited_marker_celltypes)
  # Print the updated dataframe
  head(edited_marker_celltypes)
  
  
  
  # Loop through celltypes df to add inferred cell type from marker_celltypes df
  assign_celltypes <- function(celltypes, edited_marker_celltypes) {
    celltypes["InferredCellType"] = NA
    for (i in 1:nrow(celltypes)) {
      cellrow = celltypes[i,][1:length(celltypes[i,])-1]
      
      for (p in 1:nrow(edited_marker_celltypes)) {
        phenotyperow = edited_marker_celltypes[p,][2:length(edited_marker_celltypes)]
        
        if (all(cellrow == phenotyperow)) {
          celltypes[i,"InferredCellType"] = edited_marker_celltypes[p,][1]
        }
        
      }
      if (is.na(celltypes[i,"InferredCellType"])) { # this means that nothing has changed and no column has been added
        celltypes[i,"InferredCellType"] = "Mixed"
      }
    }
    return(celltypes)
  }
  # DONT run again, takes too long
  # test_celltypes = head(celltypes, 2)
  # test_celltypes <- assign_celltypes(test_celltypes, edited_marker_celltypes)
  celltypes <- assign_celltypes(celltypes, edited_marker_celltypes)
  
  
  
  return(list(averages, celltypes, area_data, cyto, thresholds_clean))
  
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
                                "fibulin2", "lyve1", "podoplanin"),
                    thresholds_data = "export-image-thresholds.csv",
                    marker_celltypes_data = "New_Marker_Celltypes_ToEdit.csv",
                    rel_path = "R1B1ROI1",
                    cyto = TRUE)
averages <- datalist[[1]]
celltypes <- datalist[[2]]
area_data <- datalist[[3]]
averages_normalised <- datalist[[4]]
cyto <- datalist[[5]]
thresholds_clean <- datalist[[6]]





# "Madelenes selection": "actin", "cd3", "cd4",
#   "cd45", "cd45ro", "collageni", "cytokeratin",
#   "fibulin2", "lyve1", "podoplanin"














  