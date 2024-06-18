##############################################
# COUNT CELL TYPE ABUNDANCES - count_instances.R
# 
# This script takes celltypes data frames with binary expression for each marker
# and creates one data frame with total cell type counts.
# The output data frame (combined_df) also shows the biomarker counts that contribute
# to each cell type.
#
# Usage: 
# Run the script either line-by-line in RStudio or source it (top right corner button 
# in RStudio). If you want to change the input image, change the string assigned to
# the variable "my_img" below.
# 
# Input: 
#   - Data frame with cells in rows, biomarker binary expression data in columns,
#     inferred cell type in last column from "../00_Data/IntensitiesCelltypes/" 
#     called "*_celltypesdf.csv"
#
# Output: One csv file in the folder "00_Data/" called "*celltypeCounts.csv"
###############################################

onlyOneDF = TRUE

if (onlyOneDF == FALSE) {
  #Load in data and combine to one data frame
  r1b1roi1_ct = read.csv("../00_Data/IntensitiesCelltypes/R1B1ROI1_celltypesdf.csv")
  r1b1roi2_ct = read.csv("../00_Data/IntensitiesCelltypes/R1B1ROI2_celltypesdf.csv")
  r1c1roi1_ct = read.csv("../00_Data/IntensitiesCelltypes/R1C1ROI1_celltypesdf.csv")
  r1c1roi2_ct = read.csv("../00_Data/IntensitiesCelltypes/R1C1ROI2_celltypesdf.csv")
  r1d1roi1_ct = read.csv("../00_Data/IntensitiesCelltypes/R1D1ROI1_celltypesdf.csv")
  r1d1roi2_ct = read.csv("../00_Data/IntensitiesCelltypes/R1D1ROI2_celltypesdf.csv")
  celltypes = rbind(r1b1roi1_ct, r1b1roi2_ct, r1c1roi1_ct, r1c1roi2_ct, r1d1roi1_ct, r1d1roi2_ct)
} else {
  img_name = "R1B1ROI1"
  r1b1roi1_ct = read.csv(paste0("../00_Data/IntensitiesCelltypes/", img_name, "_celltypesdf.csv"))
  celltypes = r1b1roi1_ct
}


# Make table for counts of celltypes
celltype_counts <- table(celltypes$InferredCellType)

# Make table for biomarker expressions
markers <- names(celltypes[-length(celltypes)]) # all except inf.celltype

# Function to count positive occurrences of markers for each inferred cell type
count_positives <- function(celltypes, marker) {
  aggregate(celltypes[[marker]], by = list(celltypes$InferredCellType), FUN = sum)
}

# Create an empty list to store tables for each marker
marker_table <- list()

# Iterate over each marker and count positive occurrences for each inferred cell type
for (marker in markers) {
  marker_table[[marker]] <- count_positives(celltypes, marker)
}

# Print the tables
for (marker in markers) {
  cat("Table for marker:", marker, "\n")
  marker_table[[marker]]$InferredCellType = marker_table[[marker]]$Group.1 
  marker_table[[marker]]$Group.1 = NULL
  print(marker_table[[marker]])
  cat("\n")
}

# Initialise data frame
combined_df <- NULL
for (marker in markers) { # Loop through markers and add to data frame
  if (is.null(combined_df)) {
    combined_df <- marker_table[[marker]]
  } else {
    combined_df <- merge(combined_df, marker_table[[marker]], by = "InferredCellType", all = TRUE)
    colnames(combined_df)[which(colnames(combined_df) == "x")] <- marker
  }
}

# Correct columns (first two have the wrong names)
combined_df$actin = combined_df$x.y
combined_df$x.x = NULL
combined_df$x.y = NULL
combined_df$total = celltype_counts # Add column counts

# Make cell types into rownames 
rownames(combined_df) = combined_df$InferredCellType
combined_df$InferredCellType = NULL

# Add total counts row that sums up each marker column
totals <- colSums(combined_df) #Calculate the sum for each column
total_row <- as.data.frame(t(totals)) #Convert the sums to a data frame with one row
row.names(total_row) <- "Total" # This will be the row name 
combined_df <- rbind(combined_df, total_row) # Add row to total data frame

# Open df in new tab
View(combined_df)

# Save to folder
write.csv(combined_df, file = paste0("../00_Data/",img_name, "celltypeCounts.csv"), row.names = TRUE)
