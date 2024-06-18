##############################################
# CREATE HEATMAPS - heatmaps_celltypeintensities.R
# 
# This script takes intensities data frames and creates heatmaps showing 
# intensities, organised by markers horizontally and cell types vertically.
# Given the user has already created normalised intensity data frames with 
# data_exploration.R, they can also change the variable "norm" to TRUE and thereby
# use the normalised data as input instead.
# Additionally, the user can also choose to only include lymphocytes in the heatmap.
# This will make it easier to see the individual cells of the lymphocytes.
#
# Usage: 
# Run the script line-by-line in RStudio. If you want to change the input image, change the string assigned to
# the variable "my_img" below. If you want to use the normalised data as input, change the variable
# "norm" to TRUE. If you want to only plot lymphocytes in the heatmap, change plot_lymph_separate to TRUE.
#
# Required packages: ComplexHeatmap (v. 2.18.0)
# 
# Input: 
#   - Raw data frame with cells in rows, biomarker intensity expression data in columns, area,
#     x and y coordinates, and inferred cell type in last column from "../00_Data/IntensitiesCelltypes/" 
#     called "*_intensitiesdf.csv"
#   - Normalised data frame with the cells in rows, biomarker intensity expression data in columns, area,
#     x and y coordinates, and inferred cell type in last column from "../00_Data/Normalised/" 
#     called "*_norm.csv"
#
# Output: Nothing, depending on settings chosen by user, heatmaps can be inspected. Can be saved by
# clicking "Export" in RStudio above the plot.
###############################################
library(ComplexHeatmap)
set.seed(20)


onlyOneDF = TRUE # Leave this option TRUE if you only have one input data frame available
norm = TRUE # Change this option if you want to see the normalised heatmap instead
plot_lymph_separate = FALSE # Change this option if you want to plot only the lymphocytes
                            # separately in order to see them better (they become bigger)


if (norm == FALSE) {
  title_norm = "raw" # This is for the heatmap title later
  if (onlyOneDF == FALSE) {
    #Load in data and combine to one data frame
    r1b1roi1 = read.csv("../Old00_Data/IntensitiesCelltypes/R1B1ROI1_intensitiesdf.csv")
    r1b1roi2 = read.csv("../Old00_Data/IntensitiesCelltypes/R1B1ROI2_intensitiesdf.csv")
    r1c1roi1 = read.csv("../Old00_Data/IntensitiesCelltypes/R1C1ROI1_intensitiesdf.csv")
    r1c1roi2 = read.csv("../Old00_Data/IntensitiesCelltypes/R1C1ROI2_intensitiesdf.csv")
    r1d1roi1 = read.csv("../Old00_Data/IntensitiesCelltypes/R1D1ROI1_intensitiesdf.csv")
    r1d1roi2 = read.csv("../Old00_Data/IntensitiesCelltypes/R1D1ROI2_intensitiesdf.csv")
    
    df = rbind(r1b1roi1, r1b1roi2)
    
  } else {
    img_name = "R1B1ROI1"
    df = read.csv(paste0("../Old00_Data/IntensitiesCelltypes/", img_name, "_intensitiesdf.csv"))
    df$X = NULL
  }
} else { # if norm == TRUE, so it will use normalised data frames
    title_norm = "normalised" # This is for the heatmap title later
    if (onlyOneDF == FALSE) {
      r1b1roi1_norm = read.csv("../Old00_Data/IntensitiesCelltypes/R1B1ROI1_intensitiesdf.csv")
      r1b1roi2_norm = read.csv("../Old00_Data/IntensitiesCelltypes/R1B1ROI2_intensitiesdf.csv")
      r1c1roi1_norm = read.csv("../Old00_Data/IntensitiesCelltypes/R1C1ROI1_intensitiesdf.csv")
      r1c1roi2_norm = read.csv("../Old00_Data/IntensitiesCelltypes/R1C1ROI2_intensitiesdf.csv")
      r1d1roi1_norm = read.csv("../Old00_Data/IntensitiesCelltypes/R1D1ROI1_intensitiesdf.csv")
      r1d1roi2_norm = read.csv("../Old00_Data/IntensitiesCelltypes/R1D1ROI2_intensitiesdf.csv")
      
    df = rbind(r1b1roi1_norm, r1b1roi2_norm )
  } else {
    img_name = "R1B1ROI1"
    df = read.csv(paste0("../00_Data/Normalised/", img_name, "_norm.csv"))
    df$X = NULL
  }
}



# Exclude cell type "Other" because it's not very interesting
df = df[df$InferredCellType != "Other",]

# Data needs to be a matrix for the heatmap package
data_matrix <- as.matrix(df[1:12]) # Make data frame into matrix, only including the 12 markers

# Order markers
ordered_columns <- sort(colnames(data_matrix)) # Get the column names and order them alphabetically

# Order cell types
celltype_order <- order(df$InferredCellType)
cellType = df$InferredCellType

# Get the corresponding cell type labels in the same order
celltype_labels <- df$InferredCellType[celltype_order]


# Plot the COMPLETE reordered heatmap with cell type labels
Heatmap(data_matrix, name ="intensity", row_order = celltype_order, 
        row_split = cellType, row_title_rot = 0, column_title=paste0("Heatmap for ", img_name, "(", title_norm, " intensities)"),  
        show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE,
        row_gap = unit(3, "mm"), column_order = ordered_columns)  # Use cell type labels as row side color




###### Plot lymphocytes separately
if (plot_lymph_separate == TRUE) {
  
  # Make lymphocytes vector
  lymphocytes = c("T-cells Memory", "T-cells CD4+ (helper)", "T-cells Mixed", "T-cells Naive", "B-cells Plasma")
  
  # Separate lymphocytes
  # First, take all that are not lymphocytes
  df_notl = df[!(df$InferredCellType %in% lymphocytes), ]
  celltype_order_notl <- order(df_notl$InferredCellType) # order cell types
  cellType_notl = df_notl$InferredCellType
  
  title_subset = "(not lymphocytes)"
  # Create heatmap - Run this command (after all the others above) to see the heatmap for all non-lymphocytes
  Heatmap(as.matrix(df_notl[1:12]), name ="intensity", row_order = celltype_order_notl, 
          row_split = cellType_notl, row_title_rot = 0, column_title=paste0("Heatmap for ", img_name, "(", title_norm, " intensities)"),  
          show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE)  # Use cell type labels as row side order
  
  
  # Now, same for only lymphocytes
  df_l = df[df$InferredCellType %in% lymphocytes, ]
  celltype_order_l <- order(df_l$InferredCellType)
  cellType_l = df_l$InferredCellType
  
  title_subset = "(only lymphocytes)"
  # Create heatmap - Run this command (after all the others above) to see the heatmap for only lymphocytes
  Heatmap(as.matrix(df_l[1:12]), name ="intensity", row_order = celltype_order_l, 
          row_split = cellType_l, row_title_rot = 0, column_title=paste0("Heatmap for ", img_name, " ", title_norm, " ", title_subset),  
          show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE)  # Use cell type labels as row side order
  
}




