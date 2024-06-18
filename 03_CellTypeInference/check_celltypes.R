##############################################
# CREATE CELLTYPES SCATTERPLOTS - check_celltypes.R
# 
# This script takes files with the inferred cell types of cells as well as their
# coordinates and creates scatterplots coloured by cell types.
# The user can either choose to only create a plot for one image by having the 
# variable "onlyOneDF" == TRUE, or to create a grid of multiple plots for
# multiple images. For just one image, the image name must be changed in "img_name",
# and for multiple images, all paths and filenames need to be checked.
# The plot axes are in micrometres.
#
# Usage: 
# Run the script line-by-line in RStudio or click "Source" in the top right corner. 
# If you want to change the input image, change the string assigned to
# the variable "my_img" below. 
#
# Required packages & versions: 
#   - ggplot2: 3.5.0
#   - dplyr: 1.1.4
#   - patchwork: 1.2.0
# 
# Input: 
#   - Raw data frame with cells in rows, biomarker intensity expression data in columns, area,
#     x and y coordinates, and inferred cell type in last column from "../00_Data/IntensitiesCelltypes/" 
#     called "*_intensitiesdf.csv"
#
#
# Output: Nothing, depending on settings chosen by user, scatterplots can be inspected. Can be saved by
# clicking "Export" in RStudio above the plot.
###############################################

# Load the ggplot2 package
library(ggplot2)
library(dplyr)
library(patchwork)

onlyOneDF = TRUE
img_name = "R1B1ROI1"

# Define the color palette for the cell types
color_palette <- c("Epithelial cells"= "orange",     
                   "Extracellular matrix"= "cyan",
                   "Lymphatic vessel"= "yellow",
                   "Mixed"= "darkgrey",
                   "Myocytes"= "red",
                   "Other"= "lightgrey",
                   "T-cells Memory"= "magenta",
                   "T-cells CD4+ (helper)"= "blue",
                   "T-cells Naive"= "navy",
                   "T-cells Mixed"= "purple",
                   "B-cells Plasma"= "lightgreen")


if (onlyOneDF == FALSE) {
  r1b1roi1 = read.csv("../Old00_Data/IntensitiesCelltypes/R1B1ROI1_intensitiesdf.csv")
  r1b1roi2 = read.csv("../Old00_Data/IntensitiesCelltypes/R1B1ROI2_intensitiesdf.csv")
  r1c1roi1 = read.csv("../Old00_Data/IntensitiesCelltypes/R1C1ROI1_intensitiesdf.csv")
  r1c1roi2 = read.csv("../Old00_Data/IntensitiesCelltypes/R1C1ROI2_intensitiesdf.csv")
  r1d1roi1 = read.csv("../Old00_Data/IntensitiesCelltypes/R1D1ROI1_intensitiesdf.csv")
  r1d1roi2 = read.csv("../Old00_Data/IntensitiesCelltypes/R1D1ROI2_intensitiesdf.csv")
  
  # List of data frames with their names
  data_frames <- list(R1B1ROI1 = r1b1roi1, 
                      R1B1ROI2 = r1b1roi2, 
                      R1C1ROI1 = r1c1roi1, 
                      R1C1ROI2 = r1c1roi2, 
                      R1D1ROI1 = r1d1roi1, 
                      R1D1ROI2 = r1d1roi2)
  
  
  # Create a list to store plots without legends
  plots <- lapply(names(data_frames), function(name) {
    df <- data_frames[[name]]
    ggplot(df, aes(x = X_coord, y = Y_coord, color = InferredCellType)) +
      geom_point(size = 0.8) +  # Adjust point size as needed
      scale_color_manual(values = color_palette) +  # Assign colors based on the palette
      scale_x_continuous(labels = function(x) x * 0.17) + # Scale to micrometres instead of pixels
      scale_y_reverse(labels = function(y) y * 0.17) +  # Reverse and change y-axis labels
      labs(x = NULL, y = NULL, title = name) +  # Set title to the name of the data frame
      theme_minimal() +
      theme(legend.position = "none") 
  })
  
  
  # Create a plot to extract the legend
  legend_plot <- ggplot(r1b1roi1, aes(x = X_coord, y = Y_coord, color = InferredCellType)) +
    geom_point(size = 2.5) +  # Adjust point size as needed
    scale_color_manual(values = color_palette) +  # Assign colors based on the palette
    theme(legend.position = "right", legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 8)))  # Adjust legend point size
  
  # Extract the legend
  legend <- cowplot::get_legend(legend_plot)
  
  # Combine plots and add the legend at the bottom
  final_plot <- (plots[[1]] | plots[[2]]) /
    (plots[[3]] | plots[[4]]) /
    (plots[[5]] | plots[[6]]) +
    plot_layout(guides="collect") &
    theme(legend.position = "right")

} else if (onlyOneDF == TRUE) { # If just one input data frame
  data_frame <- read.csv(paste0("../00_Data/IntensitiesCelltypes/", img_name, "_intensitiesdf.csv"))
  
  # Create the scatter plot - JUST ONE
  final_plot <- ggplot(data_frame, aes(x = X_coord, y = Y_coord, color = InferredCellType)) +
    geom_point(size = 0.8) +  # Adjust point size as needed
    scale_color_manual(values = color_palette) +  # Assign colors based on the palette
    labs(x = NULL, y = NULL, title = img_name) +
    theme_minimal() +
    scale_x_continuous(labels = function(x) x * 0.17) + # Scale to micrometres instead of pixels
    scale_y_reverse(labels = function(y) y * 0.17) +  # Reverse and change y-axis labels
    theme(legend.position = "right")  # Adjust legend position as needed
  
}


# print plot
print(final_plot)



