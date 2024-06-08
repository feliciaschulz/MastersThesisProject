# Assuming you have your data frame "geo" with coordinates and a column "InferredCellType"
# Replace these with your actual data
geo <- read.csv("../00_Data/CellShapeData/R1B1ROI1_Cellshape.csv")
# Load the ggplot2 package
library(ggplot2)
library(dplyr)
library(patchwork)

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

grey_palette <- c("Epithelial cells"= "grey",     
                  "Extracellular matrix"= "grey",
                  "Lymphatic vessel"= "grey",
                  "Mixed"= "grey",
                  "Myocytes"= "green",
                  "Other"= "grey",
                  "T-cells Memory"= "grey",
                  "T-cells CD4+ (helper)"= "magenta",
                  "T-cells Naive"= "purple",
                  "T-cells Mixed"= "blue" )

#geo <- geo1
#geo$InferredCellType = celltypes$InferredCellType


data_frames <- list(r1b1roi1, r1b1roi2, r1c1roi1, r1c1roi2, r1d1roi1, r1d1roi2)


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
    geom_point(size = 1) +  # Adjust point size as needed
    scale_color_manual(values = color_palette) +  # Assign colors based on the palette
    labs(x = "X", y = "Y", title = name) +  # Set title to the name of the data frame
    scale_y_reverse() + # REMOVE THIS LATER if problem is fixed
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank())  # Remove legends
})




# Create a plot to extract the legend
legend_plot <- ggplot(r1b1roi1, aes(x = X_coord, y = Y_coord, color = InferredCellType)) +
  geom_point(size = 1.5) +  # Adjust point size as needed
  scale_color_manual(values = color_palette) +  # Assign colors based on the palette
  theme(legend.position = "right", legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust legend point size

# Extract the legend
legend <- cowplot::get_legend(legend_plot)

# Combine plots and add the legend at the bottom
final_plot <- (plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plots[[6]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

final_plot <- (plots[[1]] | plots[[2]]) /
                 (plots[[3]] | plots[[4]]) /
                 (plots[[5]] | plots[[6]]) +
    plot_layout(guides="collect") &
  theme(legend.position = "right")


print(final_plot)

# Create the scatter plot
ggplot(r1c1roi1, aes(x = X_coord, y = Y_coord, color = InferredCellType)) +
  geom_point(size = 1.5) +  # Adjust point size as needed
  scale_color_manual(values = color_palette) +  # Assign colors based on the palette
  labs(x = "X", y = "Y", title = "Colored Scatter Plot with Strings") +
  theme_minimal() +
  theme(legend.position = "right")  # Adjust legend position as needed























####### Colour by a marker ########
geo <- read.csv("../00_Data/CellShapeData/R1B1ROI1_Cellshape.csv")
averages_tomerge <- averages
averages_tomerge$ID = rownames(averages)
geo$ID = rownames(geo)
merged_geo <- inner_join(geo, averages_tomerge, by = "ID")

# Create the scatter plot with color based on intensity values
ggplot(merged_geo, aes(x = Cell.Center.X, y = Cell.Center.Y, color = cd45ro)) +  # Replace IntensityColumn with the actual column name
  geom_point(size = 1) +  # Adjust point size as needed
  scale_color_gradient(low = "blue", high = "red") +  # Adjust the gradient colors as needed
  labs(x = "X", y = "Y", title = "Colored Scatter Plot Based on Intensity Values") +
  theme_minimal() +
  theme(legend.position = "right")  # Adjust legend position as needed
 


colnames(averages)



