# Assuming you have your data frame "geo" with coordinates and a column "InferredCellType"
# Replace these with your actual data
geo <- read.csv("../00_Data/CellShapeData/R1B1ROI1_Cellshape.csv")
# Load the ggplot2 package
library(ggplot2)
library(dplyr)

# Define the color palette for the cell types
color_palette <- c("Epithelial cells"= "red",     
                   "Extracellular matrix"= "green",
                   "Lymphatic vessel"= "blue",
                   "Mixed"= "yellow",
                   "Myocytes"= "maroon",
                   "Other"= "cyan",
                   "T-cells"= "magenta",
                   "T-cells CD4+ (helper)"= "darkgreen",
                   "T-cells Naive"= "navy",
                   "T-cells Mixed"= "gray" )

grey_palette <- c("Epithelial cells"= "magenta",     
                  "Extracellular matrix"= "grey",
                  "Lymphatic vessel"= "grey",
                  "Mixed"= "grey",
                  "Myocytes"= "grey",
                  "Other"= "grey",
                  "T-cells"= "grey",
                  "T-cells CD4+ (helper)"= "grey",
                  "T-cells Naive"= "grey",
                  "T-cells Mixed"= "grey" )

geo$InferredCellType = celltypes$InferredCellType
# Create the scatter plot
ggplot(geo, aes(x = Cell.Center.X, y = Cell.Center.Y, color = InferredCellType)) +
  geom_point(size = 1) +  # Adjust point size as needed
  scale_color_manual(values = grey_palette) +  # Assign colors based on the palette
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
ggplot(merged_geo, aes(x = Cell.Center.X, y = Cell.Center.Y, color = cytokeratin)) +  # Replace IntensityColumn with the actual column name
  geom_point(size = 1) +  # Adjust point size as needed
  scale_color_gradient(low = "blue", high = "red") +  # Adjust the gradient colors as needed
  labs(x = "X", y = "Y", title = "Colored Scatter Plot Based on Intensity Values") +
  theme_minimal() +
  theme(legend.position = "right")  # Adjust legend position as needed
 






