#to rm(list=ls())

geo <- read.csv("../00_Data/CellShapeData/R1B1ROI1_Cellshape.csv")
head(geo)
View(geo)

# find image coordinates
max(geo$Cell.Center.X) #9422.321
min(geo$Cell.Center.X) #9.728
max(geo$Cell.Center.Y) #10328.55
min(geo$Cell.Center.Y) #9.152



################
# TRYING SPIAT #
################

library(SPIAT)
head(averages)
averages_matrix = as.matrix(averages)

# Matrix has to be transposed to fit SPIAT object structure
intensity_matrix = matrix(0, nrow = ncol(averages_matrix), ncol = nrow(averages_matrix))
# Calculate the transpose using nested for loops
for (j in 1:nrow(averages_matrix)) {
  for (i in 1:ncol(averages_matrix)) {
    intensity_matrix[i, j] <- averages_matrix[j, i]
  }
}
rownames(intensity_matrix) <- colnames(averages)
colnames(intensity_matrix) <- rownames(averages)

# Define other metadata
cellType <- celltypes$InferredCellType # OBS: Not pathologist-verified phenotype, only inferred cell type
coord_x <- geo$Cell.Center.X
coord_y <- geo$Cell.Center.Y

# SPIAT object
general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix,
                                            phenotypes = cellType,
                                            coord_x = coord_x, coord_y = coord_y)



### Inspect the SpatialExperiment object
colData(general_format_image)[1:5,] # shows phenotype and cell properties
assay(general_format_image)[,1:5 ] # shows intensity level of every marker for every cell
spatialCoords(general_format_image)[1:5, ] # shows cell coordinates
unique(general_format_image$Phenotype) # see which ones there are




##############
# NEXT STEPS #
##############

edited_marker_celltypes <- read.csv("New_Marker_Celltypes_ToEdit.csv", sep=";", header=TRUE)
categories <- c()
names <- c()

# Loop through each row of the data frame
for (i in 1:nrow(edited_marker_celltypes)) {
  # Get the row
  row <- edited_marker_celltypes[i, ]
  # Extract the name and markers
  name <- row$X
  markers <- row[-c(1)]
  
  # Extract the markers with value 1
  positive_markers <- names(markers)[markers == 1]
  print(positive_markers)
  
  # Append to categories and names vectors
  if (length(positive_markers) == 0) {
    categories <- c(categories, "OTHER")
  } else {
    categories <- c(categories, paste(positive_markers, collapse = ","))
  }
  names <- c(names, name)
}
result <- data.frame(names, categories)



# Define cell types
formatted_image <- define_celltypes(
  general_format_image, 
  categories = categories, 
  category_colname = "Phenotype", 
  names = names,
  new_colname = "Cell.Type")

marker_intensity_boxplot(formatted_image, "actin")


# Cell distance analysis
distances <- calculate_pairwise_distances_between_celltypes(
  spe_object = formatted_image, 
  cell_types_of_interest = c("T-cells", "Extracellular matrix", "Myocytes"),
  feature_colname = "Cell.Type")

#### not working, changing averages data frame:
thresholded_averages <- averages
# for (i in 1:nrow(averages)) {
#   
#   for (j in 1:ncol(averages)) {
#     
#     marker = colnames(averages[j])
#     
#     if (averages[i,j] < thresholds_clean[marker,"DispRangeMin"]) {
#       thresholded_averages[i,j] <- 0
#     } 
#   }
# }

intensity_matrix_th = matrix(0, nrow = ncol(averages_matrix), ncol = nrow(averages_matrix))
# Calculate the transpose using nested for loops
for (i in 1:ncol(averages_matrix)) {
  for (j in 1:nrow(averages_matrix)) {
    intensity_matrix_th[i, j] <- averages_matrix[j, i]
  }
}
rownames(intensity_matrix_th) <- colnames(averages)
colnames(intensity_matrix_th) <- rownames(averages)

# Define other metadata
cellType <- celltypes$InferredCellType # OBS: Not pathologist-verified phenotype, only inferred cell type
coord_x <- geo$Cell.Center.X
coord_y <- geo$Cell.Center.Y

# SPIAT object
general_format_image <- format_image_to_spe(format = "general",
                                            intensity_matrix = intensity_matrix_th,
                                            phenotypes = cellType,
                                            coord_x = coord_x, coord_y = coord_y)


# Define cell types
my_formatted_image <- define_celltypes(
  general_format_image, 
  categories = categories, 
  category_colname = "Phenotype", 
  names = names,
  new_colname = "Cell.Type") # still not working


data("simulated_image")

# define cell types
formatted_image <- define_celltypes(
  simulated_image, 
  categories = c("Tumour_marker","Immune_marker1,Immune_marker2", 
                 "Immune_marker1,Immune_marker3", 
                 "Immune_marker1,Immune_marker2,Immune_marker4", "OTHER"), 
  category_colname = "Phenotype", 
  names = c("Tumour", "Immune1", "Immune2", "Immune3", "Others"),
  new_colname = "Cell.Type")

