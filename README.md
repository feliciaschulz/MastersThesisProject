# README Bioinformatics Thesis Project

## Introduction to this repository

### Folder structure of this repository

## Abstract


## 1. Segmentation

### 1.1 Stardist

### 1.2 Cellpose

## 2. Initial data investigation

## 3. Cell Type Inference

## 4. Dimensionality Reduction

## 5. Cell Visualisation Application
This is a Shiny application coded in R with which a UMAP plot can be created and the result visualised. The app was created out of the realisation that there would be much to benefit from being able to relate the clusters formed in a dimensionality reduction plot back to their original location and properties.

The app works as follows: The data is loaded in and the original image is displayed as a  scatterplot with each cell as a point. 
A umap can be created out of different subsets of the data, selected from a sidebar menu.
The umap plot appears on the interface and the user can hover over
points to see them displayed in the original image.
The user can also select a certain region in the umap plot, and the selected cells will be displayed on the original image. Double-click to remove the selected region.

#### Usage
Open the script 05_CellVisualisationApplication/app.R in RStudio. Run the app by clicking "Run App" in the upper right corner in RStudio. Alternatively, run this command in the RStudio console:
```R
runApp()
```

OBS: Make sure that the "intensities" data frame has the correct path. The variable my_img can be changed to inspect different images, given the format of the intensities data frame is the same.

#### Known bugs / weaknesses
This application can have long loading times. Depending on the size of the subset selected for the umap, the plot creation can take up to a few minutes. This is solely because the umap() function is computationally complex. 

Hover slowly. Due to the live computations, there will be a slight delay.

When selecting a region, after drawing the circle around the chosen cells, leave the cursor in place and don't keep hovering. Otherwise, the application will detect your movement as hovering again and the selected region will disappear. 
The displaying of the selected region can also take a little while, because identifying the cells and referring to the original data frame with coordinates is computationally complex as well.

#### Input 
A data frame with cells in rows, biomarker intensity data in columns, two columns for x and y coordinates (X_coord, Y_coord), area of cell (area), cell phenotype (InferredCellType). The input data should be accessible from the folder "../00_Data/IntensitiesCelltypes/".

At the moment, the app only works for input data that has the format of the example in this repository.

In this repository, the example data can be found in 00_Data/IntensitiesCelltypes/R1B1ROI1_intensitiesdf.csv.
With this input data, nothing in the code of app.R has to be changed and it can be run as is.

#### Output
When running the app as specified in "Usage", it will automatically open in browser. Follow the above instructions to use the app.

## 6. Distance Analysis



















