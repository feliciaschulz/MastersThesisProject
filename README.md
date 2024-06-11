# README Bioinformatics Thesis Project

## Introduction to this repository

### Folder structure of this repository

## Abstract


## 1. Segmentation

### 1.1 Stardist

### 1.2 Cellpose

## 2. Create intensity frames
In this script, raw (stitched) images and their respective segmentation masks are loaded. The script reads the images into the notebook and their dimensions are checked. 
If the dimensions are correct, the images are measured and their data is saved to data frames. The intensity values for each cell (min, mean, max) are saved, as well as cell area (in pixels) and their centroids. 

If the dimensions of the images and the mask are not identical, the cells cannot be measured. 
For the data in this project, this was the case a few times. However, there was only padding missing in the bottoms of some images. This means horizontally at the bottom of the images, rows of black pixels had to be added.
This script allows for that to happen, if necessary.
If the dimensions are not the same for any other reason, the padding code in this script is not appropriate for fixing the problem.

After creating the data frames, they are saved to the data folder.

#### Usage
If you don't have jupyter notebook yet, install it in the terminal like this:
```bash
pip install jupyter notebook
```
To run the notebook, run this command in the directory in which you have the scripts (or any earlier directory):
```bash
jupyter notebook
```
Jupyter notebook will be opened in your browser. You can navigate to your script file and open it.

This notebook can be run as-is, cell-by-cell. 

#### Required packages
The code was written in Python v. 3.10.12. The following packages are required for and automatically loaded in the notebook:
- scikit-image: 0.22.0
- numpy: 1.26.4
- pandas: 2.2.0

#### Input
This script requires a labelled segmentation mask as input, which should be saved in "00_Data/NucleiCytoMasks". It also takes stitched preprocessed multiplexed stained images as input from the folder "00_Data/PP_Images".
All images and masks should be the exact same shape. The input files should be .tif files.

This script is initialised with the variable roi, which decides the name of the sample image to be measured.
In this script, the correct name of the image available as example data is already assigned to the roi variable.

#### Output
This script saves one data frame for each biomarker-stained input image to "00_Data/IntensityDataFrames/<img name>" as csv files.
Make sure to have created the output folder already before saving the images.


## 3. Initial data investigation

## 4. Cell Type Inference


## 5. Dimensionality Reduction


## 6. Cell Visualisation Application
This is a Shiny application coded in R with which a UMAP plot can be created and the result visualised. The app was created out of the realisation that there would be much to benefit from being able to relate the clusters formed in a dimensionality reduction plot back to their original location and properties.

The app works as follows: The data is loaded in and the original image is displayed as a  scatterplot with each cell as a point. 
A umap can be created out of different subsets of the data, selected from a sidebar menu.
The umap plot appears on the interface and the user can hover over
points to see them displayed in the original image.
The user can also select a certain region in the umap plot, and the selected cells will be displayed on the original image. Double-click to remove the selected region.

#### Usage
Open the script 06_CellVisualisationApplication/app.R in RStudio. Run the app by clicking "Run App" in the upper right corner in RStudio. Alternatively, run this command in the RStudio console:
```R
runApp()
```

OBS: Make sure that the "intensities" data frame has the correct path. The variable my_img can be changed to inspect different images, given the format of the intensities data frame is the same.

#### Required packages
The app uses the following packages and versions:
- Shiny: 1.8.1.1
- plotly: 4.10.4
- bslib: 0.7.0
- Rtsne: 0.7.0 !!
- umap: 0.2.10.0
- dplyr: 1.1.4

Check your package versions in the R console with 
```R
packageVersion("<package name>")
```


#### Known bugs / weaknesses
This application can have long loading times. Depending on the size of the subset selected for the umap, the plot creation can take up to a few minutes. This is solely because the umap() function is computationally complex. 

Hover slowly. Due to the live computations, there will be a slight delay.

When selecting a region, after drawing the circle around the chosen cells, leave the cursor in place and don't keep hovering. Otherwise, the application will detect your movement as hovering again and the selected region will disappear. 
The displaying of the selected region can also take a little while, because identifying the cells and referring to the original data frame with coordinates is computationally complex as well.

#### Input 
A data frame csv file with cells in rows, biomarker intensity data in columns, two columns for x and y coordinates (X_coord, Y_coord), area of cell (area), cell phenotype (InferredCellType). The input data should be accessible from the folder "../00_Data/IntensitiesCelltypes/".

At the moment, the app only works for input data that has the format of the example in this repository.

In this repository, the example data can be found in 00_Data/IntensitiesCelltypes/R1B1ROI1_intensitiesdf.csv.
With this input data, nothing in the code of app.R has to be changed and it can be run as is.

#### Output
When running the app as specified in "Usage", it will automatically open in browser. Follow the above instructions to use the app.

## 7. Distance Analysis
With the script 07_DistanceAnalysis/DistanceMatrix.ipynb, distance violin plots can be made as well as neighbour sankey diagrams.

#### Usage
If you don't have jupyter notebook yet, install it in the terminal like this:
```bash
pip install jupyter notebook
```
To run the notebook, run this command in the directory in which you have the scripts (or any earlier directory):
```bash
jupyter notebook
```
Jupyter notebook will be opened in your browser. You can navigate to your script file and open it.

This notebook can be run as-is, cell-by-cell. 
Small user modifications can be made in the notebook itself if desired. These options are indicated by "### USER-INPUT ###". More about this in the notebook itself.

#### Required packages
The code was written in Python v. 3.10.12. The following packages are required for and automatically loaded in the notebook:
- os (in Python's standard utility packages)
- statistics (in Python's standard utility packages)
- matplotlib: 1.26.4
- scikit-image: 0.22.0
- pandas: v. 2.2.0
- scipy: 1.12.0
- plotly: 5.19.0


#### Input
As input data, this script takes the intensity data frames created in 2. that are saved in 00_Data/IntensitiesCelltypes. They should be csv files.

#### Output
The script creates data frames and plots, none of which are automatically saved, because there are many different user-specific options. 
The plots are created in the notebook and are displayed there as well.
Instructions as to how to save the plots if desired can be found in the notebook itself.


















