# README Bioinformatics Thesis Project

## Introduction to this repository

### Folder structure of this repository

## Abstract
The aim of this research project was to create a bioinformatic image analysis pipeline for the investigation of lymphocytes in healthy human lungs. 
The data consisted of three formalin-fixed paraffin-embedded (FFPE) tissue samples, divided into six regions of interest (ROIs) which were imaged using a multiplexed cyclic fluorescence-based antibody staining approach with 26 antibodies. 

State-of-the-art downstream analysis methods and tools were applied to stitch the images, segment the cells, perform dimensionality reduction, and investigate the geographical properties of the data.

A manual gating approach was used to infer cell types based on marker intensity expression levels and binary thresholds.
The cell type inference was tested using multiple statistical plotting methods and proved to be overall successful.
A cell visualisation app was created to provide an interface for further analysing dimensionality reduction plots, with which cell type misclassifications could be examined and explained.

The locations of T-cells in the tissues were investigated using distance analysis methods and no particular trends were found regarding which structural regions of the tissue the cells were located in.

Regardless, this methodology was a successful approach to analysing multiplexed imaging data.
In the future, it can be applied to further classifying lymphocytes in healthy lungs, for example to distinguish circulating cells from tissue-resident memory cells.
It can also be used to further investigate diseased tissue.


## 1. Segmentation

### 1.1 StarDist
With the tool StarDist, a cell segmentation can be predicted on image .tif files. Here, a pretrained StarDist model 2D_versatile_fluo was used to make the predictions.

Find the StarDist paper here:
Schmidt, U., Weigert, M., Broaddus, C., & Myers, G. (2018). Cell detection with star-convex polygons. In Medical Image Computing and Computer Assisted Intervention–MICCAI 2018: 21st International Conference, Granada, Spain, September 16-20, 2018, Proceedings, Part II 11 (pp. 265-273). Springer International Publishing.

Find the GitHub repository here: https://github.com/stardist/stardist

The script used here called **StardistPretrained2DPrediction.ipynb** was adapted from the example usage scripts on StarDist's GitHub.
In a Jupyter Notebook, images and the segmentation model are loaded. The prediction is then carried out, the image with the segmentation mask is plotted and the results can be saved to folder. The script can be found in 01_Segmentation.

#### Usage
An environment file can also be found in 01_Segmentation to facilitate usage.
This is how you can use the stardist_env.yml file to create the conda environment, and have Jupyter Notebook run on the specific environment:

```bash
# You should be in the folder 01_Segmentation
# Creating environment from yml file
conda env create -f stardist_env.yml

# Adding kernel to jupyter
conda install -c anaconda ipykernel
python -m ipykernel install --user --name stardist_env --display-name "Python (stardist_env)"

# Opening jupyter notebook
jupyter notebook
# After having opened a jupyter notebook, choose "Python (stardist_env)" as the kernel
```

Once all the required packages have been installed and can be imported, the notebook can be run cell-by-cell. Make sure the input data is in the folder "prediction_on_images_stardist2d" and that the images have the extension .tif. 

#### Required packages
StarDist was used with version 0.8.5. All other required packages can be found in the stardist_env.yml file.

#### Input
Images in .tif file format from the folder "prediction_on_images_stardist2d" that is in the same folder as this Notebook.

#### Output
A segmentation mask called "labels*.tif", the original image called "image*.tif" and a zip file containing segmented cell ROIs called "rois*.zip", saved to the same folder in which the Notebook is. Images are also plotted within the notebook itself.

### 1.2 Cellpose
With Cellpose, image segmentations can also be created either using a pretrained model from the Cellpose model zoo, or by refining the prediction and retraining one's own network with a human-in-the-loop approach using the Cellpose GUI (with Cellpose 2.0).

Find the general Cellpose paper here:
Stringer, C., Wang, T., Michaelos, M., & Pachitariu, M. (2021). Cellpose: a generalist algorithm for cellular segmentation. Nature methods, 18(1), 100-106.

Find the trainable GUI Cellpose paper here:
Pachitariu, M., & Stringer, C. (2022). Cellpose 2.0: how to train your own model. Nature methods, 19(12), 1634-1641.

Find the Cellpose GitHub repository here: https://github.com/MouseLand/cellpose

#### Usage
An environment file can also be found in 01_Segmentation to facilitate usage.
This is how you can use the cellpose_env.yml file to create the conda environment, and have Jupyter Notebook run on the specific environment:

```bash
# You should be in the folder 01_Segmentation
# Creating environment from yml file
conda env create -f cellpose_env.yml
conda activate cellpose_env
```

Then, you can run the Cellpose GUI with this command:
```bash
python -m cellpose
```
Find out more about how to use the Cellpose GUI on the GitHub repository linked above.

#### Required packages
Cellpose was used with version 2.2.3. All other required packages can be found in the cellpose_env.yml file.

#### Input
The GUI has a variety of input options. Most importantly, you can input a .tif image, either by drag-and-drop or by choosing the import option via the menu. Find out more about how to use the Cellpose GUI on the GitHub repository linked above.

#### Output
The GUI has a variety of output options. Most importantly, you can save a .tif segmentation mask or a .zip file with all of the cell ROIs. Find out more about how to use the Cellpose GUI on the GitHub repository linked above.


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





## 3. Cell Type Inference
create_dataDONE, check_celltypes


### 3.1 Create cell type data frame
The next step in this analysis is the computational cell type assignment. For this purpose, the script **create_data.R** in 03_CellTypeInference was made. This script takes cell intensity data frames as input csv files and predicts cell types for each cell based on its intensity values for different markers. 

It uses an input csv file with thresholds for each marker that decide at which intensity value the signal is seen as true and positive. Based on that, it creates a binary data frame with positive / negative expressions for each cell for each marker. Then, a binary marker expression profiles file is used as input to determine which cell type profile fits each cell. The final inferred cell type is then written to a data frame along with area, centroids and intensity information. It is then saved.

#### Usage: 
Run the script either line-by-line in RStudio or source it (top right corner button in RStudio). If you want to change the input image, change the string assigned to the variable "my_img" below.

#### Input: 
- Data frames with cells in rows, biomarker intensity data in columns, two columns for x and y coordinates (X_coord, Y_coord), area of cell (area) from folder "../00_Data/IntensityDataFrames/". One data frame for each marker.
- Thresholds csv file with biomarkers as rows and threshold minimum values as columns from "../00_Data/Thresholds/"
- Marker Expression Profiles csv file with cell types as rows, biomarkers as columns and binary values indicating positive or negative expression from "../00_Data/"

#### Output: 
Two csv files in the folder "../00_Data/IntensitiesCelltypes/"
- ..._intensitiesdf.csv: intensity data as well as centroids, area, inferred cell type
- ..._celltypesdf.csv: binary marker expression data as well as inferred cell type



!! try out with cleared environment to see if any packages are necessary

## 4. Data investigation & Normalisation
data_exploration.RDONE, analyse_thresholds.RDONE, count_instances.RDONE, heatmaps, violin

### 4.1 Analysing thresholds
The intensity thresholds that were chosen for each marker as indicators of positive marker expression vs negative marker expression (true signal vs noise) are saved as csv files in "../00_Data/Thresholds".

These data frames can be summarised and cleaned up with **analyse_thresholds.R**, which can be found in 04_DataInvestigationNormalisation.

#### Usage:
Open the script in RStudio and run, either line-by-line or by using the "Source" button.
BUT: Seeing as the purpose of this script is to summarise all threshold data frames, it requires all threshold files as input. However, only the first one is provided as example input data here, therefore, this script cannot be run this way.

#### Input: 
Thresholds data frames from "../00_Data/Thresholds"

#### Output: 
One data frame including all thresholds values for all images called filtered_thresholds_all.csv, saved in "00_Data/Thresholds"

### 4.2 Data Exploration & Normalisation
This script, **data_exploration.R** in "04_DataInvestigationNormalisation" takes the intensities data frame as input and lets the user run a variety of functions to inspect the data, such as marker expression levels.
In addition to that, multiple normalisation methods are available as functions.
These can again be plotted with the plotting functions and the results can be compared.

#### Usage: 
Due to the different functionalities of this script, it is up to the user to decide which functions to run and which markers to inspect and plot.
Therefore, this script is best run line-by-line (rather than sourcing).
The user can also choose the image to be investigated. Here, it is already set to R1B1ROI1 because that is the example image provided in this repository.

#### Required packages & versions: 
- tidyverse: 2.0.0
- DescTools: 0.99.54
- ggplot2: 3.5.0

#### Input: 
An intensity data frame from "00_Data/IntensitiesCelltypes/"

#### Output: 
Depends on each function. Some output text in the console, some plots, some return data frames.
If the script is run as-is, for example sourced, the loaded intensity data frame will be winsorized, scaled and then saved to "../00_Data/Normalised/".

### 4.3 Count cell type abundances
This script called **count_instances.R** takes celltypes data frames with binary expression for each marker and creates one data frame with total cell type counts.
The output data frame also shows the biomarker counts that contribute to each cell type.

#### Usage: 
Run the script either line-by-line in RStudio or source it (top right corner button in RStudio). If you want to change the input image, change the string assigned to the variable "my_img".

#### Input: 
Data frame .csv file with cells in rows, biomarker binary expression data in columns, inferred cell type in last column from "../00_Data/IntensitiesCelltypes/" called "*_celltypesdf.csv".
The script can also be used with multiple input data frames that then get combined. To do that, change the variable "onlyOneDF" in the script to FALSE and change the paths to point to the correct input files.

#### Output: 
One csv file in the folder "00_Data/" called "*celltypeCounts.csv". The data frame also opens automatically in the script when it's created.

### 4.4 Heatmaps
The script **heatmap_celltypeintensities.R** takes intensities data frames and creates heatmaps showing intensities, organised by markers horizontally and cell types vertically.
Given the user has already created normalised intensity data frames with data_exploration.R, they can also change the variable "norm" to TRUE and thereby use the normalised data as input instead.
Additionally, the user can also choose to only include lymphocytes in the heatmap. This will make it easier to see the individual cells of the lymphocytes.

#### Usage: 
Run the script line-by-line in RStudio. If you want to change the input image, change the string assigned to the variable "my_img" below. If you want to use the normalised data as input, change the variable "norm" to TRUE. If you want to only plot lymphocytes in the heatmap, change plot_lymph_separate to TRUE.

#### Required packages: 
- ComplexHeatmap: 2.18.0

#### Input: 
- Raw data frame with cells in rows, biomarker intensity expression data in columns, area, x and y coordinates, and inferred cell type in last column from "../00_Data/IntensitiesCelltypes/" called "*_intensitiesdf.csv"
- Normalised data frame with the cells in rows, biomarker intensity expression data in columns, area, x and y coordinates, and inferred cell type in last column from "../00_Data/Normalised/" called "*_norm.csv"

#### Output: 
Nothing, depending on settings chosen by user, heatmaps are shown and can be inspected. Can be saved by clicking "Export" in RStudio above the plot.

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
- Rtsne: 0.7.0
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

When running the application, warnings might appear, they can be ignored.

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


















