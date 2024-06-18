##############################################
# DIMENSIONALITY REDUCTION PLOTS - dimensionality_reduction.R
# 
# This script takes files with the inferred cell types of cells as well as their
# coordinates and creates dimensionality reduction plots.
# The user can choose between the dimensionality reduction methods ("t-SNE" or "UMAP").
# They can also choose which inferred cell types to include in the plots, or rather, 
# they get four options: 1. Remove "mixed" cells. 2. Remove "other" cells. 
# 3. Remove both "other" and "mixed" cells. 4. Only keep "mixed" cells.
# After creating the dimensionality reduction data frames, they are also plotted.
# Here, if the user has chosen the UMAP method, they can also choose whether to
# plot it normally, colouring the points according to their inferred cell types,
# or whether to colour them according to their expression levels of a biomarker.
#
# Reminder of the available biomarkers that can be used when plotting expression levels:
#  "actin", "cd3", "cd4", "cd45", "cd45ro", "collageni", "cytokeratin",
#  "fibulin2", "lyve1", "podoplanin", "cd38", "cd138"
#
# Usage: 
# Run the script line-by-line in RStudio. 
# If you want to change the input image, change the string assigned to
# the variable "my_img" below. Also change desired parameters where it says "### USER ###".
#
# Required packages & versions: 
#   - tidyverse: 2.0.0
#   - Rtsne: 0.17
#   - umap: 0.2.10.0
#   - spatialEco: 2.0.2
# 
# Input: 
#     Raw data frame with cells in rows, biomarker intensity expression data in columns, area,
#     x and y coordinates, and inferred cell type in last column from "../00_Data/IntensitiesCelltypes/" 
#     called "*_intensitiesdf.csv"
#
#
# Output: Nothing, depending on settings chosen by user, dimensionality reduction plots can be inspected. 
# They can also be saved by clicking "Export" in RStudio above the plot.
###############################################

# Load packages
library(tidyverse)
library(Rtsne)
library(umap)
library(spatialEco)


# Set seed for consistency
set.seed(20)

# Load in data
img_name = "R1B1ROI1"
df = read.csv(paste0("../00_Data/IntensitiesCelltypes/", img_name, "_intensitiesdf.csv"))
df$X = NULL

### USER ###

# Creating dimensionality reduction:
method = "UMAP" # method = "UMAP" OR method = "t-SNE
exclude_mixed=FALSE # Put this to TRUE if you don't want "mixed" cells in the plot
exclude_other=TRUE # Put this to TRUE if you don't want "other" cells in the plot
only_mixed=FALSE # Put this to TRUE if you ONLY want "mixed" cells in the plot

# Plotting:
umap_expression_plot = FALSE # This works only with umap, creates a plot with only one marker and its expression levels
marker = "cytokeratin" # This will only be used if umap_expression_plot == TRUE



################ Creating functions

#### t-sne ####
tsne_function <- function(df, pca=TRUE, perplexity=30, max_iter=1000, epoch=100, exclude_mixed=FALSE, exclude_other=FALSE, only_mixed=FALSE) {  
  
  if (exclude_mixed==TRUE) { 
    df <- df[df$InferredCellType != "Mixed", ] # Remove all "mixed" cells from df
  }
  
  if (exclude_other==TRUE) {
    df <- df[df$InferredCellType != "Other", ] # Remove all "other" cells from df
  }
  
  if (only_mixed==TRUE) {
    df <- df[df$InferredCellType == "Mixed", ] # Remove all NOT "mixed" cells from df
  }
  
  cellType = df$InferredCellType # save cell types
  df = df[1:12] # remove non-intensity columns
  
  set.seed(20) # make sure its always the same
  # run tsne
  tsne_out <- Rtsne(df, 
                    pca=pca,
                    perplexity=perplexity,
                    k=2,
                    max_iter=max_iter,
                    epoch=epoch)
  
  # save output as data frame
  Y = as.data.frame(tsne_out$Y)
  
  # return result
  return(list(Y, cellType))
}


#### UMAP ####
umap_function <- function(df, exclude_mixed=TRUE, exclude_other=TRUE, only_mixed=FALSE) {
  set.seed(20)
  if (exclude_mixed==TRUE) {
    df <- df[df$InferredCellType != "Mixed", ] # Remove all "mixed" cells from df
  }
  
  if (exclude_other==TRUE) {
    df <- df[df$InferredCellType != "Other", ] # Remove all "other" cells from df
  }
  
  if (only_mixed==TRUE) {
    df <- df[df$InferredCellType == "Mixed", ] # Remove all NOT "mixed" cells from df
  }
  
  umap_celltypes <- df$InferredCellType # save cell types
  umap_intensities <- df[1:12] # save intensity values in case user wants to plot marker expression
  
  df <- df[1:12] # remove all non-intensity columns
  # run umap
  umap_out <- umap(df)
  
  # return
  return(list(umap_out, umap_celltypes, umap_intensities))
}


#### umap marker expression #####
umap_plot_onlyone <- function(umap_out, umap_averages, marker) {
  
  # umap output
  layout <- umap_out$layout
  
  # initialise plot with title, colour palette, font sizes
  main=paste(marker)
  palette = colorRampPalette(c('grey','blue'), bias=1)
  rev_palette =  colorRampPalette(c('blue','grey'), bias=1)
  pad=0.1 
  cex=0.6 # point size, can make smaller, originally 0.6
  pch=19
  add=FALSE
  cex.main=4
  
  # for legend: break data up into quantiles (10 parts) so that the legend
  # is not linear but represents density of the data
  quantiles <- quantile(umap_averages[[marker]], probs = seq(0, 1, by = 0.1), na.rm = TRUE, type = 6)
  breaks <- unique(round(unname(quantiles)))
  # also create colour palette for the quantile breaks
  umap_averages$Col <- palette(length(breaks))[as.numeric(cut(umap_averages[[marker]], breaks = breaks))]
  
  # define plot layout (x and y limits / plot size)
  xylim <- range(layout)
  xylim_small <- xylim + c(-1, 1) # for padding
  xylim_big <- xylim + c(-5, 1)
  
  # create plots with layout defined above
  par(mar=c(0.2,0.7,1.2,0.7), ps=10) # this is for formatting; mar is margin and ps is point size (font size)
  plot(xylim_big, xylim_small, type="n", axes=F, frame=F)
  rect(xylim_small[1], xylim_small[1], xylim_small[2], xylim_small[2], border="#aaaaaa", lwd=0.25)  # rectangle lining plot
  
  # plot points
  points(layout[,1], layout[,2], col=umap_averages$Col,
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main, line=-1.2) # title 
  
  # create legend
  legend_text <- round(breaks) # quantile value labels for legend
  legend_image <- as.raster(matrix(rev_palette(length(breaks)), ncol=1)) # define raster palette
  text(x=xylim_big[1]+1.8, y = seq(-7,9,l=length(breaks)), labels = legend_text, cex=1.5) # position labels
  rasterImage(legend_image, xylim_big[1],-7, xylim_big[1]+0.5,9) # create raster palette
  
}



########### Running functions with user-defined parameters

#### run dim. red. functions ####
if (method == "t-SNE") {
  # run tsne
  tsne_list <- tsne_function(df, exclude_mixed=exclude_mixed, exclude_other=exclude_other, only_mixed=only_mixed)
  Y <- tsne_list[[1]]
  cellType <- tsne_list[[2]]
} else if (method == "UMAP") {
  # run umap
  umap_list <- umap_function(df, exclude_mixed=exclude_mixed, exclude_other=exclude_other, only_mixed=only_mixed)
  layout <- umap_list[[1]]$layout
  cellType <- umap_list[[2]]
  umap_intensities <- umap_list[[3]]
  Y <- as.data.frame(layout)
  Y$cellType <- cellType
}



#### plot ####
if (method == "tsne" | umap_expression_plot == FALSE) {
  # plot dimred
  ggplot(Y,
         aes(x = V1, y = V2, col = cellType)) +
    
    geom_point() +
    
    labs(x = "UMAP - 1",
         y = "UMAP - 2") +
    
    ggtitle(paste0(method, " for ", img_name)) +
    
    theme_classic()
} else {
  # plot umap specific for marker
  umap_plot_onlyone(umap_out, umap_intensities, marker)
}

