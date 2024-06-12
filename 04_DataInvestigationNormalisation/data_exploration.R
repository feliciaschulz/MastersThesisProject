###########################################
# EXPLORE DATA AND NORMALISE
#
# This script takes the intensities data frame as input and lets the
# user run a variety of functions to inspect the data, such as marker expression levels.
# In addition to that, multiple normalisation methods are available as functions.
# These can again be plotted with the plotting functions and the results can be compared.
#
# Usage: Due to the different functionalities of this script, it is up to the
# user to decide which functions to run and which markers to inspect and plot.
# Therefore, this script is best run line-by-line (rather than sourcing).
# The user can also choose the image to be investigated. Here, it is already set to
# R1B1ROI1 because that is the example image provided in this repository.
#
# Required packages & versions: 
#   - tidyverse: 2.0.0
#   - DescTools: 0.99.54
#   - ggplot2: 3.5.0
#
# Input: An intensity data frame from "00_Data/IntensitiesCelltypes/"
#
# Output: Depends on each function. Some output text in the console, some plots,
# some return data frames.
#
# Functions:
#   - outlierKD(): Takes intensities df and a marker name as input, checks for outliers.
#     Plots a histogram and a box plot of the distribution for this marker, with 
#     and without outliers. Also prints outlier statistics to console. Function has
#     been adapted from https://www.r-bloggers.com/2016/04/identify-describe-plot-and-remove-the-outliers-from-the-dataset/
#
#   - plot_marker(): Plots just simple histogram of marker expression of a specified marker
#   - minmax_func(): Function for creating min-max normalisation, returns normalised df
#   - winsorize_func(): Function for creating winsorization, returns winsorized df
#   - standardize_func(): Function for standardization, returns stuandardized df
#######################################

# Load packages
library(tidyverse)
library(DescTools) # for winsorize()
library(ggplot2)

# Load data
my_img = "R1B1ROI1"
intensities_df <- read.csv(paste0("../00_Data/IntensitiesCelltypes/", my_img, "_intensitiesdf.csv"))
intensities_df$X <- NULL
intensities <- intensities_df[1:12] # remove area data, coords, inf.c.t.

# Reminder of available markers to choose from:
# names(intensities)
# "actin"       "cd3"         "cd4"         "cd45"        "cd45ro"      "collageni"  
# "cytokeratin" "fibulin2"    "lyve1"       "podoplanin"  "cd38" 


#####################
# outlier detection #
#####################

# code from https://www.r-bloggers.com/2016/04/identify-describe-plot-and-remove-the-outliers-from-the-dataset/
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", breaks=50, xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", breaks=50, xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "\n")
  cat("Proportion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "\n")
  cat("Mean of the outliers:", round(mo, 2), "\n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "\n")
  cat("Mean if we remove outliers:", round(m2, 2), "\n")
}

outlierKD(intensities, fibulin2)



#######################
# marker distribution #
#######################

plot_marker <- function(df, marker) {
  # marker = df[[marker]]
  # Create a histogram
  markercol = df[[marker]]
  plot <- ggplot(df, aes(x = markercol)) +
    geom_histogram(binwidth = max(markercol)/100, fill = "skyblue", color = "black") +
    labs(title = paste("Intensity Distribution", marker, sep=" "),
         x = "Intensity",
         y = "Frequency")
  
  return(plot)
}

marker_plot <- plot_marker(intensities, "fibulin2")
marker_plot



#################
# NORMALISATION #
#################



# min-max normalisation
minmax_func <- function(intensities) {
  intensities_minmax <- data.frame(matrix(nrow = length(rownames(intensities)), ncol = length(colnames(intensities))))
  rownames(intensities_minmax) <- rownames(intensities)
  colnames(intensities_minmax) <- colnames(intensities)
  for (j in 1:ncol(intensities)) {
    biggest = max(intensities[[j]])
    smallest = min(intensities[[j]])
    
    for (i in 1:nrow(intensities)) {
      intensity = intensities[i,j]
      new_intensity = (intensity-smallest)/(biggest-smallest)
      intensities_minmax[i,j] = new_intensity
    }
  }
  return(intensities_minmax)
}
intensities_minmax = minmax_func(intensities)


# winsorization
winsorize_func <- function(intensities) {
  intensities_winsorized <- as.data.frame(matrix(NA, nrow = nrow(intensities), ncol = ncol(intensities)))
  rownames(intensities_winsorized) <- rownames(intensities)
  colnames(intensities_winsorized) <- colnames(intensities)
  
  for (i in colnames(intensities)) {
    outlier <- boxplot.stats(intensities[[i]])$out # outliers identified with the interquartile range (IQR) criterion
    m = mean(intensities[[i]])
    num_outliers <- length(outlier)
    percent_outliers <- round(num_outliers / sum(!is.na(intensities[[i]]))*100, 1)
    percent_outliers_low <- round(sum(outlier < m) / sum(!is.na(intensities[[i]]))*100, 1)
    percent_outliers_high <- round(sum(outlier > m) / sum(!is.na(intensities[[i]]))*100, 1)
    print(i)
    print(paste("There are ", num_outliers, "outliers in total, making up", percent_outliers, "% of the total values."))
    print(paste(percent_outliers_low, "% are below the mean and", percent_outliers_high, "% are above the mean."))
    print(paste("Mean",m))
    print(paste("Num of outliers:",sum(outlier > m)))
    print(paste("Outlier intensity value after winsorization:", min(outlier[outlier>m])))
    cat("\n")
    
    
    intensities_winsorized[[i]] <- DescTools::Winsorize(intensities[[i]], probs = c((percent_outliers_low/100) , 1-(percent_outliers_high/100)))
  }
  return(intensities_winsorized)
}


# standardization
standardize_func <- function(intensities) {
  intensities_scaled = as.data.frame(scale(intensities))
  return(intensities_scaled)
}



### Saving final normalised data chosen in this pipeline
# Chosen methods: winsorizationi & scaling
int_wins <- winsorize_func(intensities)
int_wins_scal <- standardize_func(int_wins)
# Add cell type and coords back on
cols_to_add <- intensities_df[, c("area", "X_coord", "Y_coord", "InferredCellType")]
intensities_norm <- cbind(int_wins_scal, cols_to_add)
write.csv(intensities_norm, paste0("../00_Data/Normalised/", my_img, "_norm.csv"))



