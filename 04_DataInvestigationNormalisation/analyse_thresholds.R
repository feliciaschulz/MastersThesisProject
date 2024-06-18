######################################
# ANALYSING THRESHOLDS
#
# In this script, the intensity thresholds are investigated
# and combined into one data frame for comparison.
# The data is loaded, cleaned, merged and filtered.
#
# Usage: This script can be loaded and run as it is.
# !!BUT!! seeing as the purpose of this script is to summarise all threshold data 
# frames, it requires all threshold files as input. However, only the first one
# is provided as example input data here, therefore, this script cannot be run.
# 
# Input: Thresholds data frames from "../00_Data/Thresholds"
#
# Output: One data frame including all thresholds values for all images
#
# Functions: 
#   - rename_markers() for cleaning data and renaming rows
#   - rename_columns() for renaming columns
#######################################






path = "/Users/fschulz/Documents/Uni/LundUniversity/Year2/Thesis/00_Data/Thresholds"
th_r1b1roi1 = read.csv(paste0(path, "/R1B1ROI1", "_Thresholds.csv"))
th_r1b1roi2 = read.csv(paste0(path, "/R1B1ROI2", "_Thresholds.csv"))
th_r1c1roi1 = read.csv(paste0(path, "/R1C1ROI1", "_Thresholds.csv"))
th_r1c1roi2 = read.csv(paste0(path, "/R1C1ROI2", "_Thresholds.csv"))
th_r1d1roi1 = read.csv(paste0(path, "/R1D1ROI1", "_Thresholds.csv"))
th_r1d1roi2 = read.csv(paste0(path, "/R1D1ROI2", "_Thresholds.csv"))

# Change threshold marker names so that I can loop through the thresholds easily
rename_markers <- function(th) {
  th <- subset(th, !grepl("DAPI", Name) & !grepl("background", Name) & !grepl("AF ", Name)) # remove DAPIs, bleached images and autofluorescence
  for (i in seq_len(nrow(th))) {
    words <- strsplit(th$Name[i], " ")[[1]]   # Split the name into words as a list
    if (length(words)>1) {
      words <- words[-length(words)]   # Remove the last word
    }
    words <- tolower(words)   # Convert all words to lowercase
    new_name <- paste(words, collapse = "")   # Concatenate remaining words
    
    # Assign the new name to the corresponding row in the dataframe
    th$Name[i] <- new_name
  }
  return(th)
}

th_r1b1roi1 <- rename_markers(th_r1b1roi1)
th_r1b1roi2 <- rename_markers(th_r1b1roi2)
th_r1c1roi1 <- rename_markers(th_r1c1roi1)
th_r1c1roi2 <- rename_markers(th_r1c1roi2)
th_r1d1roi1 <- rename_markers(th_r1d1roi1)
th_r1d1roi2 <- rename_markers(th_r1d1roi2)

# Function to rename columns and keep only the relevant columns
rename_columns <- function(df, suffix) {
  df <- df[, c("Name", "RangeMin", "RangeMax")]
  names(df)[2:3] <- paste0(c("RangeMin_", "RangeMax_"), suffix)
  return(df)
}


# Rename columns
th_r1b1roi1 <- rename_columns(th_r1b1roi1, "R1B1ROI1")
th_r1b1roi2 <- rename_columns(th_r1b1roi2, "R1B1ROI2")
th_r1c1roi1 <- rename_columns(th_r1c1roi1, "R1C1ROI1")
th_r1c1roi2 <- rename_columns(th_r1c1roi2, "R1C1ROI2")
th_r1d1roi1 <- rename_columns(th_r1d1roi1, "R1D1ROI1")
th_r1d1roi2 <- rename_columns(th_r1d1roi2, "R1D1ROI2")

# Merge data frames
merged_thresholds <- Reduce(function(x, y) merge(x, y, by = "Name", all = TRUE), 
                    list(th_r1b1roi1, th_r1b1roi2, th_r1c1roi1, th_r1c1roi2, th_r1d1roi1, th_r1d1roi2))

# Filtering for only the right markers
filtered_th <- merged_thresholds[merged_thresholds$Name %in% markers, ]


# Select columns with "Min" in their names
min_columns <- filtered_th[, grep("Min", names(filtered_th))]

# Compute the mean and standard deviation for each row
mean_values <- apply(min_columns, 1, mean, na.rm = TRUE)
sd_values <- apply(min_columns, 1, sd, na.rm = TRUE)

# Add the computed values to the data frame
filtered_th$Mean_Min <- mean_values
filtered_th$SD_Min <- sd_values

# View thresholds
View(filtered_th)

# Save data frame as csv file
write.csv(filtered_th, paste0(path, "/filtered_thresholds_all.csv"))



