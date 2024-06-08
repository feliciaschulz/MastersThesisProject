#rm(list=ls())
library(tidyverse)
library(Rtsne)
library(umap)
library(spatialEco)
library(patchwork)
library(cowplot)

#source("create_data.R")

set.seed(20)

###############
#### t-sne ####
###############

tsne_function <- function(averages, celltypes, area_data, pca=TRUE, perplexity=30, max_iter=1000, epoch=100, just_plot=FALSE, exclude_mixed=FALSE, exclude_other=FALSE, only_mixed=FALSE) {  
  if (just_plot==FALSE) { # Actually computing tsne
    if (!("Area" %in% names(averages))) {
      averages["Area"] = area_data # add area
    }
    
    if (exclude_mixed==TRUE) {
      averages <- averages[celltypes$InferredCellType != "Mixed", ]
      celltypes <- celltypes[celltypes$InferredCellType != "Mixed", ]
    }
    
    if (exclude_other==TRUE) {
      averages <- averages[celltypes$InferredCellType != "Other", ]
      celltypes <- celltypes[celltypes$InferredCellType != "Other", ]
    }
    
    if (only_mixed==TRUE) {
      averages <- averages[celltypes$InferredCellType == "Mixed", ]
      celltypes <- celltypes[celltypes$InferredCellType == "Mixed", ]
    }
    
    set.seed(10) # make sure its always the same
    tsne_out <- Rtsne(averages,
                      pca=pca,
                      perplexity=perplexity,
                      k=2,
                      max_iter=max_iter,
                      epoch=epoch)
    
    Y = as.data.frame(tsne_out$Y)
    cellType = celltypes$InferredCellType
  }
  
  
  # plot tsne
  ggplot(Y,
         aes(x = V1, y = V2, col = cellType)) +
    
    geom_point() +
    
    labs(x = "tSNE - 1",
         y = "tSNE - 2") +
    
    ggtitle("t-SNE for R1B1ROI1 with nuclei and cyto segmentation") +
    
    theme_classic()

}

# run tsne
tsne_function(averages, celltypes, area_data, pca=TRUE, perplexity=30, max_iter=1000, epoch=100, just_plot=FALSE, exclude_mixed=FALSE, exclude_other=FALSE, only_mixed=TRUE)



##############
#### UMAP ####
##############

umap_function <- function(averages, celltypes, exclude_mixed=TRUE, exclude_other=FALSE, only_mixed=FALSE) {
  set.seed(20)
  
  if (exclude_mixed==TRUE) {
    averages <- averages[celltypes$InferredCellType != "Mixed", ]
    celltypes <- celltypes[celltypes$InferredCellType != "Mixed", ]
  }
  
  if (exclude_other==TRUE) {
    averages <- averages[celltypes$InferredCellType != "Other", ]
    celltypes <- celltypes[celltypes$InferredCellType != "Other", ]
  }
  
  if (only_mixed==TRUE) {
    averages <- averages[celltypes$InferredCellType == "Mixed", ]
    celltypes <- celltypes[celltypes$InferredCellType == "Mixed", ]
  }
  print("works")
  umap_out <- umap(averages)
  umap_averages <- averages
  umap_celltypes <- celltypes
  
  return(list(umap_out, umap_averages, umap_celltypes))
}

# just for now #
averages <- r1b1roi1[1:13]
celltypes <- r1b1roi1

umap_list <- umap_function(averages, celltypes, exclude_mixed=FALSE, exclude_other=FALSE, only_mixed=TRUE)
umap_out <- umap_list[[1]]
umap_averages <- umap_list[[2]]
umap_celltypes <- umap_list[[3]]


########
# plot #
########

umap_plot <- function(umap_out, umap_averages, umap_celltypes) {

  layout <- umap_out$layout
  labels <- umap_celltypes$InferredCellType
  
  
  main="UMAP for R1B1ROI1 Mixed Cells"
  palette = rainbow(length(unique(labels)))
  cex=0.6
  pch=19
  add=FALSE
  legend.suffix=""
  cex.main=1.5 # size of title
  cex.legend=0.85 
  
  
  xylim <- range(layout)
  xylim_big <- xylim + c(-6, 1)
  xylim_small <- xylim + c(-1, 1)
  
  
  par(mar=c(0.2,0.7,1.2,0.7), ps=10) # this is for formatting; mar is margin and ps is point size (font size)
  plot(xylim_big, xylim_small, type="n", axes=F, frame=F)
  rect(xylim_small[1], xylim_small[1], xylim_small[2], xylim_small[2], border="#aaaaaa", lwd=0.25)  
  
  unique_labels <- sort(unique(labels))
  label_colors <- palette[as.numeric(factor(labels, levels = unique_labels))]
  
  points(layout[,1], layout[,2], col=label_colors,
         cex=cex, pch=pch) # drawing the points
  mtext(side=3, main, cex=cex.main) # adding the main title

  legend.text <- as.character(unique_labels)
  legend("topleft", legend=legend.text,
         col=palette,
         bty="n", pch=pch, cex=cex.legend)
}

# plot whole umap
umap_plot(umap_out, umap_averages, umap_celltypes)
print("Computing umap with the following inferred celltypes:"); print(unique(umap_celltypes$InferredCellType))




umap_plot_onlyone <- function(umap_out, umap_averages, marker) {

  layout <- umap_out$layout
  
  #main=paste(marker, "expression levels on R1B1ROI1 mixed")
  main=paste(marker)
  palette = colorRampPalette(c('grey','blue'), bias=1)
  rev_palette =  colorRampPalette(c('blue','grey'), bias=1)
  pad=0.1 
  cex=0.6 # point size, can make smaller, originally 0.6
  pch=19
  add=FALSE
  cex.main=2
  
  quantiles <- quantile(umap_averages[[marker]], probs = seq(0, 1, by = 0.1), na.rm = TRUE, type = 6)
  breaks <- unique(round(unname(quantiles)))
  umap_averages$Col <- palette(length(breaks))[as.numeric(cut(umap_averages[[marker]], breaks = breaks))]
  #umap_averages$Col <- palette(10)[as.numeric(cut(umap_averages[[marker]],breaks = 10))]
  print(breaks)

  xylim <- range(layout)
  #xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  xylim_small <- xylim + c(-1, 1)
  xylim_big <- xylim + c(-5, 1)

  
  par(mar=c(0.2,0.7,1.2,0.7), ps=10) # this is for formatting; mar is margin and ps is point size (font size)
  plot(xylim_big, xylim_small, type="n", axes=F, frame=F)
  rect(xylim_small[1], xylim_small[1], xylim_small[2], xylim_small[2], border="#aaaaaa", lwd=0.25)  
  
  
  points(layout[,1], layout[,2], col=umap_averages$Col,
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)
  
  
  #percentiles <- seq(0, 1, by = 0.1)
  #quantiles <- quantile(umap_averages[[marker]], probs = seq(0, 1, by = 0.1), na.rm = TRUE, type = 6)
  #intervals <- seq(min(umap_averages[[marker]]), max(umap_averages[[marker]]), length.out = 11)
  
  legend_text <- round(breaks)
  #legend_text <- round(intervals)
  legend_image <- as.raster(matrix(rev_palette(length(breaks)), ncol=1))
  text(x=xylim_big[1]+1.5, y = seq(-7,9,l=length(breaks)), labels = legend_text)
  rasterImage(legend_image, xylim_big[1],-7, xylim_big[1]+0.5,9)

}

# plot umap specific for marker
umap_plot_onlyone(umap_out, umap_averages, "actin")

# "Madelenes selection": "actin", "cd3", "cd4",
#   "cd45", "cd45ro", "collageni", "cytokeratin",
#   "fibulin2", "lyve1", "podoplanin", "cd38", "cd138"





r1b1 <- rbind(r1d1roi1, r1d1roi2)

### ADJUSTED ONE


tsne_function <- function(df, pca=TRUE, perplexity=30, max_iter=1000, epoch=100, exclude_mixed=FALSE, exclude_other=FALSE, only_mixed=FALSE) {  
    
    if (exclude_mixed==TRUE) {
      df <- df[df$InferredCellType != "Mixed", ]
    }
    
    if (exclude_other==TRUE) {
      df <- df[df$InferredCellType != "Other", ]
    }
    
    if (only_mixed==TRUE) {
      df <- df[df$InferredCellType == "Mixed", ]
    }
  
    cellType = df$InferredCellType
    df = df[1:12]
    
    set.seed(20) # make sure its always the same
    tsne_out <- Rtsne(df, 
                      pca=pca,
                      perplexity=perplexity,
                      k=2,
                      max_iter=max_iter,
                      epoch=epoch)
    
    Y = as.data.frame(tsne_out$Y)
  
    return(list(Y, cellType))
  
}

# run tsne
outlist <- tsne_function(r1d1roi2, pca=TRUE, perplexity=30, max_iter=1000, epoch=100, exclude_mixed=TRUE, exclude_other=TRUE, only_mixed=FALSE)
Y <- outlist[[1]]
cellType <- outlist[[2]]

r1d1roi2_tsne <- Y
r1d1roi2_tsne$cellType <- cellType



##### UMAP #####

umap_function <- function(df, exclude_mixed=TRUE, exclude_other=TRUE, only_mixed=FALSE) {
  set.seed(20)
  if (exclude_mixed==TRUE) {
    df <- df[df$InferredCellType != "Mixed", ]
  }
  
  if (exclude_other==TRUE) {
    df <- df[df$InferredCellType != "Other", ]
  }
  
  if (only_mixed==TRUE) {
    df <- df[df$InferredCellType == "Mixed", ]
  }
  
  umap_celltypes <- df$InferredCellType
  print(nrow(df))
  print(length(umap_celltypes))
  
  df <- df[1:12]
  umap_out <- umap(df)

  return(list(umap_out, umap_celltypes))
}

df <- r1b1roi1
umap_list <- umap_function(df, exclude_mixed=FALSE, exclude_other=TRUE, only_mixed=FALSE)
layout <- umap_list[[1]]$layout
cellType <- umap_list[[2]]
Y <- as.data.frame(layout)
Y$cellType <- cellType

r1d1roi2_umap <- Y


# plot dimred
ggplot(Y,
       aes(x = V1, y = V2, col = cellType)) +
  
  geom_point() +
  
  labs(x = "UMAP - 1",
       y = "UMAP - 2") +
  
  ggtitle("UMAP for R1B1ROI1 with Mixed cells") +
  
  theme_classic()


# List of data frames with their names
data_frames_tsne <- list(R1B1ROI1 = r1b1roi1_tsne, 
                    R1B1ROI2 = r1b1roi2_tsne, 
                    R1C1ROI1 = r1c1roi1_tsne, 
                    R1C1ROI2 = r1c1roi2_tsne, 
                    R1D1ROI1 = r1d1roi1_tsne, 
                    R1D1ROI2 = r1d1roi2_tsne)

data_frames_umap <- list(R1B1ROI1 = r1b1roi1_umap, 
                         R1B1ROI2 = r1b1roi2_umap, 
                         R1C1ROI1 = r1c1roi1_umap, 
                         R1C1ROI2 = r1c1roi2_umap, 
                         R1D1ROI1 = r1d1roi1_umap, 
                         R1D1ROI2 = r1d1roi2_umap)

# Create a list to store plots without legends
plots <- lapply(names(data_frames_umap), function(name) {
  df <- data_frames_umap[[name]]
  ggplot(df, aes(x = V1, y = V2, col = cellType)) +
         geom_point(size=0.5) +
         labs(x = "UMAP - 1",
              y = "UMAP - 2",
              title = name) +
        theme_minimal() +
         theme_classic() +
    theme(legend.position = "none")  # Remove legends
})


# Create a plot to extract the legend
legend_plot <- ggplot(Y, aes(x = V1, y = V2, col = cellType)) +
                      geom_point() +
                      labs(x = "UMAP - 1",
                           y = "UMAP - 2") +
                      ggtitle("t-SNE for R1B1ROI2") +
                      theme_classic()
                      
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


















