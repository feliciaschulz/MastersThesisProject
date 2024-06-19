##############################################
# DIMENSIONALITY REDUCTION VISUALISATION APP 
# 
# This is a Shiny application with which a UMAP plot can be created
# and the result visualised.
# The data is loaded in and the original image is displayed as a 
# scatterplot with each cell as a point. 
# The umap can be created out of different subsets of the data.
# The umap plot appears on the interface and the user can hover over
# points to see them displayed in the original image.
# The user can also select a certain region in the umap plot, and the
# selected cells will be displayed on the original image. Double-click
# to remove the selected region.
#
# Usage: Run the app by clicking "Run App" in the upper right corner in
#   RStudio. Alternatively, run the command runApp() in the RStudio console.
#
# Required packages & versions:
#   - shiny: 1.8.1.1
#   - plotly: 4.10.4
#   - bslib: 0.7.0
#   - Rtsne: 0.17
#   - umap: 0.2.10.0
#   - dplyr: 1.1.4
#
# OBS: Make sure that the "intensities" data frame has the correct path.
#   The variable my_img can be changed to inspect different images, given the
#   format of the intensities data frame is the same.
#
# Known bugs / weaknesses: This application can have long loading times. Depending on the 
#   size of the subset selected for the umap, the plot creation can take up to
#   a few minutes. This is solely because the umap() function is computationally
#   complex. Also, hover slowly. When selecting a region, after drawing the 
#   circle around the chosen cells, leave the cursor in place and don't keep
#   hovering. Otherwise, the application will detect your movement as hovering
#   again and the selected region will disappear. The displaying of the selected
#   region can also take a little while, because identifying the cells and
#   referring to the original data frame with coordinates is computationally
#   complex as well.
# 
# Input: A data frame with cells in rows, biomarker intensity data in columns,
#   two columns for x and y coordinates (X_coord, Y_coord), area of cell (area),
#   cell phenotype (InferredCellType). 
#   Accessible in folder "../00_Data/IntensitiesCelltypes/"
###############################################



library(shiny)
library(plotly)
library(bslib)
library(Rtsne)
library(umap)
library(dplyr)

    
set.seed(20) # Make sure the umap result is always the same 

#### Load in data
my_img = "R1B1ROI1"
intensities = read.csv(paste0("../00_Data/IntensitiesCelltypes/", my_img, "_intensitiesdf.csv"))
intensities$X = NULL
geo = intensities 
geo$colour = "nohover"


#### Define the User Interface
ui <- fluidPage(
  titlePanel("Cell visualisation app"),
  sidebarLayout(
    sidebarPanel(
      # Input dropdown dimensionality reduction method
      selectInput(
        "method",
        "Select dimensionality reduction method",
        choices = list("UMAP"),
        selected = 1
      ),
      # Input selection cell type subset for umap
      checkboxGroupInput(
        inputId = "celltype",
        label = "Select all cell types that should be included",
        choices=list("Mixed", "Myocytes", "Other", "Lymphatic vessel", "Extracellular matrix", "T-cells Memory",
                     "T-cells CD4+ (helper)", "T-cells Mixed", "Epithelial cells", "T-cells Naive", "B-cells Plasma"),
        selected = 19
      ),
      # "Run" button
      actionButton(inputId="submit_button", label="Run"),
      width = 2
    ),
    mainPanel(
        plotlyOutput("dimredplot"), # umap plot
        textOutput("vector"),
        plotlyOutput("scatterplot"), # scatterplot
        textOutput("more_text")
        
        
    )
    
  )
)



#### Define the server
server <- function(input, output, session) {
  
  # Create input intensity data from subset of selection from selection button
  data_int <- eventReactive(input$submit_button, {
    celltype_choices = input$celltype # cell types chosen by user
    intensities_df <- intensities[intensities$InferredCellType %in% celltype_choices, ] # subset df with choices
    intensities_df <- intensities_df %>% select(-InferredCellType, -X_coord, -Y_coord) # remove non-intensity columns
    return(intensities_df)
  })
  
  # Create input data frame for celltypes 
  data_ct <- eventReactive(input$submit_button, {
    celltype_choices = input$celltype # cell types chosen by user
    celltypes_df <- intensities[intensities$InferredCellType %in% celltype_choices, ] # subset df with choices
    celltypes_df$colour = "nohover" # set as default, later, selected points will be individually changed to "hover" and then have a different colour
    celltypes_df$pointsize <- 2 # set as default, later, selected points will be individually changed to a larger size
    return(celltypes_df)
  })
  
  # create umap data frame
  compute <- eventReactive(input$submit_button, {
    method = input$method # method chosen by user
    set.seed(20) # just in case :)
    
    # run dimensionality reduction with umap
    umap_out <- umap(data_int())
    df <- data.frame(col1=umap_out$layout[,1], col2=umap_out$layout[,2]) # clean up output data frame
    df$cell_index <- rownames(df) # add cell_index so we can later refer back to the original row index and find coordinates
  
    return(df)
  })
  
  # Display text on umap plot while hovering over
  hover_text <- eventReactive(input$submit_button, {
    paste("Cell number:", 1:nrow(data_int()),"<br>", "Cell size:", data_int()$area,"<br>", "<br>", "Actin:         ", data_int()$actin, "<br>", 
                        "cd3:            ", data_int()$cd3, "<br>", "cd4:            ", data_int()$cd4, "<br>", "cd45:          ", data_int()$cd45, "<br>", 
                        "cd45ro:      ", data_int()$cd45ro, "<br>", "Collagen1:  ", data_int()$collageni, "<br>", "Cytokeratin:", data_int()$cytokeratin, "<br>", 
                        "Fibulin2:     ", data_int()$fibulin2, "<br>", "Podoplanin:", data_int()$podoplanin, "<br>",
                        "Lyve1:        ", data_int()$lyve1, "<br>", data_int()$cd38, "<br>", data_int()$cd138, "<br>")
  })
  

  # Display information about how many cells and how to use
  output$vector <- renderText({
    paste("There are", nrow(compute()), "cells in your selection.\nYou have used the UMAP method.\nClick and drag events (i.e., select/lasso) appear below (double-click to clear).\n
          Hover events also appear below.")
  })

  
  # Dimensionality reduction plot
  output$dimredplot <- renderPlotly({
    p <- plot_ly(data = compute(), x = ~col1, y = ~col2, color = ~data_ct()$InferredCellType, colors = "Spectral", type = "scatter", 
            mode = "markers", text=hover_text(), source = "plot1") %>%
      layout(
        xaxis = list(title = ("UMAP - 1")),
        yaxis = list(title = ("UMAP - 2")),
        title = paste("UMAP for ", my_img, " with nuclei and cyto segmentation"),
        showlegend = TRUE, 
        dragmode = "lasso"
      )
    #fig['layout']['yaxis']['autorange'] = "reversed"
    event_register(p, "plotly_hover") # Register user events
    event_register(p, "plotly_selected")
    return(p)
  })
  
  # Scatterplot with coordinates
  output$scatterplot <- renderPlotly({
    plot_ly(geo, x = ~X_coord, y = ~Y_coord, type = "scatter", color = ~colour, colors = c("grey"), 
            mode = "markers", marker = list(color = "grey", size = 2)) %>%
      layout(title = "Original image", yaxis = list(autorange = "reversed"))
  })

  
  # Code for registering lasso selection points and displaying on scatterplot
  observeEvent(event_data("plotly_selected", source = "plot1"), {
    selected_points <- event_data("plotly_selected", source = "plot1") # the data registered by the app
    if (!is.null(selected_points)) { # if a selection has been made
      indices <- selected_points$pointNumber + 1  # Adjust for 1-based indexing
      dimred_x <- selected_points$x # data frame with umap x coords
      dimred_y <- selected_points$y # data frame with umap y coords

      
      # Load in umap data frame again, but just coordinates
      # The coordinates registered by the event_data can't be directly accessed in the 
      # origin dimensionality reduction df, because the values are slightly changed
      # by the application, so we have to find them by looping
      # through and finding close matches
      dim_red_df <- compute() %>% select(-cell_index)
      data_rounded <- round(dim_red_df, 2) # getting the data from here
      tolerance <- 1e-6 # Defining a tolerance for imperfect matches
      
      # Initialise vector for collecting cell indexes
      cell_indexes_selected <- c()
      
      # Loop through x and y coordinates
      for (i in seq_along(dimred_x)) {
        row_to_check <- c(round(dimred_x[i], 2), round(dimred_y[i], 2))
        row_index <- which(apply(data_rounded, 1, function(row) all(abs(row - row_to_check) < tolerance))) # Check if any row matches within tolerance
        
        # Might find multiple close matches, choose first one
        if (length(row_index) == 1) {
          row_index <- row_index
        } else if (length(row_index ) > 1) {
          row_index <- row_index[1]
        } else {
          row_index <- 0
        }
        
        # If a matching row has been found
        if (row_index != 0) {
          cell_index <- compute()[row_index, ]$cell_index # access original cell index from row in dim red df
          cell_indexes_selected <- c(cell_indexes_selected, cell_index) # add cell index to vector
        }
      }
      
      # If the cell index vector has cells found, display in scatterplot
      if (!is.null(cell_indexes_selected)) {
        geo$colour <- "lightgrey" # set all points back to grey
        geo[cell_indexes_selected, ]$colour <- "black" # change only selected cell indexes to black
        
        geo$pointsize <- 2 # set all points back to small
        geo[cell_indexes_selected, ]$pointsize <- 7 # change only selected cell indexes to big
        
        # Change original df with new columns and update plot
        plotlyProxy("scatterplot", session) %>%
          plotlyProxyInvoke("restyle", list(marker = list(color = geo$colour, size = geo$pointsize)))
      } 
      
    

    }
  })

  # Code for registering lasso selection points and displaying on scatterplot
  observeEvent(event_data("plotly_hover", source = "plot1"), {
    hover_info <- event_data("plotly_hover", source = "plot1") # the data registered by the app
    x <- round(hover_info$x, 2) # x coord
    y <- round(hover_info$y, 2) # y coord
    row_to_check <- c(x, y) # total coord as vector

    
    # Load in umap data frame again, but just coordinates
    # The coordinates registered by the event_data can't be directly accessed in the 
    # origin dimensionality reduction df, because the values are slightly changed
    # by the application, so we have to find them by looping
    # through and finding close matches
    dim_red_df <- compute() %>% select(-cell_index)
    data_rounded <- round(dim_red_df, 2) # getting the data from here
    tolerance <- 1e-6 # Defining a tolerance for imperfect matches
    row_index <- which(apply(data_rounded, 1, function(row) all(abs(row - row_to_check) < tolerance))) # Check if any row matches within tolerance
    # getting the row index of my geographical coordinate data frame to find the hover cell coord ^^^^^

    if (length(row_index) == 1) {
      row_index = row_index
    } else if (length(row_index ) > 1) {
      row_index = row_index[1]
    }  else {
      row_index = 0
    }
  
    # if it has found a cell, display in scatterplot
    if (!is.null(row_index)) {
      cell_index = compute()[row_index, ]$cell_index # access original cell index from row in dim red df
      
      geo$colour <- "lightgrey" # set all points back to grey
      geo[cell_index, ]$colour <- "black" # change selected point to black

      geo$pointsize <- 2 # set all points back to small
      geo[cell_index, ]$pointsize <- 7 # change selected point to big

      # Change original df with new columns and update plot
      plotlyProxy("scatterplot", session) %>%
        plotlyProxyInvoke("restyle", list(marker  = list(color = geo$colour, size = geo$pointsize)))
    }

})

  
}


shinyApp(ui = ui, server = server)



