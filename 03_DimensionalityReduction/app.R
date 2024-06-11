library(shiny)
library(plotly)
library(bslib)
library(Rtsne)
library(umap)
library(dplyr)

    
set.seed(20) # make sure its always the same
my_img = "R1B1ROI1"
intensities = read.csv(paste0("../00_Data/IntensitiesCelltypes/", my_img, "_intensitiesdf.csv"))
geo2 = intensities #%>% select(X_coord, Y_coord)
geo2$colour = "nohover"


# Define UI ----
ui <- fluidPage(
  titlePanel("Cell visualisation app"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "method",
        "Select dimensionality reduction method",
        choices = list("t-SNE", "UMAP"),
        selected = 1
      ),
      checkboxGroupInput(
        inputId = "celltype",
        label = "Select all cell types that should be included",
        choices=list("Mixed", "Myocytes", "Other", "Lymphatic vessel", "Extracellular matrix", "T-cells Memory",
                     "T-cells CD4+ (helper)", "T-cells Mixed", "Epithelial cells", "T-cells Naive", "B-cells Plasma"),
        selected = 1
      ),
      actionButton(inputId="submit_button", label="Run"),
      width = 2
    ),
    mainPanel(
      
        plotlyOutput("dimredplot"), 
        textOutput("vector"),
        plotlyOutput("scatterplot"), 
        textOutput("more_text"),
        #imageOutput("coloured_img")
        verbatimTextOutput("brush_info")
        
      
    # ),
    # fluidRow(
    #   column(width = 6,
    #          h4("Brushed points"),
    #          verbatimTextOutput("brush_info")
      #)
    )
    
  )
)



# Define server logic ----
server <- function(input, output, session) {
  data_int <- eventReactive(input$submit_button, {
    celltype_choices = input$celltype
    intensities_df <- intensities[intensities$InferredCellType %in% celltype_choices, ]
    intensities_df <- intensities_df %>% select(-InferredCellType, -X_coord, -Y_coord)
    return(intensities_df)
  })
  
  data_ct <- eventReactive(input$submit_button, {
    celltype_choices = input$celltype
    celltypes_df <- intensities[intensities$InferredCellType %in% celltype_choices, ]
    celltypes_df$colour = "nohover"
    celltypes_df$pointsize <- 2
    return(celltypes_df)
  })
  
  
  compute <- eventReactive(input$submit_button, {
    method = input$method
    set.seed(20)
    if (method=="t-SNE") {
      tsne_out <- Rtsne(data_int(),
                    pca=TRUE,
                    perplexity=30,
                    k=2,
                    max_iter=500,
                    epoch=100)
      Y <- as.data.frame(tsne_out$Y)
      df <- data.frame(col1=Y$V1, col2=Y$V2)
      df$cell_index <- rownames(df)

    } else { #method==UMAP
      umap_out <- umap(data_int())
      df <- data.frame(col1=umap_out$layout[,1], col2=umap_out$layout[,2])
      df$cell_index <- rownames(df)
    }
    return(df)
  })
  
  hover_text <- eventReactive(input$submit_button, {
    paste("Cell number:", 1:nrow(data_int()),"<br>", "Cell size:", data_int()$area,"<br>", "<br>", "Actin:         ", data_int()$actin, "<br>", 
                        "cd3:            ", data_int()$cd3, "<br>", "cd4:            ", data_int()$cd4, "<br>", "cd45:          ", data_int()$cd45, "<br>", 
                        "cd45ro:      ", data_int()$cd45ro, "<br>", "Collagen1:  ", data_int()$collageni, "<br>", "Cytokeratin:", data_int()$cytokeratin, "<br>", 
                        "Fibulin2:     ", data_int()$fibulin2, "<br>", "Podoplanin:", data_int()$podoplanin, "<br>",
                        "Lyve1:        ", data_int()$lyve1, "<br>", data_int()$cd38, "<br>", data_int()$cd138, "<br>")
  })
  

  
  # output$vector <- renderText({
  #   paste("There are", nrow(compute()), "cells in your selection.\n You have used the", 
  #         input$method, "method.")
  # })
  # 
  

  output$dimredplot <- renderPlotly({
    p <- plot_ly(data = compute(), x = ~col1, y = ~col2, color = ~data_ct()$InferredCellType, colors = "Spectral", type = "scatter", 
            mode = "markers", text=hover_text(), source = "plot1") %>%
      layout(
        xaxis = list(title = paste(input$method, "- 1")),
        yaxis = list(title = paste(input$method, "- 2")),
        title = paste(input$method, " for ", my_img, " with nuclei and cyto segmentation"),
        showlegend = TRUE, 
        dragmode = "lasso"
      )
    event_register(p, "plotly_hover")
    event_register(p, "plotly_selected")
    return(p)
  })
  
  output$scatterplot <- renderPlotly({
    plot_ly(geo2, x = ~X_coord, y = ~Y_coord, type = "scatter", color = ~colour, colors = c("grey"), 
            mode = "markers", marker = list(color = "grey", size = 2)) %>%
      layout(title = "Original image")
  })
  
  # Output brushed points
  output$brush_info <- renderPrint({
    d <- event_data("plotly_selected")
    if (is.null(d)) {
      "Click and drag events (i.e., select/lasso) appear here (double-click to clear)"
    } else {
      d
    }
  })
  
  observeEvent(event_data("plotly_selected", source = "plot1"), {
    selected_points <- event_data("plotly_selected", source = "plot1")
    if (!is.null(selected_points)) {
      indices <- selected_points$pointNumber + 1  # Adjust for 1-based indexing
      dimred_x <- selected_points$x
      dimred_y <- selected_points$y

      
      dim_red_df <- compute() %>% select(-cell_index)
      data_rounded <- round(dim_red_df, 2) # getting the data from here
      tolerance <- 1e-6 # Defining a tolerance for imperfect matches
      
      cell_indexes_selected <- c()
      
      for (i in seq_along(dimred_x)) {
        row_to_check <- c(round(dimred_x[i], 2), round(dimred_y[i], 2))
        row_index <- which(apply(data_rounded, 1, function(row) all(abs(row - row_to_check) < tolerance))) # Check if any row matches within tolerance
        
        if (length(row_index) == 1) {
          row_index <- row_index
        } else if (length(row_index ) > 1) {
          row_index <- row_index[1]
        } else {
          row_index <- 0
        }
        
        if (row_index != 0) {
          cell_index <- compute()[row_index, ]$cell_index
          
          cell_indexes_selected <- c(cell_indexes_selected, cell_index)

        }
      }
      
      if (!is.null(cell_indexes_selected)) {
        geo2$colour <- "lightgrey"
        geo2[cell_indexes_selected, ]$colour <- "black"
        
        geo2$pointsize <- 2
        geo2[cell_indexes_selected, ]$pointsize <- 7
        
        plotlyProxy("scatterplot", session) %>%
          plotlyProxyInvoke("restyle", list(marker = list(color = geo2$colour, size = geo2$pointsize)))
      }
      
    

    }
  })


  observeEvent(event_data("plotly_hover", source = "plot1"), {
    hover_info <- event_data("plotly_hover", source = "plot1")
    x <- round(hover_info$x, 2)
    y <- round(hover_info$y, 2)
    row_to_check <- c(x, y)

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
  


    if (!is.null(row_index)) {
      cell_index = compute()[row_index, ]$cell_index
      
      x_coord <- geo2[cell_index, "X_coord"] # geo has my coordinates
      y_coord <- geo2[cell_index, "Y_coord"]


      geo2$colour <- "lightgrey"
      geo2[cell_index, ]$colour <- "black"

      geo2$pointsize <- 2
      geo2[cell_index, ]$pointsize <- 7

      plotlyProxy("scatterplot", session) %>%
        plotlyProxyInvoke("restyle", list(marker  = list(color = geo2$colour, size = geo2$pointsize)))
    }

})



  
}


shinyApp(ui = ui, server = server)


# # Render the image and draw a point at the extracted coordinates
# output$geoPlot <- renderPlot({
#   # Render the image
#   img <- readJPEG("C-000_S-000_S_DAPI_R-01_W-B-1_ROI-01_A-DAPI.jpg")
#   rasterImage(img, 0, 0, 1, 1)
#   
#   # Draw the point on the image
#   points(x_coord, y_coord, col = "red", pch = 16)
# }, res = 96, width = 800, height = 600)


