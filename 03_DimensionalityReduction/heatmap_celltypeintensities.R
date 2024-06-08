# Code creating heatmap showing cell intensities and inferred cell types



##### !!! Try to do it so that the T-cells labels don't stick together #######



library(ComplexHeatmap)
set.seed(20)

r1d1 = rbind(r1d1roi1, r1d1roi2)
df = r1d1

df_norm = rbind(r1d1roi1_norm, r1d1roi2_norm)
df = df_norm
# Exclude Other
df = df[df$InferredCellType != "Other",]

# Exclude Mixed
#df = df[df$InferredCellType != "Mixed",]

celltype_order <- order(df$InferredCellType)
cellType = df$InferredCellType

# Get the corresponding cell type labels in the same order
celltype_labels <- df$InferredCellType[celltype_order]

# Plot the COMPLETE reordered heatmap with cell type labels
Heatmap(as.matrix(df[1:12 ]), name ="intensity", row_order = celltype_order, 
        row_split = cellType, row_title_rot = 0, column_title="Heatmap for R1D1 (normalised intensities)",  
        show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE)  # Use cell type labels as row side color


# Plot lymphocytes separately
lymphocytes = c("T-cells Memory", "T-cells CD4+ (helper)", "T-cells Mixed", "T-cells Naive", "B-cells Plasma")

# separate t cells:
df_t = df[df$InferredCellType %in% lymphocytes, ]
celltype_order_t <- order(df_t$InferredCellType)
cellType_t = df_t$InferredCellType

Heatmap(as.matrix(df_t[1:12]), name ="intensity", row_order = celltype_order_t, 
        row_split = cellType_t, row_title_rot = 0, column_title="Heatmap for R1B1 (raw intensities)",  
        show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE)  # Use cell type labels as row side order


df_nott = df[!(df$InferredCellType %in% lymphocytes), ]
celltype_order_nott <- order(df_nott$InferredCellType)
cellType_nott = df_nott$InferredCellType

Heatmap(as.matrix(df_nott[1:12]), name ="intensity", row_order = celltype_order_nott, 
        row_split = cellType_nott, row_title_rot = 0, column_title="Heatmap for R1B1 (raw intensities)",  
        show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE)  # Use cell type labels as row side order






###### OLD CODE ######

lymphocytes = c("T-cells Memory", "T-cells CD4+ (helper)", "T-cells Mixed", "T-cells Naive", "B-cells Plasma")

celltype_order <- order(celltypes$InferredCellType)
cellType = celltypes$InferredCellType

# Get the corresponding cell type labels in the same order
celltype_labels <- celltypes$InferredCellType[celltype_order]


# Plot the COMPLETE reordered heatmap with cell type labels
Heatmap(as.matrix(averages), name ="intensity", row_order = celltype_order, 
        row_split = cellType, row_title_rot = 0, column_title="Heatmap for R1B1ROI1 (raw intensities)",  
        show_column_dend = FALSE)  # Use cell type labels as row side color

# Normalisation or not?
averages_df = averages_scaled
#averages_df = averages


# separate t cells:
averages_t = averages_df[celltypes$InferredCellType %in% t_cells, ]
celltypes_t = celltypes[celltypes$InferredCellType %in% t_cells, ]
celltype_order_t <- order(celltypes_t$InferredCellType)
cellType_t = celltypes_t$InferredCellType

Heatmap(as.matrix(averages_t), name ="intensity", row_order = celltype_order_t, 
        row_split = cellType_t, row_title_rot = 0, column_title="Heatmap for R1B1ROI1 (raw intensities)",  
        show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE)  # Use cell type labels as row side order


averages_nott = averages_df[!(celltypes$InferredCellType %in% t_cells), ]
celltypes_nott = celltypes[!(celltypes$InferredCellType %in% t_cells), ]
celltype_order_nott <- order(celltypes_nott$InferredCellType)
cellType_nott = celltypes_nott$InferredCellType

Heatmap(as.matrix(averages_nott), name ="intensity", row_order = celltype_order_nott, 
        row_split = cellType_nott, row_title_rot = 0, column_title="Heatmap for R1B1ROI1 (raw intensities)",  
        show_column_dend = FALSE, show_row_dend = FALSE, show_row_names = FALSE)  # Use cell type labels as row side order




