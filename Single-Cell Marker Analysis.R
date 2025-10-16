# Create  synthetic dataset
cells <- data.frame(
  cell_id = paste0("Cell", 1:12),
  CD3D  = c(10.5, 9.2, 11.3, 8.7, 1.1, 0.8, 1.5, 0.6, 0.9, 1.2, 0.7, 1.0),
  MS4A1 = c(0.5, 0.2, 1.0, 0.8, 10.1, 9.5, 11.0, 8.8, 0.4, 0.9, 0.3, 0.6),
  LYZ   = c(1.0, 0.5, 0.8, 0.6, 0.9, 1.1, 0.7, 0.5, 9.8, 11.2, 10.5, 8.9),
  MKI67 = c(6.2, 1.8, 0.9, 3.4, 7.7, 9.4, 2.5, 5.1, 1.3, 6.0, 2.2, 9.8),
  cluster = rep(c("Tcell", "Bcell", "Myeloid"), each = 4),
  stringsAsFactors = FALSE
)

cells
# This step mimics real single-cell RNA data in a simplified form.

# Inspect the dataset structure and summary statistics
str(cells)
summary(cells)
# Check for missing values
colSums(is.na(cells))

# Compute Per-Cell Total Expression
cells$total_expression <- apply(cells[, c("CD3D", "MS4A1", "LYZ", "MKI67")], 1, sum)
cells$total_expression
# Display the updated data frame with total expression
cells

# Label cells as "Active" or "Resting"
median_expression <- median(cells$total_expression)
cells$activity_status <- ifelse(cells$total_expression > median_expression, "Active", "Resting")
cells$activity_status
# Display the updated data frame with activity status
cells

# Compute Cluster-Wise Mean Expression Using tapply()
# Mean CD3D expression by cluster
mean_CD3D_by_cluster <- tapply(cells$CD3D, cells$cluster, mean)
mean_CD3D_by_cluster
# Mean MS4A1 expression by cluster
mean_MS4A1_by_cluster <- tapply(cells$MS4A1, cells$cluster, mean)
mean_MS4A1_by_cluster
# Mean LYZ expression by cluster
mean_LYZ_by_cluster <- tapply(cells$LYZ, cells$cluster, mean)
mean_LYZ_by_cluster
# Mean MKI67 expression by cluster
mean_MKI67_by_cluster <- tapply(cells$MKI67, cells$cluster, mean)
mean_MKI67_by_cluster
# Combine all mean expressions into a data frame
mean_expression_by_cluster <- data.frame(
  cluster = names(mean_CD3D_by_cluster),
  mean_CD3D = as.numeric(mean_CD3D_by_cluster),
  mean_MS4A1 = as.numeric(mean_MS4A1_by_cluster),
  mean_LYZ = as.numeric(mean_LYZ_by_cluster),
  mean_MKI67 = as.numeric(mean_MKI67_by_cluster)
)
mean_expression_by_cluster

# Identify the cells with the highest MKI67 expression
max_MKI67 <- which.max(cells$MKI67) # To find the row index of the maximum value
highest_MKI67_cell <- cells[max_MKI67, c("cell_id", "MKI67")]
highest_MKI67_cell

# For each cell, compute max gene (the marker with the highest value). 
# Produce a vector of marker names per cell.
marker_names <- apply(cells[, c("CD3D", "MS4A1", "LYZ", "MKI67")], 1, function(x) {
  names(x)[which.max(x)]
})
marker_names
# Add this information to the data frame
cells$max_marker <- marker_names
cells

# Count how many cells per cluster and how many Active vs Resting (table).
cluster_counts <- table(cells$cluster)
cluster_counts
activity_counts <- table(cells$activity_status)
activity_counts
# Combine cluster and activity status counts
combined_counts <- table(cells$cluster, cells$activity_status)
combined_counts
# Convert to data frame for better readability
combined_counts_df <- as.data.frame(combined_counts)
combined_counts_df

# Create a matrix expr_mat (12×4) from the four gene columns
expr_mat <- as.matrix(cells[, c("CD3D", "MS4A1", "LYZ", "MKI67")])
rownames(expr_mat) <- cells$cell_id
colnames(expr_mat) <- c("CD3D", "MS4A1", "LYZ", "MKI67")
expr_mat
# Confirm that the rownames match the cell_id column
all(rownames(expr_mat) == cells$cell_id)
# Confirm that the colnames match the gene names
all(colnames(expr_mat) == c("CD3D", "MS4A1", "LYZ", "MKI67"))

# Compute the mean expression of each gene across all cells using colMeans()
mean_expression_genes <- colMeans(expr_mat)
mean_expression_genes
# Identify the marker(gene) with the highest mean expression
highest_mean_expression_gene <- names(which.max(mean_expression_genes))
highest_mean_expression_gene

# Compute row ranges (max–min) , which cell is most uneven across markers?
row_ranges <- apply(expr_mat, 1, function(x) max(x) - min(x))
row_ranges
# Identify the cell with the highest range (most uneven expression)
most_uneven_cell <- names(which.max(row_ranges))
most_uneven_cell

# Subset the matrix to include only T cells (cluster == "Tcell")
tcell_indices <- which(cells$cluster == "Tcell")
tcell_expr_mat <- expr_mat[tcell_indices, ]
tcell_expr_mat
# Recompute the mean expression of CD3D across T cells
mean_CD3D_tcells <- mean(tcell_expr_mat[, "CD3D"])
mean_CD3D_tcells
# Compare the overall mean CD3D expression across all cells
overall_mean_CD3D <- mean(expr_mat[, "CD3D"])
overall_mean_CD3D

# Rename column MS4A1 to MS4A1_B and verify
colnames(expr_mat)[colnames(expr_mat) == "MS4A1"] <- "MS4A1_B"
colnames(expr_mat)
# Verify the change
all(colnames(expr_mat) == c("CD3D", "MS4A1_B", "LYZ", "MKI67"))
# Display the updated column names
colnames(expr_mat)

# Get the ascending order index for MKI67
ascending_index <- order(expr_mat[, "MKI67"])
# Reverse the index to get descending order
descending_index <- rev(ascending_index)
# Reorder the data frame based on descending MKI67
cells_desc_MKI67 <- cells[descending_index, ]
cells_desc_MKI67
# Verify that the MKI67 values are in descending order
all(diff(cells_desc_MKI67$MKI67) <= 0)

# Create a list called scRNA containing cells,expr_mat and markers(CD3D, MS4A1_B, LYZ, MKI67)).
scRNA <- list(
  cells = cells,
  expr_mat = expr_mat,
  markers = c("CD3D", "MS4A1_B", "LYZ", "MKI67")
)
scRNA
names(scRNA)
# Display the structure of the list
str(scRNA)

# From scRNA, extract the LYZ values for all Myeloid cells
myeloid_lyz <- scRNA$cells$LYZ[scRNA$cells$cluster == "Myeloid"]
myeloid_lyz

# Create a simple barplot of cluster-wise means for MKI67
# Compute mean MKI67 expression per cluster
mean_MKI67_by_cluster <- tapply(scRNA$cells$MKI67, scRNA$cells$cluster, mean)
mean_MKI67_by_cluster
# Create the barplot
barplot(mean_MKI67_by_cluster, 
        main = "Mean MKI67 Expression by Cluster", 
        xlab = "Cluster", 
        ylab = "Mean MKI67 Expression", 
        col = c("lightblue", "lightgreen", "lightcoral"))

# Use sapply() on gene columns to return standard deviations per gene
# Use sapply() to compute standard deviation for each gene column
gene_sd <- sapply(scRNA$cells[, c("CD3D", "MS4A1", "LYZ", "MKI67")], sd)
gene_sd

# Which cluster is most active on average?
mean_total_expression_by_cluster <- tapply(scRNA$cells$total_expression, scRNA$cells$cluster, mean)
mean_total_expression_by_cluster
most_active_cluster <- names(which.max(mean_total_expression_by_cluster))
most_active_cluster

# Make a copy of expr_mat
expr_mat_copy <- scRNA$expr_mat
# Count how many entries are below 1 before thresholding
num_below_threshold <- sum(expr_mat_copy < 1)
num_below_threshold
# Apply thresholding: set values below 1 to 0
expr_mat_copy[expr_mat_copy < 1] <- 0
# Count how many entries are below 1 after thresholding
num_below_threshold_after <- sum(expr_mat_copy < 1)
num_below_threshold_after

# The marker expression patterns strongly support the assigned cluster labels.
# The T cell cluster shows high expression of CD3D, a known T-cell marker;
# the B cell cluster shows high expression of MS4A1, a canonical B-cell marker;
# and the Myeloid cluster shows elevated LYZ, a myeloid cell marker.

# Additionally, the proliferation marker MKI67 varies across clusters, indicating 
# differences in cell activity but not mislabeling.
# Overall, the observed gene expression trends are consistent with known 
# biological signatures of these immune cell types, confirming that the cluster 
# labels are biologically valid and well supported by the marker data.


