# -------------------------------
# Create explicit synthetic dataset
# -------------------------------

cells <- data.frame(
  cell_id = paste0("Cell", 1:12),
  CD3D  = c(10.5, 9.2, 11.3, 8.7, 1.1, 0.8, 1.5, 0.6, 0.9, 1.2, 0.7, 1.0),
  MS4A1 = c(0.5, 0.2, 1.0, 0.8, 10.1, 9.5, 11.0, 8.8, 0.4, 0.9, 0.3, 0.6),
  LYZ   = c(1.0, 0.5, 0.8, 0.6, 0.9, 1.1, 0.7, 0.5, 9.8, 11.2, 10.5, 8.9),
  MKI67 = c(6.2, 1.8, 0.9, 3.4, 7.7, 9.4, 2.5, 5.1, 1.3, 6.0, 2.2, 9.8),
  cluster = rep(c("Tcell", "Bcell", "Myeloid"), each = 4),
  stringsAsFactors = FALSE
)

print(cells)

# Inspect the Dataset
str(cells)
summary(cells[,2:5])   # numeric only
table(cells$cluster)    # categorical only

# Compute per-cell total expression (sum across the 4 genes)
cells$total <- apply(cells[, c("CD3D","MS4A1","LYZ","MKI67")], 1, sum)
cells$total

# Label cells "Active" if total > median(total), otherwise "Resting"
median_total <- median(cells$total)
cells$state <- ifelse(cells$total > median_total, "Active", "Resting")
table(cells$state)
print(median_total)
print(cells)

# Using tapply() (or aggregate) to compute cluster-wise mean for each gene.
aggregate(cells[, c("CD3D","MS4A1","LYZ","MKI67")], 
          by=list(cluster=cells$cluster), FUN=mean)

# To check which single cell has highest MKI67? Return cell_id and value.
idx <- which.max(cells$MKI67)
cells$cell_id[idx]; cells$MKI67[idx]

# For each cell compute the marker with the highest value (max per row)
cells$max_gene <- apply(cells[, c("CD3D","MS4A1","LYZ","MKI67")], 1, function(x) names(x)[which.max(x)])
cells$max_gene

# Tabulate the number of cells with each max_geneulate 
# the number of cells with each max_gene per cluster     
table(cells$cluster, cells$max_gene)

# Count how many cells per cluster and Active vs Resting
table(cells$cluster, cells$state)
# Proportion of Active vs Resting per cluster
prop.table(table(cells$cluster, cells$state), margin=1)

# Create expr_mat (12×4) from gene columns only; 
# set rownames = cell_id; confirm identical rownames
expr_mat <- as.matrix(cells[, c("CD3D","MS4A1","LYZ","MKI67")])
rownames(expr_mat) <- cells$cell_id
all(rownames(expr_mat) == cells$cell_id)

# Compute column means on expr_mat; 
# check which marker has the highest overall mean
colMeans(expr_mat)

# Compute row ranges (max - min). Check which cell is most uneven across markers
row_ranges <- apply(expr_mat, 1, function(x) max(x) - min(x))
which.max(row_ranges); rownames(expr_mat)[which.max(row_ranges)]

# Subset to only Tcell and recompute mean of CD3D and compare to overall mean
mean_Tcell_CD3D <- mean(cells$CD3D[cells$cluster == "Tcell"])
mean_overall_CD3D <- mean(cells$CD3D)
mean_Tcell_CD3D
mean_overall_CD3D

# Rename column MS4A1 → MS4A1_B and verify
colnames(expr_mat)[colnames(expr_mat) == "MS4A1"] <- "MS4A1_B"
"MS4A1_B" %in% colnames(expr_mat)

# Reorder cells by descending MKI67 without using order(-x) 
# (get descending index manually)
idx_asc <- order(cells$MKI67)
idx_desc <- rev(idx_asc)
cells_desc <- cells[idx_desc, ]
cells_desc$cell_id

# Create a list scRNA containing cells, expr_mat, and markers 
# c("CD3D","MS4A1_B","LYZ","MKI67")
scRNA <- list(cells = cells, expr_mat = expr_mat, markers = c("CD3D","MS4A1_B","LYZ","MKI67"))
names(scRNA)

# From scRNA, extract the LYZ values for all Myeloid cells (nested subsetting)
scRNA$cells$LYZ[scRNA$cells$cluster == "Myeloid"]

# Create a simple barplot of cluster-wise means for MKI67
mk_means <- tapply(cells$MKI67, cells$cluster, mean)
barplot(mk_means, main="Cluster-wise MKI67 means", 
        ylab="Mean MKI67", xlab="Cluster")

# Create a boxplot of MKI67 values per cluster
boxplot(MKI67 ~ cluster, data=cells, 
        main="Cluster-wise MKI67 distribution", 
        ylab="MKI67", xlab="Cluster")

# Use sapply() on gene columns to return standard deviations per gene
sapply(cells[, c("CD3D","MS4A1","LYZ","MKI67")], sd)

# Check which cluster is most active on average
tapply(cells$total, cells$cluster, mean)

# Check which cluster has the most MKI67-high cells (MKI67 > 5)
table(cells$cluster[cells$MKI67 > 5])

# Zero all expression values < 1 (threshold) on a copy of expr_mat
expr_copy <- expr_mat
changed_idx <- expr_copy < 1
num_changed <- sum(changed_idx)
expr_copy[changed_idx] <- 0
num_changed
expr_copy
# Number of entries changed (set to 0): 18

# Interpretation:
# The marker expression patterns clearly validate the assigned cluster labels. 
# T cells are distinguished by strong CD3D expression, 
# B cells by MS4A1, and Myeloid cells by LYZ, 
# while MKI67 provides information on proliferation across clusters. 
# The per-cell and cluster-level summaries, maximum marker identification, 
# and distribution plots consistently show that each cluster predominantly 
# expresses its characteristic marker with minimal overlap. 
# Together, these results confirm that the cluster labels align well 
# with the biological identity of the cells.
