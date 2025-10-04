Single-Cell Marker Analysis
===========================

**Author:** Merylin Ogunlola

**Date:** 2025-10-04

Introduction
------------

This analysis uses a synthetic single-cell dataset with expression levels of four marker genes: **CD3D, MS4A1, LYZ, MKI67**. Twelve cells are divided into three clusters (**Tcell, Bcell, Myeloid**). The goal is to verify whether marker expression patterns support the given cluster labels and to examine proliferation and activity states.

Dataset
-------

`cells <- data.frame(cell_id = paste0("Cell", 1:12),
CD3D  = c(10.5, 9.2, 11.3, 8.7, 1.1, 0.8, 1.5, 0.6, 0.9, 1.2, 0.7, 1.0),    MS4A1 = c(0.5, 0.2, 1.0, 0.8, 10.1, 9.5, 11.0, 8.8, 0.4, 0.9, 0.3, 0.6),    LYZ   = c(1.0, 0.5, 0.8, 0.6, 0.9, 1.1, 0.7, 0.5, 9.8, 11.2, 10.5, 8.9),    MKI67 = c(6.2, 1.8, 0.9, 3.4, 7.7, 9.4, 2.5, 5.1, 1.3, 6.0, 2.2, 9.8),    cluster = rep(c("Tcell", "Bcell", "Myeloid"), each = 4),    stringsAsFactors = FALSE  )   `

*   **12 cells**
    
*   **4 marker genes**
    
*   **3 clusters**, balanced (4 cells each)
    

Per-Cell Expression and Activity
--------------------------------

`cells$total <- apply(cells[, c("CD3D","MS4A1","LYZ","MKI67")], 1, sum)  median_total <- median(cells$total)  cells$state <- ifelse(cells$total > median_total, "Active", "Resting")   `

*   Cells with total > median = **Active**
    
*   Cells with total ≤ median = **Resting**
    

About half the cells fall into each group.

Cluster-Level Gene Means
------------------------

`aggregate(cells[, c("CD3D","MS4A1","LYZ","MKI67")],             by=list(cluster=cells$cluster), FUN=mean)   `

*   **T cells** → CD3D-high
    
*   **B cells** → MS4A1-high
    
*   **Myeloid** → LYZ-high
    

👉 Marker patterns **match known immune cell biology**.

Proliferation Marker MKI67
--------------------------

`which.max(cells$MKI67)   # highest MKI67 cell   `

*   Highest MKI67 = **Cell12 (9.8)**
    
*   Some B cells and Myeloid cells also show strong proliferation
    

Plots:

`barplot(tapply(cells$MKI67, cells$cluster, mean),           main="Cluster-wise MKI67 means")`
<img width="613" height="351" alt="37ce1c3d-af28-4abc-9369-457ed3a757c2" src="https://github.com/user-attachments/assets/306d7859-1c97-4e22-99a9-9b2c05b6e37c" />




`boxplot(MKI67 ~ cluster, data=cells,           main="Cluster-wise MKI67 distribution")`
<img width="613" height="351" alt="00fd6272-3f08-4895-91b6-62ce0e9317eb" src="https://github.com/user-attachments/assets/da04d5de-6f5a-4dca-b1de-2ac834275d41" />




Marker Dominance per Cell
-------------------------

`cells$max_gene <- apply(cells[, c("CD3D","MS4A1","LYZ","MKI67")], 1,                           function(x) names(x)[which.max(x)])  table(cells$cluster, cells$max_gene)   `

*   Most cells are dominated by their **canonical marker**
    
*   A few cells dominated by **MKI67**, reflecting proliferative dominance
    

Matrix-Based Analysis
---------------------

`expr_mat <- as.matrix(cells[, c("CD3D","MS4A1","LYZ","MKI67")])  rownames(expr_mat) <- cells$cell_id  colMeans(expr_mat)  # gene-wise means  sapply(expr_mat, sd) # standard deviations   `

*   **LYZ** and **MKI67** have the highest overall means
    
*   **MS4A1** and **LYZ** show highest variability between clusters
    

Row range analysis:

`row_ranges <- apply(expr_mat, 1, function(x) max(x) - min(x))  rownames(expr_mat)[which.max(row_ranges)]   `

*   Most uneven expression = a **Myeloid cell** with strong LYZ compared to other markers
    

Thresholding Low Expression
---------------------------

`expr_copy <- expr_mat  changed_idx <- expr_copy < 1  num_changed <- sum(changed_idx)  expr_copy[changed_idx] <- 0  num_changed   `

*   **18 low-expression values (<1)** were set to 0 (noise filtering)
    

Final Interpretation
--------------------

The results **strongly support cluster labels**:

*   **T cells** are CD3D-high
    
*   **B cells** are MS4A1-high
    
*   **Myeloid cells** are LYZ-high
    
*   **MKI67** captures proliferation, with some cells dominated by this marker
    

This alignment confirms that marker expression is consistent with immune cell identity. The presence of MKI67-dominant cells highlights biological heterogeneity: some cells are strongly proliferative, regardless of lineage.
