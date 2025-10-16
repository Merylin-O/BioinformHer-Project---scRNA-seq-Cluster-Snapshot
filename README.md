
# Single-Cell Marker Analysis (Base R Project)

### **Author:** Merylin Ogunlola

### **Date:** October 2025

---

## **Project Overview**

This project analyses a **synthetic single-cell dataset** of 12 cells, each quantified for 4 marker genes — **CD3D**, **MS4A1**, **LYZ**, and **MKI67**.
Cells were manually grouped into three clusters (**Tcell**, **Bcell**, and **Myeloid**) to verify that marker expression patterns align with expected cell identities.

All analysis was performed using base R only.

---

## **Objectives**

* Create and inspect a small synthetic single-cell dataset.
* Calculate total expression and classify cells as *Active* or *Resting*.
* Compute cluster-wise mean expression for each gene using `tapply()`.
* Identify cell-specific and cluster-specific trends.
* Visualize proliferation (MKI67) across clusters using **base R plots**.
* Assess whether marker expression patterns support cluster labels.

---

## **Methods Summary**

1. **Data Creation** – Generated a 12×4 dataset using `data.frame()`.
2. **Exploration** – Used `str()` and `summary()` for structure and distribution.
3. **Computation** – Applied `apply()`, `tapply()`, `colMeans()`, and `sapply()` for summaries.
4. **Classification** – Labelled cells “Active” or “Resting” based on total expression.
5. **Visualization** – Created a **base R barplot** showing mean MKI67 by cluster.
6. **Verification** – Compared marker trends with expected cluster identities.

---

## **Key Findings**

* **CD3D** is highest in **Tcells**.
* **MS4A1** is highest in **Bcells**.
* **LYZ** is highest in **Myeloid cells**.
* **MKI67** expression varies, showing differential proliferation.
  → Marker expression patterns **support the assigned cluster labels**.

---

## **Files Included**

* `Single-Cell Marker Analysis.R` – Base R script containing all analysis code and comments.
* `Single-Cell Marker Analysis.docx` – Detailed report of the project.

---

## **How to Run**

1. Open the R script in RStudio.
2. Run each code chunk sequentially.
3. View output tables and base R plots directly in the console or plot pane.

---

