# ML Pipeline Flowcharts

This directory contains the methodology flowcharts supporting the donor-stratified machine learning pipeline used to identify and evaluate a compact diagnostic gene signature for systemic lupus erythematosus (SLE).

The diagrams provide a visual summary of the main analytical workflow, from the prepared DEG-based input matrix and donor-level data split, through algorithm benchmarking and final gene-panel definition, to sealed held-out test-set evaluation.

---

## Folder Structure

```text
[04] Flowchart/
├── `README.md`
├── `Overall_Pipeline_Flowchart.pdf` = includes the whole pipeline methodology in pdf format for greater detail
├── 
`Overall_Pipeline_Flowchart.png` = includes the whole pipeline methodology in png format
├── `Flowchart_Section1.png` = includes only Section 1 "Data Input & Donor Split"
├── `Flowchart_Section2.png` = includes only Section 2 "Algorithm Benchmarking on Development Set"
├── `Flowchart_Section3.png` = includes only Section 3 "Final gene-signature panel definition"
└── `Flowchart_Section4.png` = includes only Section 4 "Held-out test set evaluation"