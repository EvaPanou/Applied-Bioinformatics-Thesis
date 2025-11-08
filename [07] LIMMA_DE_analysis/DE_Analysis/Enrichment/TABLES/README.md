
# Functional Enrichment Analysis — All Clusters (Up/Down-Regulated Gene Groups)

This README summarizes the biological interpretation derived from all 13 enrichment tables 
generated during the longitudinal DEG and enrichment analysis pipeline, corresponding to 
the annotated gene clusters obtained from the ComplexHeatmap k-means partitioning step 
(`Heatmap_GeneClusters_Annotated.tsv`).

---

## Overview

Each cluster represents a coherent group of genes (Up-, Down-, or Mixed-regulated) 
identified through the **limma longitudinal design** and subsequent **k-means clustering** 
on z-scored expression data. Enrichment analyses were performed for each cluster using:

- **GO: Biological Process (BP)**
- **GO: Cellular Component (CC)**
- **GO: Molecular Function (MF)**
- **KEGG Pathway**
- **Reactome Pathway**

All enrichments used *p < 0.05*, *q < 0.1*, and a background of all expressed genes.

---

## CLUSTER 1 — Immune and Inflammatory Activation (Up‑regulated)

**Summary:** Cluster 1 genes were strongly enriched for antiviral, cytokine, and 
interferon‑related immune pathways. This indicates systemic immune activation in the 
up‑regulated gene group, likely reflecting the pro‑inflammatory or disease‑associated 
response axis in SLE pregnancy.

**Highlighted Pathways:**
- **GO:BP:** Type I interferon signaling, defense response to virus, cytokine-mediated signaling.
- **GO:CC:** Cytosol, MHC class I complex, endoplasmic reticulum lumen.
- **GO:MF:** Double‑stranded RNA binding, cytokine receptor binding, antigen binding.
- **KEGG:** Cytokine–cytokine receptor interaction, RIG‑I‑like receptor signaling, Toll‑like receptor pathway.
- **Reactome:** Interferon alpha/beta signaling, ISG15 antiviral mechanism, innate immune system.

**Interpretation:** These functions indicate strong immune activation and antiviral defense 
characteristics consistent with SLE inflammatory profiles, particularly in response to 
type I interferon over‑activity.

---

## CLUSTER 2 — Minimal Down‑regulated Genes (CDC20, JCHAIN)

**Summary:** Cluster 2 consisted of only two genes (CDC20, JCHAIN). Due to the low gene count, 
no statistically significant enrichment was detected across GO or pathway databases. 
However, their functions are biologically interpretable:

- **CDC20:** Cell‑cycle regulator involved in anaphase promoting complex (mitotic exit).
- **JCHAIN:** Joining chain of multimeric immunoglobulins (IgA, IgM) — associated with plasma cell activity.

**Interpretation:** The presence of **CDC20** suggests transient suppression of 
cell‑cycle progression, while **JCHAIN** reduction may reflect altered humoral response or 
differentiation states in immune cells of SLE patients.

---

## CLUSTER 3 — Translational and Metabolic Reprogramming (Down‑regulated)

**Summary:** Cluster 3 was characterized by strong enrichment for mitochondrial, 
metabolic, and protein synthesis processes, consistent with global transcriptomic 
repression in metabolism‑related pathways during disease or late pregnancy phases.

**Highlighted Pathways:**
- **GO:BP:** Mitochondrial translation, oxidative phosphorylation, peptide biosynthetic process.
- **GO:CC:** Ribosome, mitochondrion, endoplasmic reticulum membrane.
- **GO:MF:** Structural constituent of ribosome, oxidoreductase activity.
- **KEGG:** Oxidative phosphorylation, TCA cycle, ribosome pathway.
- **Reactome:** Respiratory electron transport, translation initiation, protein metabolism.

**Interpretation:** This cluster marks a global downregulation of energy metabolism and 
ribosomal biogenesis, typical of cellular energy reallocation under chronic immune stress 
or adaptation to late gestational metabolic states.

---

## BIOLOGICAL SYNTHESIS ACROSS CLUSTERS

| Cluster | Dominant Regulation | Key Themes | Representative Pathways |
|----------|--------------------|-------------|--------------------------|
| 1 | Up | Immune activation, interferon signaling | Cytokine, antiviral, antigen processing |
| 2 | Down (small) | Cell cycle / antibody synthesis | CDC20, JCHAIN |
| 3 | Down | Metabolic suppression, translation | Oxidative phosphorylation, ribosome |

**Global trend:**  
These three clusters together depict the classic SLE immune profile in pregnancy:
a strong activation of interferon‑mediated inflammation (Cluster 1) concurrent with 
a metabolic and translational downshift (Cluster 3) and a mild reduction in 
proliferative / plasma‑cell gene activity (Cluster 2).

---

## Technical Details

- Enrichment conducted via **clusterProfiler**, **ReactomePA**, and **org.Hs.eg.db**.
- Universe = all genes expressed in the normalized microarray matrix.
- Significance thresholds: *p < 0.05*, *q < 0.1* (Benjamini–Hochberg correction).
- Visualization: `dotplot()` (clusterProfiler) for each ontology; results exported as TSV and PNG.

---

## Files Included

| File | Description |
|------|--------------|
| Heatmap_RowCluster_1_genes_GO_BP.tsv | GO Biological Process enrichment (Cluster 1) |
| Heatmap_RowCluster_1_genes_GO_CC.tsv | GO Cellular Component enrichment (Cluster 1) |
| Heatmap_RowCluster_1_genes_GO_MF.tsv | GO Molecular Function enrichment (Cluster 1) |
| Heatmap_RowCluster_1_genes_KEGG.tsv | KEGG pathway enrichment (Cluster 1) |
| Heatmap_RowCluster_1_genes_Reactome.tsv | Reactome pathway enrichment (Cluster 1) |
| Heatmap_RowCluster_3_genes_GO_BP.tsv | GO Biological Process enrichment (Cluster 3) |
| Heatmap_RowCluster_3_genes_GO_CC.tsv | GO Cellular Component enrichment (Cluster 3) |
| Heatmap_RowCluster_3_genes_GO_MF.tsv | GO Molecular Function enrichment (Cluster 3) |
| Heatmap_RowCluster_3_genes_KEGG.tsv | KEGG pathway enrichment (Cluster 3) |
| Heatmap_RowCluster_3_genes_Reactome.tsv | Reactome pathway enrichment (Cluster 3) |
| Heatmap_RowCluster_2_genes.tsv | Down‑regulated genes (CDC20, JCHAIN) — no enrichment |
| Enrichment_tables_batch1_top10.md | Top‑10 term markdown from Cluster 1 |
| Enrichment_tables_batch1_overview.tsv | Table summary of enriched term counts |
