# Enrichment Analysis — Dotplot Visualization (All Clusters)

*Generated automatically using R (clusterProfiler v4.10.0, Reactom (clusterProfiler v4.10.0, ReactomePA v1.42.0).
This Folder presents the **functional enrichment visualization** results obtained using `clusterProfiler` and `ReactomePA` for each expression cluster identified in the SLE pregnancy dataset. Each dotplot highlights overrepresented Gene Ontology (GO) terms and pathways across five ontologies — **GO Biological Process (BP)**, **GO Cellular Component (CC)**, **GO Molecular Function (MF)**, **KEGG**, and **Reactome**.
All clusters were derived from hierarchical clustering of normalized expression data (Row-clustered heatmap). Functional annotation was performed to identify the predominant biological themes contributing to the molecular phenotype of each cluster.

---

### Summary

- **Cluster 1:** Represents **innate and neutrophil-driven inflammatory activation**, linking lupus flares to infection-like signatures.  
- **Cluster 2:** Minimal enrichment — reflects isolated genes without clear pathway co-regulation.  
- **Cluster 3:** Exhibits a **strong antiviral and interferon-stimulated profile**, the transcriptional hallmark of SLE pathogenesis, relevant to pregnancy immune dysregulation.

Together, these clusters outline the **dual immune axis** active in SLE pregnancy: an inflammatory/neutrophilic arm (Cluster 1) and an antiviral/interferon-mediated arm (Cluster 3), emphasizing how maternal immune adaptation and autoimmunity intersect.

---

## 📁 File Index

| File Name | Ontology | Cluster | Description |
|------------|-----------|----------|--------------|
| Heatmap_RowCluster_1_genes_GO_BP_dotplot.png | GO: Biological Process | Cluster 1 | Top enriched biological processes (immune & inflammatory signaling) |
| Heatmap_RowCluster_1_genes_GO_CC_dotplot.png | GO: Cellular Component | Cluster 1 | Cellular compartments enriched in Cluster 1 genes |
| Heatmap_RowCluster_1_genes_GO_MF_dotplot.png | GO: Molecular Function | Cluster 1 | Enriched molecular activities of Cluster 1 proteins |
| Heatmap_RowCluster_1_genes_KEGG_dotplot.png | KEGG Pathway | Cluster 1 | Overrepresented KEGG pathways (infection and immune-related) |
| Heatmap_RowCluster_1_genes_Reactome_dotplot.png | Reactome Pathway | Cluster 1 | Reactome pathways enriched in Cluster 1 |
| Heatmap_RowCluster_3_genes_GO_BP_dotplot.png | GO: Biological Process | Cluster 3 | Top enriched biological processes (antiviral and interferon-related) |
| Heatmap_RowCluster_3_genes_GO_CC_dotplot.png | GO: Cellular Component | Cluster 3 | Main subcellular localizations of Cluster 3 gene products |
| Heatmap_RowCluster_3_genes_GO_MF_dotplot.png | GO: Molecular Function | Cluster 3 | Enriched molecular functions in Cluster 3 |
| Heatmap_RowCluster_3_genes_KEGG_dotplot.png | KEGG Pathway | Cluster 3 | KEGG pathways linked to viral infection and immune response |
| Heatmap_RowCluster_3_genes_Reactome_dotplot.png | Reactome Pathway | Cluster 3 | Reactome modules enriched in Cluster 3 genes |

---

## Cluster 1 — Inflammatory & Neutrophil-Driven Responses

### GO Biological Process (GO_BP)
**File:** `Heatmap_RowCluster_1_genes_GO_BP_dotplot.png`  
This dotplot reveals strong enrichment in **antibacterial humoral response**, **cell killing**, and **response to lipopolysaccharide**, reflecting a robust **innate immune activation** and **inflammatory cell recruitment**. The high GeneRatio and low adjusted p-values indicate that these genes participate extensively in **defense mechanisms**, **neutrophil-mediated immunity**, and **lipid response signaling**, consistent with heightened inflammation.

**Interpretation:** These functions suggest a **pro-inflammatory activation state**, possibly representing the **innate immune component** of SLE pathophysiology, which becomes further intensified during pregnancy where immune modulation is delicate.

---

### GO Cellular Component (GO_CC)
**File:** `Heatmap_RowCluster_1_genes_GO_CC_dotplot.png`  
Enrichment in **azurophil granule** and **primary lysosome** indicates localization of Cluster 1 proteins in **neutrophil secretory granules** and **phagolysosomal compartments**. This pattern suggests active **degranulation** and **enzymatic release** processes during immune activation.

**Interpretation:** These subcellular localizations align with **neutrophil extracellular trap (NET) formation**, a hallmark mechanism implicated in SLE flare-ups and vascular inflammation, which may contribute to adverse pregnancy outcomes.

---

### GO Molecular Function (GO_MF)
**File:** `Heatmap_RowCluster_1_genes_GO_MF_dotplot.png`  
Highlighted molecular functions include **serine-type endopeptidase activity**, **protease binding**, and **heparin binding**. These activities correspond to **proteolytic enzymes** involved in extracellular matrix degradation and immune cell migration.

**Interpretation:** Upregulation of these enzymatic functions reinforces the concept of **tissue remodeling** and **leukocyte infiltration**, integral to both normal pregnancy adaptation and SLE-associated inflammation.

---

### KEGG Pathway
**File:** `Heatmap_RowCluster_1_genes_KEGG_dotplot.png`  
Top enriched KEGG pathways include **Systemic lupus erythematosus**, **NOD-like receptor signaling**, and **Staphylococcus aureus infection**, reflecting strong overlaps between **autoimmune** and **infectious-like inflammatory pathways**.

**Interpretation:** Cluster 1 represents a transcriptional signature consistent with **autoimmune activation** and **innate immune dysregulation**, bridging microbial response pathways and lupus-specific mechanisms.

---

### Reactome Pathway
**File:** `Heatmap_RowCluster_1_genes_Reactome_dotplot.png`  
Enriched pathways such as **Neutrophil degranulation** and **Antimicrobial peptide production** underline the **innate effector arm of inflammation**. Reactome terms like *Defensins* and *Extracellular matrix organization* suggest an active tissue defense response.

**Interpretation:** This pattern likely corresponds to the **effector phase of inflammation**, dominated by neutrophil activation and tissue turnover — a recurrent SLE hallmark that can be exacerbated under pregnancy-associated immune tolerance shifts.

---

## Cluster 2 — Minimal Enrichment Signature

**File:** `Heatmap_RowCluster_2_Genes.tsv`  
Cluster 2 contained only **two genes (CDC20, JCHAIN)**, which was **insufficient for enrichment analysis**.  
- *CDC20* is linked to **cell cycle progression**, while *JCHAIN* encodes the **joining chain of multimeric immunoglobulins** (IgA/IgM).  
Together, these genes suggest a **limited regulatory signature**, potentially representing a **transitional or low-variance cluster** without robust pathway-level coherence.

---

## Cluster 3 — Antiviral & Interferon-Mediated Immunity

### GO Biological Process (GO_BP)
**File:** `Heatmap_RowCluster_3_genes_GO_BP_dotplot.png`  
This plot is dominated by **response to virus**, **type I interferon signaling**, and **defense response to biotic stimulus**, all key elements of **antiviral immunity** and **cytokine-mediated defense**.

**Interpretation:** The interferon-centered enrichment reflects the **canonical SLE interferon signature**, which remains active even during pregnancy. This cluster likely represents **plasmacytoid dendritic cell activity** or interferon-stimulated gene expression contributing to chronic immune stimulation.

---

### GO Cellular Component (GO_CC)
**File:** `Heatmap_RowCluster_3_genes_GO_CC_dotplot.png`  
Only *P granule* and *germ plasm* are enriched, pointing to **RNA-protein complexes** involved in **RNA metabolism** and **post-transcriptional regulation**.

**Interpretation:** This localization may correspond to **interferon-induced RNA-binding proteins** or **stress granule components**, consistent with antiviral signaling cascades activated in SLE.

---

### GO Molecular Function (GO_MF)
**File:** `Heatmap_RowCluster_3_genes_GO_MF_dotplot.png`  
Functional enrichment includes **RNA helicase activity**, **double-stranded RNA binding**, and **ATP-dependent chromatin remodeling**, all of which play roles in **viral RNA recognition** and **immune gene regulation**.

**Interpretation:** These functions underscore **innate immune sensing** mechanisms — particularly **RIG-I-like helicases (DDX58, IFIH1)** — that trigger type I interferon production in response to viral or self-derived nucleic acids, a central mechanism in SLE immunopathogenesis.

---

### KEGG Pathway
**File:** `Heatmap_RowCluster_3_genes_KEGG_dotplot.png`  
Enriched KEGG terms include **Hepatitis C**, **Influenza A**, **COVID-19**, and **Herpesvirus infection**, all converging on **viral recognition and interferon activation pathways**.

**Interpretation:** Although labeled as infection pathways, these represent **pattern-recognition signaling routes** reused by the immune system during lupus flare states, mimicking antiviral responses in absence of infection — a phenomenon well-described in SLE molecular pathology.

---

### Reactome Pathway
**File:** `Heatmap_RowCluster_3_genes_Reactome_dotplot.png`  
Major terms include **Interferon alpha/beta signaling**, **ISG15 antiviral mechanism**, and **TRAF6-mediated IRF7 activation** — all hallmark steps of **innate immune sensing and interferon-stimulated response**.

**Interpretation:** This Reactome profile consolidates the **type I interferon axis** as the defining molecular signature of Cluster 3. During pregnancy, sustained activation of these pathways could contribute to **placental inflammation**, **fetal growth restriction**, or **disease exacerbation** in SLE patients.

---

## Figure Legend

Each dotplot represents the **top enriched terms** for a given ontology:
- **x-axis (GeneRatio):** Fraction of cluster genes annotated to a given term.  
- **y-axis:** Enriched GO term or pathway name.  
- **Bubble size:** Number of genes contributing to the term.  
- **Color scale (p.adjust):** Adjusted p-value (Benjamini–Hochberg), with red indicating higher significance.  

All plots were generated using the `clusterProfiler` (GO, KEGG) and `ReactomePA` (Reactome) packages in R.  
Significant terms were defined at *p.adjust < 0.05*.
