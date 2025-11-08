# Visualisation of the Differential Expression Analysis in SLE Pregnancy

This section summarizes and interprets the principal transcriptomic patterns visualized through the general plots and tables generated from the longitudinal differential expression (DE) analysis of **Systemic Lupus Erythematosus (SLE) pregnancies** versus healthy pregnancies. Using the *limma* framework, gene expression levels were modeled as a function of **Condition (SLE vs. Healthy)**, **Time** (five gestational stages and postpartum), and their interaction, incorporating within-subject correlation to account for repeated measures. The resulting figures illustrate how gene expression trajectories evolve across gestation and postpartum, highlighting both stable and time-dependent disease-associated modules.

---

## 1. Global Expression Structure (Heatmap Analysis)

*Figure: Heatmap_A_DEG_Union_FDR0.05.png*  

The heatmap visualizes the global transcriptional configuration of all genes differentially expressed (FDR < 0.05) between SLE and healthy pregnancies across timepoints. Rows correspond to genes, columns to samples arranged by Condition (SLE → Healthy) and gestational progression (<16 weeks → postpartum). Color intensity represents standardized expression (z-score per gene).  

It displays an overarching **bimodal organization**, where most SLE samples cluster toward **sustained upregulation** (red shading), particularly during mid-to-late gestation, while healthy samples show reciprocal **downregulation** (blue shading) of the same gene set. This pattern reflects a pervasive **pro-inflammatory and interferon-stimulated transcriptional signature**, consistent with the characteristic SLE molecular phenotype described in prior studies.  

In the heatmap, gene rows were partitioned into three major clusters—**Up-regulated genes**, **Mixed/time-dependent genes**, and **Down-regulated genes**—using k-means clustering (k = 3). The “Up-regulated” cluster corresponds to innate immune activation, dominated by interferon-stimulated genes (ISGs) such as *IFI44L*, *MX1*, and *OAS1*, while the “Down-regulated” cluster includes genes involved in lymphocyte signaling and metabolic maintenance. The “Mixed/time-dependent” group shows transient modulation around the late gestational and postpartum stages, reflecting physiological rebalancing mechanisms characteristic of maternal immune tolerance.  

Collectively, it reveals reveal that the **SLE transcriptional landscape during pregnancy is marked by sustained immune activation and incomplete postpartum resolution**, pointing to a maladaptive persistence of inflammatory signals.

---

## 2. Intersecting Differential Gene Sets Across Contrasts (UpSet Analysis — Program A)

*Figure: UpSet_A_UP_FDR0.05.png*  
*Figure: UpSet_A_DN_FDR0.05.png*  

The UpSet plots summarize intersections of DE gene lists across multiple pairwise contrasts (SLE vs. Healthy at each timepoint). In **UpSet_A_UP_FDR0.05.png**, the largest intersection, comprising approximately 30 genes, represents a **core set of genes consistently upregulated in SLE** throughout gestation. Functional annotation of this core (from enrichment analysis) highlights **type I interferon signaling**, **viral defense responses**, and **neutrophil activation**, reinforcing previous reports of sustained ISG expression in active SLE and during complicated pregnancies.  

Conversely, **UpSet_A_DN_FDR0.05.png** shows sparse and heterogeneous intersections, indicating that **downregulation events are more timepoint-specific** and transient. These patterns may reflect context-dependent modulation of adaptive or metabolic pathways, which are temporarily suppressed in response to fluctuating immune or hormonal states. The imbalance between broad recurrent upregulation and limited downregulation supports the model of **persistent innate immune hyperactivation** in SLE pregnancies.

---

## 3. Temporal Expression Dynamics Within Conditions (UpSet Analysis — Program B)

*Figure: UpSet_B_UP_FDR0.05.png*  
*Figure: UpSet_B_DN_FDR0.05.png*  

The second pair of UpSet plots illustrates **temporal differential expression within each condition**, comparing adjacent gestational intervals (“adjacent” mode). In **UpSet_B_DN_FDR0.05.png**, large intersections of downregulated genes are evident between late gestation and postpartum stages, in both SLE and healthy cohorts. This widespread downregulation likely corresponds to a **global postpartum immune contraction**, as maternal physiology transitions back from a pregnancy-adapted immune state. The magnitude of this contraction, however, appears more pronounced in the SLE group, suggesting an **overcorrection or incomplete normalization** of immune gene expression.  

In contrast, **UpSet_B_UP_FDR0.05.png** reveals smaller, condition-specific intersections of upregulated genes, reflecting **localized activation events** during gestation. These may capture bursts of transcriptional activity related to **placental signaling**, **cytokine-driven inflammation**, or **stress responses**, which occur sporadically across trimesters. Together, the Program B plots reveal that while healthy pregnancies exhibit a controlled cyclical modulation of immune genes, SLE pregnancies retain aberrant activation signatures that extend beyond the expected postpartum retraction phase.


## 4. Integrative Discussion

Taken together, these visualizations delineate a **dual-axis immune disturbance** underlying SLE pregnancies. The first axis represents a **stable, interferon-dominated innate immune signature**, present from early gestation through postpartum, characterized by antiviral and pro-inflammatory transcripts. The second axis captures **time-dependent modulation**, reflecting the partial adaptation of immune expression to physiological gestational demands.  

While healthy pregnancies exhibit predictable immunological transitions—attenuation of inflammatory pathways in mid-gestation followed by normalization postpartum—SLE pregnancies demonstrate **persistent upregulation** of both interferon and neutrophil-related programs, coupled with incomplete restoration of adaptive gene expression after delivery. These findings align with the broader literature on immune dysregulation in SLE pregnancies, and provide transcriptomic evidence for a **blunted immune rebalancing mechanism** in affected mothers.  

In summary, the collective data presented in these figures support the conclusion that **SLE pregnancies are characterized by sustained innate immune activation, insufficient adaptive compensation, and delayed postpartum recovery**, reinforcing the notion of chronic transcriptional disequilibrium within the maternal immune system.
