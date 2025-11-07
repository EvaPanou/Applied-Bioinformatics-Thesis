# Folder: `A_SLE_vs_Healthy/` — What these files mean and how to read them

This folder contains **per-timepoint differential expression** results from the model:

- **Design:** a no-intercept model with one column per (Condition × Time) cell.

- **Contrast family A:** *SLE vs Healthy* performed **separately at each timepoint** (<16, 16–23, 24–31, 32–40 weeks, and Postpartum).

- **Blocking:** repeated measures handled by `duplicateCorrelation` using `Donor_id`.



Each `A_FULL_...tsv` file is the **complete limma topTable** for that single contrast (no filtering). Columns:

- `logFC` → SLE minus Healthy (positive = **higher in SLE**).

- `AveExpr` → average expression across samples.

- `t`, `P.Value` → moderated t-statistic and its raw p-value.

- `adj.P.Val` → FDR (Benjamini–Hochberg).

- `B` → log-odds of differential expression.

- `Gene` → gene symbol.



## Quick summary across timepoints

| Contrast                          |   Genes tested |   UP @ FDR 0.05 (logFC>1) |   DOWN @ FDR 0.05 (logFC<-1) |   UP @ FDR 0.01 (logFC>1) |   DOWN @ FDR 0.01 (logFC<-1) |
|:----------------------------------|---------------:|--------------------------:|-----------------------------:|--------------------------:|-----------------------------:|
| SLE vs Healthy at <16 weeks       |           9037 |                        49 |                            1 |                        49 |                            1 |
| SLE vs Healthy at 16–23 weeks     |           9037 |                        49 |                            1 |                        49 |                            1 |
| SLE vs Healthy at 24–31 weeks     |           9037 |                        44 |                            0 |                        44 |                            0 |
| SLE vs Healthy at 32–40 weeks     |           9037 |                        32 |                            0 |                        32 |                            0 |
| SLE vs Healthy at Postpartum (PP) |           9037 |                        40 |                            0 |                        40 |                            0 |


> **How to read this table**  

- We apply your thresholds: FDR ∈ {0.05, 0.01} and |logFC| > 1 (≈2× change).  

- **UP** means over-expressed in SLE vs Healthy; **DOWN** means under-expressed in SLE.  



### Take‑home pattern

- A **robust SLE up‑regulation signature** is present at *all* timepoints, with **very few or zero down‑regulated genes** under the strict |logFC|>1 threshold.

- The counts are **identical at FDR 0.05 and 0.01**, implying these SLE‑up genes are **highly significant**.



## The biological signal you’re seeing (beginner‑friendly)

When we look at the **top genes** by p‑value in each file, we repeatedly see classic **Type I interferon‑stimulated genes (ISGs)**. Examples include:

**SLE vs Healthy at <16 weeks**

- IFI44L (logFC 2.96, adj.P.Val 2.1e-32)
- IFI27 (logFC 4.04, adj.P.Val 2.8e-31)
- EPSTI1 (logFC 2.22, adj.P.Val 1.6e-29)
- RSAD2 (logFC 3.28, adj.P.Val 4.2e-29)
- LY6E (logFC 2.26, adj.P.Val 7.9e-28)
- ISG15 (logFC 2.45, adj.P.Val 8.7e-28)
- XAF1 (logFC 1.78, adj.P.Val 8.7e-28)
- OASL (logFC 1.87, adj.P.Val 1.4e-27)
- HERC5 (logFC 2.03, adj.P.Val 5.8e-27)
- IFI6 (logFC 1.65, adj.P.Val 7.1e-27)

**SLE vs Healthy at 16–23 weeks**

- IFI44L (logFC 2.99, adj.P.Val 1.2e-31)
- EPSTI1 (logFC 2.32, adj.P.Val 1.1e-30)
- LY6E (logFC 2.45, adj.P.Val 1.1e-30)
- OASL (logFC 2.02, adj.P.Val 4.2e-30)
- IFI6 (logFC 1.81, adj.P.Val 4.2e-30)
- ISG15 (logFC 2.61, adj.P.Val 4.2e-30)
- MX1 (logFC 2.02, adj.P.Val 1.4e-29)
- OAS1 (logFC 2.23, adj.P.Val 1.4e-29)
- HERC5 (logFC 2.18, adj.P.Val 1.5e-29)
- RSAD2 (logFC 3.30, adj.P.Val 9.7e-29)

**SLE vs Healthy at 24–31 weeks**

- OASL (logFC 1.93, adj.P.Val 1.0e-26)
- ISG15 (logFC 2.39, adj.P.Val 1.0e-24)
- LY6E (logFC 2.17, adj.P.Val 4.5e-24)
- IFI44L (logFC 2.54, adj.P.Val 1.8e-23)
- IFI6 (logFC 1.59, adj.P.Val 2.0e-23)
- OAS1 (logFC 1.99, adj.P.Val 2.0e-23)
- OAS3 (logFC 1.70, adj.P.Val 2.0e-23)
- EPSTI1 (logFC 2.00, adj.P.Val 3.2e-23)
- RSAD2 (logFC 2.97, adj.P.Val 3.6e-23)
- HERC5 (logFC 1.93, adj.P.Val 3.9e-23)

**SLE vs Healthy at 32–40 weeks**

- ZBP1 (logFC 1.19, adj.P.Val 7.3e-16)
- ISG15 (logFC 2.03, adj.P.Val 7.3e-16)
- RSAD2 (logFC 2.50, adj.P.Val 2.5e-14)
- OAS3 (logFC 1.40, adj.P.Val 5.0e-14)
- HERC5 (logFC 1.60, adj.P.Val 5.0e-14)
- IFI44L (logFC 2.07, adj.P.Val 5.0e-14)
- OASL (logFC 1.45, adj.P.Val 5.1e-14)
- HERC6 (logFC 1.23, adj.P.Val 7.3e-14)
- IFI27 (logFC 2.79, adj.P.Val 3.8e-13)
- DHX58 (logFC 1.16, adj.P.Val 4.5e-13)

**SLE vs Healthy at Postpartum (PP)**

- IFI27 (logFC 3.48, adj.P.Val 7.7e-19)
- XAF1 (logFC 1.56, adj.P.Val 6.1e-17)
- IFITM3 (logFC 1.41, adj.P.Val 1.4e-16)
- EPSTI1 (logFC 1.81, adj.P.Val 3.1e-16)
- IFI44L (logFC 2.27, adj.P.Val 3.2e-16)
- EIF2AK2 (logFC 1.26, adj.P.Val 4.3e-16)
- OASL (logFC 1.56, adj.P.Val 9.5e-16)
- ISG15 (logFC 1.96, adj.P.Val 1.0e-14)
- RSAD2 (logFC 2.55, adj.P.Val 1.1e-14)
- LY6E (logFC 1.79, adj.P.Val 1.7e-14)


**Interpretation in one sentence:** SLE samples show a strong, recurring **interferon/antiviral response program** relative to Healthy controls at every pregnancy stage and postpartum.


## Academic conclusions you can state from these files

1. **Consistency across gestation:** The SLE‑associated transcriptional program is detectable at <16 weeks and persists through 32–40 weeks and postpartum.  

2. **Directionality:** The signature is **dominated by up‑regulation in SLE** (virtually no strong down‑regulation at |logFC|>1), indicating an activated state rather than suppression.  

3. **Mechanistic theme:** The leading genes (e.g., *IFI44L, IFI27, IFIT3, RSAD2, ISG15, MX1, OASL*) point to **Type I interferon signaling** and antiviral defense pathways as key axes differentiating SLE from Healthy at each timepoint.  

4. **Statistical strength:** The fact that your UP counts are the same at FDR 0.05 and 0.01 shows the effect is **not marginal**—these genes pass very stringent multiple‑testing correction.  

5. **Clinical angle (hypothesis‑level):** A stable interferon‑high transcriptomic state through pregnancy could relate to SLE disease activity and risk; this justifies downstream **pathway enrichment** and **predictive modeling** using these features.


## How to use these outputs downstream

- Use the **`_UP_FDR0.05_` / `_UP_FDR0.01_`** files (same folder) when you want the filtered gene lists for plots or enrichment.

- Feed the **union of A‑contrast DEGs** into your ML step (your script already exports `ML_Features_AllDEGs.tsv`).

- Run GO/KEGG/Reactome enrichment on these UP genes; you should recover interferon and antiviral pathways prominently.

- Compare across timepoints with **UpSet plots** to document the stability of the signature.


## Caveats / good‑practice notes

- **Effect‑size threshold:** |logFC|>1 is strict; relaxing to 0.5 can reveal additional, biologically relevant genes while keeping FDR control.

- **Sample balance:** Always check counts per (Condition × Time) cell to ensure contrasts aren’t driven by sparse groups.

- **Blocking:** You accounted for repeated measures via `duplicateCorrelation`, which is the right approach for longitudinal donors.
