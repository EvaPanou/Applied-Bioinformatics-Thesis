# 02_DGE_Analysis / 02_Results / Gene_Clusters

## Folder Overview

Method A's 65-gene DEG union (`Method_A_Union`, from `../A_SLE_vs_Healthy`), split by k-means clustering on the z-scored heatmap data (`../Heatmap_A_DEG_Union_FDR0.05.png`) into up to 3 groups by expression pattern. Produced within Step 7 of `../../01_limma_DEG_pipeline.R`. Categorized raw gene lists for manual review, not an automated enrichment result.

## Folder Structure & File Reference

| File | Genes | Description |
|---|---|---|
| `Heatmap_Cluster_Up_regulated_in_SLE_genes.tsv` | 65 | Every gene in Method A's DEG union — identical, gene-for-gene, to `../DE_Genes/Method_A_Union.tsv`, since all 65 genes landed in this one cluster. |

**Only one file exists here — this is expected, not a gap.** The clustering code always runs `kmeans(..., centers = 3)`, producing exactly 3 clusters, and labels each cluster "Up-regulated in SLE" / "Down-regulated in SLE" / "Mixed/time-dependent genes" by comparing its mean SLE-column expression to its mean Healthy-column expression (thresholds: >0.15 / <-0.15 / between). A file is only written per label if at least one cluster actually received that label — confirmed directly against a screenshot of this folder in Eva's actual local file explorer, showing the same single file. **All 3 of the k-means clusters landed on "Up-regulated in SLE"** for this 65-gene set; none of the 65 genes pulled a cluster's average past the -0.15 threshold in the opposite direction. This is consistent with the rest of this stage's findings: only one gene (`ORM1`) shows any down-regulation anywhere in the DEG union, which isn't enough on its own to pull a whole k-means cluster's average into "Down-regulated" or even "Mixed" territory.

---

## Results

The file is a single-column gene list, identical to `Method_A_Union.tsv`:

```text
Gene
IFI44L
IFI27
EPSTI1
RSAD2
HERC5
...
ORM1
```

(65 genes total — the full list, including `ORM1` itself, since `ORM1`'s single down-regulated data point wasn't enough to break it out into its own cluster.)

**Conclusion:** this folder doesn't add new gene-selection information beyond `Method_A_Union` — its value is confirmatory: an unsupervised clustering step, run independently of the significance testing that produced the union in the first place, agrees that this gene set behaves as one coherent up-regulated block rather than splitting into meaningfully distinct sub-patterns.

**Not used downstream** by any later stage — this is an end-of-line manual-review artifact for the thesis discussion, not an input to `03_Feature_Validation` or `04_ML_Prerequisites` (those build on `DE_Genes/Method_All_Union`, not on this clustering).
