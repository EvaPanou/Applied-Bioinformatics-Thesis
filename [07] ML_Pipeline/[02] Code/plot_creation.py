#-----------------------------------------------------------------------
# PLOT CREATION SCRIPT
#-----------------------------------------------------------------------

# THIS WILL CONTAIN ALL THE PLOT CREATION FUNCTIONS FOR THE MAIN SLE DIAGNOSTIC PIPELINE SCRIPT
# It involves only figure-generation functions.
# The idea is to keep the main pipeline and helper_functions.py cleaner
# by moving all plots into one dedicated script.
# IT WILL NEED TO BE LOCATED IN THE SAME DIRECTORY AS THE MAIN PIPELINE SCRIPT "sle_pipeline.py"
# in order to be direclty called into action


from __future__ import annotations

from pathlib import Path
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import confusion_matrix

from configurations import CB_BLUE, CB_ORANGE, CB_GREY
from requirements import HAS_SHAP

# Optional plotting packages
try:
    import seaborn as sns
except Exception:
    sns = None

try:
    import shap
except Exception:
    shap = None


# --------------------------------------------------------------------
# Small internal helper functions
# --------------------------------------------------------------------

def get_log(log: logging.Logger | None = None) -> logging.Logger:
    """
    Return the provided log object or fall back to the root logger.
    """
    return log or logging.getLogger()


def ensure_figures_dir(out_dir: Path) -> Path:
    """
    Create and return the figures subfolder inside the main output folder.
    """
    fig_dir = Path(out_dir) / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    return fig_dir


def save_figure(fig: plt.Figure, save_path: Path, log: logging.Logger | None = None) -> None:
    """
    Save a matplotlib figure, close it, and record the action in the log.
    """
    log = get_log(log)
    fig.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("saved plot -> %s", save_path)


# --------------------------------------------------------------------
# Benchmark plots
# --------------------------------------------------------------------

# Benchmark Boxplot: in order to compare all the FS-classifier combination's AUC
def plot_benchmark_boxplot(
    results_df: pd.DataFrame,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Boxplot of outer-fold AUC values for each FS/classifier combination.

    Expected columns:
    - fs_method
    - classifier
    - auc
    """
    log = get_log(log)

    if results_df is None or results_df.empty:
        log.warning("plot_benchmark_boxplot skipped: results_df is empty")
        return

    needed = {"fs_method", "classifier", "auc"}
    missing = needed - set(results_df.columns)
    if missing:
        log.warning("plot_benchmark_boxplot skipped: missing columns %s", sorted(missing))
        return

    df = results_df.copy()
    df["combo"] = df["fs_method"].astype(str) + " + " + df["classifier"].astype(str)

    combo_order = (
        df.groupby("combo")["auc"]
        .mean()
        .sort_values(ascending=False)
        .index
        .tolist()
    )

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(max(10, 0.9 * len(combo_order)), 6))

    if sns is not None:
        sns.boxplot(
            data=df,
            x="combo",
            y="auc",
            order=combo_order,
            orient="v",
            color="#d9ecf7",
            ax=ax,
        )
        sns.stripplot(
            data=df,
            x="combo",
            y="auc",
            order=combo_order,
            orient="v",
            color=CB_BLUE,
            alpha=0.75,
            size=4,
            ax=ax,
        )
    else:
        grouped = [df.loc[df["combo"] == c, "auc"].values for c in combo_order]
        ax.boxplot(grouped, tick_labels=combo_order, patch_artist=True)

    ax.set_title("Outer-fold benchmark AUC distributions")
    ax.set_xlabel("Feature selector + classifier")
    ax.set_ylabel("AUC")
    ax.set_ylim(0.8, 1.0)
    ax.tick_params(axis="x", rotation=45)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    save_figure(fig, fig_dir / "04_benchmark_boxplot.png", log)

# Benchmark BARplot: in order to compare all the FS-classifier combination results' mean and AUC
def plot_benchmark_barplot(
    summary: pd.DataFrame,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Barplot of mean benchmark AUC for each FS/classifier combination.

    Expected columns:
    - fs_method
    - classifier
    - mean
    - std (optional but expected)
    """
    log = get_log(log)

    if summary is None or summary.empty:
        log.warning("plot_benchmark_barplot skipped: summary is empty")
        return

    needed = {"fs_method", "classifier", "mean"}
    missing = needed - set(summary.columns)
    if missing:
        log.warning("plot_benchmark_barplot skipped: missing columns %s", sorted(missing))
        return

    df = summary.copy()
    df["combo"] = df["fs_method"].astype(str) + " + " + df["classifier"].astype(str)
    df = df.sort_values("mean", ascending=False).reset_index(drop=True)

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(max(10, 0.9 * len(df)), 6))

    x = np.arange(len(df))
    y = df["mean"].values
    yerr = df["std"].fillna(0).values if "std" in df.columns else None

    ax.bar(
        x,
        y,
        yerr=yerr,
        color=CB_BLUE,
        alpha=0.85,
        capsize=4,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(df["combo"], rotation=45, ha="right")
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Mean AUC")
    ax.set_xlabel("Feature selector + classifier")
    ax.set_title("Mean outer-fold benchmark AUC")
    ax.grid(axis="y", alpha=0.3)

    save_figure(fig, fig_dir / "04_benchmark_barplot.png", log)

# Benchmark Boxplot: in order to compare all the FS-classifier combination's mean AUC to find which is the best one
def plot_benchmark_heatmap(
    summary: pd.DataFrame,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Heatmap of mean AUC with rows = feature selectors and columns = classifiers.

    Expected columns:
    - fs_method
    - classifier
    - mean
    """
    log = get_log(log)

    if summary is None or summary.empty:
        log.warning("plot_benchmark_heatmap skipped: summary is empty")
        return

    needed = {"fs_method", "classifier", "mean"}
    missing = needed - set(summary.columns)
    if missing:
        log.warning("plot_benchmark_heatmap skipped: missing columns %s", sorted(missing))
        return

    heat = summary.pivot(index="fs_method", columns="classifier", values="mean")

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.8 * len(heat.index))))

    if sns is not None:
        sns.heatmap(
            heat,
            annot=True,
            fmt=".3f",
            cmap="Blues",
            vmin=0.0,
            vmax=1.0,
            linewidths=0.5,
            cbar_kws={"label": "Mean AUC"},
            ax=ax,
        )
    else:
        im = ax.imshow(heat.values, aspect="auto", cmap="Blues", vmin=0.0, vmax=1.0)
        ax.set_xticks(np.arange(len(heat.columns)))
        ax.set_xticklabels(heat.columns)
        ax.set_yticks(np.arange(len(heat.index)))
        ax.set_yticklabels(heat.index)

        for i in range(len(heat.index)):
            for j in range(len(heat.columns)):
                value = heat.iloc[i, j]
                if pd.notna(value):
                    ax.text(j, i, f"{value:.3f}", ha="center", va="center", fontsize=9)

        fig.colorbar(im, ax=ax, label="Mean AUC")

    ax.set_title("Benchmark heatmap")
    ax.set_xlabel("Classifier")
    ax.set_ylabel("Feature selector")

    save_figure(fig, fig_dir / "04_benchmark_heatmap.png", log)


# --------------------------------------------------------------------
# Final signature panel plots
# --------------------------------------------------------------------

# Elbow Curve: in order to find the best gene combo
def plot_elbow_curve(
    auc_curve,
    elbow_k: int,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Plot the inner-CV AUC curve used for choosing panel size.

    auc_curve is expected to correspond to k = 2, 3, 4, ...
    """
    log = get_log(log)

    if auc_curve is None or len(auc_curve) == 0:
        log.warning("plot_elbow_curve skipped: auc_curve is empty")
        return

    y = np.asarray(auc_curve, dtype=float)
    x = np.arange(2, 2 + len(y))

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.plot(x, y, marker="o", color=CB_BLUE, linewidth=2)
    ax.axvline(elbow_k, color=CB_ORANGE, linestyle="--", linewidth=2, label=f"elbow_k = {elbow_k}")

    if elbow_k in x:
        elbow_y = y[np.where(x == elbow_k)[0][0]]
        ax.scatter([elbow_k], [elbow_y], color=CB_ORANGE, s=80, zorder=3)

    ax.set_title("Elbow curve for panel-size selection")
    ax.set_xlabel("Number of genes kept (k)")
    ax.set_ylabel("Mean inner-CV AUC")
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.3)
    ax.legend()

    save_figure(fig, fig_dir / "03_elbow_curve.png", log)

# Stability barplot: in order to find the most common genes after bootstrap
def plot_stability_barplot(
    panel_table: pd.DataFrame,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Barplot of stability values for genes in the final selected panel.

    This function tries to detect the gene column and the stability column.
    Accepted gene column names:
    - gene, Gene, feature, Feature

    Accepted stability column names:
    - bootstrap_freq, stability, frequency, fold_frequency
    """
    log = get_log(log)

    if panel_table is None or panel_table.empty:
        log.warning("plot_stability_barplot skipped: panel_table is empty")
        return

    df = panel_table.copy()

    gene_col = None
    for c in ["gene", "Gene", "feature", "Feature"]:
        if c in df.columns:
            gene_col = c
            break

    value_col = None
    for c in ["bootstrap_freq", "stability", "frequency", "fold_frequency"]:
        if c in df.columns:
            value_col = c
            break

    if gene_col is None or value_col is None:
        log.warning(
            "plot_stability_barplot skipped: could not detect gene/stability columns in %s",
            list(df.columns),
        )
        return

    df = df.sort_values(value_col, ascending=False).reset_index(drop=True)

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(8, max(4, 0.45 * len(df))))

    ax.barh(df[gene_col], df[value_col], color=CB_BLUE, alpha=0.9)
    ax.invert_yaxis()
    ax.set_title("Final panel stability")
    ax.set_xlabel(value_col)
    ax.set_ylabel("Gene")
    ax.grid(axis="x", alpha=0.3)

    if df[value_col].max() <= 1.05:
        ax.set_xlim(0, 1)

    save_figure(fig, fig_dir / "03_panel_stability_barplot.png", log)

# Stability top 20: in order to find the 20 most common genes after bootstrap
def plot_full_stability_top20(
    stability,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Plot the top 20 most stable genes across the full finalisation process.

    Accepted input:
    - pandas Series with gene names in the index
    - pandas DataFrame with one gene column and one value column
    """
    log = get_log(log)

    if stability is None or len(stability) == 0:
        log.warning("plot_full_stability_top20 skipped: stability input is empty")
        return

    if isinstance(stability, pd.Series):
        df = stability.sort_values(ascending=False).head(20).reset_index()
        df.columns = ["gene", "stability"]
    else:
        df = pd.DataFrame(stability).copy()

        gene_col = None
        for c in ["gene", "Gene", "feature", "Feature"]:
            if c in df.columns:
                gene_col = c
                break
        if gene_col is None:
            gene_col = df.columns[0]

        value_col = None
        for c in ["bootstrap_freq", "stability", "frequency", "fold_frequency"]:
            if c in df.columns:
                value_col = c
                break
        if value_col is None:
            numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
            if not numeric_cols:
                log.warning("plot_full_stability_top20 skipped: no numeric stability column found")
                return
            value_col = numeric_cols[-1]

        df = df[[gene_col, value_col]].copy()
        df.columns = ["gene", "stability"]
        df = df.sort_values("stability", ascending=False).head(20)

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(8, 7))

    ax.barh(df["gene"], df["stability"], color=CB_ORANGE, alpha=0.9)
    ax.invert_yaxis()
    ax.set_title("Top 20 most stable genes")
    ax.set_xlabel("Stability")
    ax.set_ylabel("Gene")
    ax.grid(axis="x", alpha=0.3)

    if df["stability"].max() <= 1.05:
        ax.set_xlim(0, 1)

    save_figure(fig, fig_dir / "03_full_stability_top20.png", log)


def plot_final_panel_correlation_heatmap(
    X_dev: pd.DataFrame,
    panel_genes: list[str],
    out_dir: Path,
    logger: logging.Logger | None = None,
) -> None:
    log = logger or logging.getLogger()

    if len(panel_genes) < 2:
        log.warning("Skipping final panel correlation heatmap: fewer than 2 genes.")
        return

    corr_df = X_dev[panel_genes].corr(method="pearson")

    plt.figure(figsize=(max(6, 0.6 * len(panel_genes)), max(5, 0.6 * len(panel_genes))))
    sns.heatmap(
        corr_df,
        cmap="coolwarm",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.5,
        cbar_kws={"label": "Pearson r"},
    )
    plt.title("Final panel Pearson correlation heatmap")
    plt.tight_layout()
    plt.savefig(out_dir / "03_final_panel_correlation_heatmap.png", dpi=300, bbox_inches="tight")
    plt.close()

    log.info("saved -> 03_final_panel_correlation_heatmap.png")


# --------------------------------------------------------------------
# Held-out evaluation plots
# --------------------------------------------------------------------

# ROC Curve: in order to find the ROC for the hel-out TEST dataset of donors
def plot_roc_curve(
    holdout_metrics: dict,
    lr_holdout: dict | None,
    ci_info: dict | None,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Plot the main held-out ROC curve.

    Expected keys in holdout_metrics:
    - fpr
    - tpr
    - auc
    """
    log = get_log(log)

    if holdout_metrics is None:
        log.warning("plot_roc_curve skipped: holdout_metrics is None")
        return

    needed = {"fpr", "tpr", "auc"}
    missing = needed - set(holdout_metrics.keys())
    if missing:
        log.warning("plot_roc_curve skipped: missing keys %s", sorted(missing))
        return

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(6, 6))

    label = f"AUC = {holdout_metrics['auc']:.3f}"
    if ci_info is not None and {"ci_lo", "ci_hi"}.issubset(ci_info.keys()):
        label += f" (95% CI {ci_info['ci_lo']:.3f}-{ci_info['ci_hi']:.3f})"

    ax.plot(
        holdout_metrics["fpr"],
        holdout_metrics["tpr"],
        color=CB_BLUE,
        linewidth=2.5,
        label=label,
    )
    ax.plot([0, 1], [0, 1], linestyle="--", color=CB_GREY, label="Chance")

    ax.set_title("Held-out ROC curve")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")

    save_figure(fig, fig_dir / "04_roc_curve.png", log)

# ROC + lR: in order to compare the ROC with a LR one
def plot_roc_curve_with_lr(
    holdout_metrics: dict,
    lr_holdout: dict,
    ci_info: dict | None,
    lr_comparison: dict | None,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Plot the winning model ROC and a logistic-regression comparator ROC.

    Expected keys in both holdout_metrics and lr_holdout:
    - fpr
    - tpr
    - auc
    """
    log = get_log(log)

    if holdout_metrics is None or lr_holdout is None:
        log.warning("plot_roc_curve_with_lr skipped: missing holdout_metrics or lr_holdout")
        return

    for obj_name, obj in [("holdout_metrics", holdout_metrics), ("lr_holdout", lr_holdout)]:
        needed = {"fpr", "tpr", "auc"}
        missing = needed - set(obj.keys())
        if missing:
            log.warning("plot_roc_curve_with_lr skipped: %s missing %s", obj_name, sorted(missing))
            return

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(6, 6))

    main_label = f"Winning model AUC = {holdout_metrics['auc']:.3f}"
    if ci_info is not None and {"ci_lo", "ci_hi"}.issubset(ci_info.keys()):
        main_label += f" (95% CI {ci_info['ci_lo']:.3f}-{ci_info['ci_hi']:.3f})"

    lr_label = f"LR comparator AUC = {lr_holdout['auc']:.3f}"
    if lr_comparison is not None and {"ci_lo", "ci_hi"}.issubset(lr_comparison.keys()):
        lr_label += f" (95% CI {lr_comparison['ci_lo']:.3f}-{lr_comparison['ci_hi']:.3f})"

    ax.plot(
        holdout_metrics["fpr"],
        holdout_metrics["tpr"],
        color=CB_BLUE,
        linewidth=2.5,
        label=main_label,
    )
    ax.plot(
        lr_holdout["fpr"],
        lr_holdout["tpr"],
        color=CB_ORANGE,
        linewidth=2.2,
        label=lr_label,
    )
    ax.plot([0, 1], [0, 1], linestyle="--", color=CB_GREY, label="Chance")

    ax.set_title("Held-out ROC curve with LR comparator")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")

    save_figure(fig, fig_dir / "04_roc_curve_with_lr.png", log)

# Confusion Matrices: in order to visualise the True and False Positive and Negative Rates
def plot_confusion_matrices(
    holdout_metrics: dict,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Plot the held-out confusion matrix.

    Expected keys:
    - y_true
    - y_pred
    """
    log = get_log(log)

    if holdout_metrics is None:
        log.warning("plot_confusion_matrices skipped: holdout_metrics is None")
        return

    needed = {"y_true", "y_pred"}
    missing = needed - set(holdout_metrics.keys())
    if missing:
        log.warning("plot_confusion_matrices skipped: missing keys %s", sorted(missing))
        return

    y_true = np.asarray(holdout_metrics["y_true"])
    y_pred = np.asarray(holdout_metrics["y_pred"])
    cm = confusion_matrix(y_true, y_pred)

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(5, 4.5))

    if sns is not None:
        sns.heatmap(
            cm,
            annot=True,
            fmt="d",
            cmap="Blues",
            cbar=False,
            xticklabels=["Healthy", "SLE"],
            yticklabels=["Healthy", "SLE"],
            ax=ax,
        )
    else:
        im = ax.imshow(cm, cmap="Blues")
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, str(cm[i, j]), ha="center", va="center", fontsize=12)
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(["Healthy", "SLE"])
        ax.set_yticklabels(["Healthy", "SLE"])
        fig.colorbar(im, ax=ax)

    ax.set_title("Held-out confusion matrix")
    ax.set_xlabel("Predicted label")
    ax.set_ylabel("True label")

    save_figure(fig, fig_dir / "04_confusion_matrix.png", log)

# Permutation Histogram: in order to find if the performance of the model cna be replicated by chance / not because of the distinction of SL / SLE-genes
def plot_permutation_histogram(
    permutation_results: dict,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Plot the null AUC histogram from the permutation test.

    Expected keys:
    - null_aucs
    - observed_auc
    - p_value (optional)
    """
    log = get_log(log)

    if permutation_results is None:
        log.warning("plot_permutation_histogram skipped: permutation_results is None")
        return

    needed = {"null_aucs", "observed_auc"}
    missing = needed - set(permutation_results.keys())
    if missing:
        log.warning("plot_permutation_histogram skipped: missing keys %s", sorted(missing))
        return

    null_aucs = np.asarray(permutation_results["null_aucs"], dtype=float)
    observed_auc = float(permutation_results["observed_auc"])
    p_value = permutation_results.get("p_value", None)

    if len(null_aucs) == 0:
        log.warning("plot_permutation_histogram skipped: null_aucs is empty")
        return

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.hist(null_aucs, bins=20, color=CB_GREY, alpha=0.75, edgecolor="white")
    ax.axvline(
        observed_auc,
        color=CB_BLUE,
        linestyle="--",
        linewidth=2.5,
        label=f"Observed AUC = {observed_auc:.3f}",
    )

    title = "Permutation-test null distribution"
    if p_value is not None:
        title += f" (p = {float(p_value):.4f})"

    ax.set_title(title)
    ax.set_xlabel("Null AUC")
    ax.set_ylabel("Count")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()

    save_figure(fig, fig_dir / "04_permutation_histogram.png", log)

# Probability Violin: in order to explain how predicted probabilites classify SLE vs Healthy
def plot_predicted_probability_violin(
    holdout_metrics: dict,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Plot predicted probabilities grouped by true class.

    Expected keys:
    - y_true
    - y_proba
    """
    log = get_log(log)

    if holdout_metrics is None:
        log.warning("plot_predicted_proba_violin skipped: holdout_metrics is None")
        return

    needed = {"y_true", "y_proba"}
    missing = needed - set(holdout_metrics.keys())
    if missing:
        log.warning("plot_predicted_proba_violin skipped: missing keys %s", sorted(missing))
        return

    df = pd.DataFrame({
        "true_class": np.asarray(holdout_metrics["y_true"]).astype(int),
        "predicted_probability": np.asarray(holdout_metrics["y_proba"]).astype(float),
    })
    df["true_class_label"] = df["true_class"].map({0: "Healthy", 1: "SLE"})

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(6, 5))

    if sns is not None:
        sns.violinplot(
            data=df,
            x="true_class_label",
            y="predicted_probability",
            palette=[CB_ORANGE, CB_BLUE],
            inner="box",
            cut=0,
            ax=ax,
        )
        sns.stripplot(
            data=df,
            x="true_class_label",
            y="predicted_probability",
            color="black",
            alpha=0.55,
            size=4,
            ax=ax,
        )
    else:
        groups = [
            df.loc[df["true_class_label"] == "Healthy", "predicted_probability"].values,
            df.loc[df["true_class_label"] == "SLE", "predicted_probability"].values,
        ]
        parts = ax.violinplot(groups, showmeans=False, showmedians=True, showextrema=True)
        for body, color in zip(parts["bodies"], [CB_ORANGE, CB_BLUE]):
            body.set_facecolor(color)
            body.set_alpha(0.7)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Healthy", "SLE"])

    ax.axhline(0.5, linestyle="--", color=CB_GREY, linewidth=1.5)
    ax.set_ylim(0, 1)
    ax.set_title("Held-out predicted probabilities by true class")
    ax.set_xlabel("True class")
    ax.set_ylabel("Predicted probability for SLE")
    ax.grid(axis="y", alpha=0.3)

    save_figure(fig, fig_dir / "04_predicted_probability_violin.png", log)


# --------------------------------------------------------------------
# SHAP plots
# --------------------------------------------------------------------

# SHAP barplot: in order to show what genes are more important within the signature
def plot_shap_barplot(
    shap_df: pd.DataFrame,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Plot a SHAP summary barplot from a saved SHAP summary dataframe.

    This function tries to detect:
    - gene column
    - importance column
    """
    log = get_log(log)

    if shap_df is None or shap_df.empty:
        log.warning("plot_shap_barplot skipped: shap_df is empty")
        return

    df = shap_df.copy()

    gene_col = None
    for c in ["gene", "feature", "Gene", "Feature"]:
        if c in df.columns:
            gene_col = c
            break

    value_col = None
    for c in ["mean_abs_shap", "mean_abs_value", "importance", "shap_importance", "mean_abs"]:
        if c in df.columns:
            value_col = c
            break

    if gene_col is None:
        gene_col = df.columns[0]

    if value_col is None:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        if numeric_cols:
            value_col = numeric_cols[0]

    if value_col is None:
        log.warning("plot_shap_barplot skipped: no numeric importance column found")
        return

    df = df[[gene_col, value_col]].copy()
    df.columns = ["gene", "importance"]
    df = df.sort_values("importance", ascending=False).head(20)

    fig_dir = ensure_figures_dir(out_dir)
    fig, ax = plt.subplots(figsize=(8, 7))

    ax.barh(df["gene"], df["importance"], color=CB_BLUE, alpha=0.9)
    ax.invert_yaxis()
    ax.set_title("Top SHAP feature importances")
    ax.set_xlabel("Mean absolute SHAP value")
    ax.set_ylabel("Gene")
    ax.grid(axis="x", alpha=0.3)

    save_figure(fig, fig_dir / "04_shap_barplot.png", log)

# SHAP beeswarm: in order to show how/in what degree genes affect the classification between SLE vs Healthy
def plot_shap_beeswarm(
    shap_values,
    Xs_df: pd.DataFrame,
    out_dir: Path,
    log: logging.Logger | None = None,
) -> None:
    """
    Create a SHAP beeswarm plot from per-sample SHAP values.

    Expected:
    - shap_values: array-like of shape (n_samples, n_features)
    - Xs_df: DataFrame with the same feature order
    """
    log = get_log(log)

    if not HAS_SHAP or shap is None:
        log.warning("plot_shap_beeswarm skipped: SHAP is not available")
        return

    if shap_values is None or Xs_df is None or Xs_df.empty:
        log.warning("plot_shap_beeswarm skipped: missing SHAP values or feature matrix")
        return

    shap_array = np.asarray(shap_values)

    if shap_array.ndim != 2:
        log.warning("plot_shap_beeswarm skipped: shap_values must be 2D, got shape %s", shap_array.shape)
        return

    if shap_array.shape[1] != Xs_df.shape[1]:
        log.warning(
            "plot_shap_beeswarm skipped: feature mismatch shap=%s vs X=%s",
            shap_array.shape,
            Xs_df.shape,
        )
        return

    fig_dir = ensure_figures_dir(out_dir)

    plt.figure(figsize=(8, 6))
    shap.summary_plot(
        shap_array,
        Xs_df,
        show=False,
        plot_type="dot",
        max_display=min(20, Xs_df.shape[1]),
    )
    plt.title("SHAP beeswarm plot")
    plt.tight_layout()
    plt.savefig(fig_dir / "04_shap_beeswarm.png", dpi=300, bbox_inches="tight")
    plt.close()

    log.info("saved plot -> %s", fig_dir / "04_shap_beeswarm.png")