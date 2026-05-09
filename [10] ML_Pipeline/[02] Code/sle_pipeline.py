"""
SLE Diagnostic Gene Signature Pipeline

Goal:
Build a donor-aware diagnostic classifier pipeline for SLE from whole-blood
microarray data, starting from a 200-gene candidate set produced by upstream
longitudinal limma analysis.

Main stages:
1. Setup and logging
2. Data loading and sealed donor split
3. Fold-level bootstrap RF stability filtering
4. Feature selection
5. Signature-size selection
6. Outer cross-validation benchmarking
7. Final signature construction and stability
8. Held-out evaluation
"""

# Import Step - installs prerequisites
# ------------------------------------
import requirements as req
import pandas as pd  # in case it won't load later
import numpy as np  # in case it won't load later

# Configuration Step - defines paths, data columns & parameters
# -------------------------------------------------------------
import configurations as config

# Helper Functions Step - develops all the variables, algorithms & deliverables
# -------------------------------------------------------------------------------
import helper_functions as helpers

# Plot Creation Step - develops all plots and diagrams
# -------------------------------------------------------------------------------
import plot_creation as plots

# in case shap doesn't load the previous way
try:
    import shap
    HAS_SHAP = True
except Exception:
    shap = None
    HAS_SHAP = False

# Small compatibility helper in case returned dictionary keys differ slightly
# ---------------------------------------------------------------------------
def _pick(d: dict, *keys, default=None):
    for k in keys:
        if k in d:
            return d[k]
    return default


# ===========================================================================
# main() -- end-to-end thesis run, saving every artefact as it goes
# ===========================================================================

def main() -> None:  # This function returns no actual values; it just formulates and executes the workflow
    """
    Run the full SLE diagnostic pipeline from raw input file to held-out evaluation.

    This function is a workflow controller rather than a value-returning function.
    It coordinates the major steps of the thesis pipeline, writes outputs to disk,
    and records progress in the log file.

    Bootstrap-related methodology checkpoints:
    1. Fold-level bootstrap RF stability filter -> already implemented.
    2. Final-panel bootstrap stability column -> placeholder, revisit later.
    3. Held-out bootstrap AUC evaluation -> placeholder, revisit later.
    """

    # Create the output directory if it does not already exist.
    out_dir = helpers.ensure_output_dir(config.OUTPUT_DIR)  # Ensure [03] Output (and its parents) exist on disk.

    # Start logging to both terminal and log file.
    log = helpers.setup_logging(config.LOG_PATH)  # Configure a log variable that writes to the console and to pipeline_log.txt

    # Fix random seeds for reproducibility.
    helpers.set_global_seed(config.SEED)  # Set Python, NumPy, and hash seeds so repeated runs are reproducible

    log.info("=" * 70)
    log.info("SLE Pipeline -- thesis run")
    log.info("=" * 70)

    log.info(
        "seed=%d holdout_frac=%.2f outer/inner=%d/%d n_repeats=%d",
        config.SEED,
        config.HOLDOUT_FRAC,
        config.N_OUTER_FOLDS,
        config.N_INNER_FOLDS,
        config.N_REPEATS,
    )
    log.info("FS methods: %s", ", ".join(config.FS_METHODS))
    log.info("Classifiers: %s", ", ".join(config.CLASSIFIERS))

    # Bootstrap step 1: this one is definitely active in the current code
    log.info(
        "Fold bootstrap stability: %d iterations, %d trees, top-%d%%",
        config.FOLD_BOOTSTRAP_ITERATIONS,
        config.FOLD_BOOTSTRAP_TREES,
        config.FOLD_BOOTSTRAP_TOP_PERCENTILE,
    )

    # Bootstrap step 3: planned for held-out uncertainty estimation.
    log.info(
        "Held-out resampling plan: %d permutations, %d AUC bootstraps",
        config.N_PERMUTATIONS,
        config.N_HOLDOUT_BOOTSTRAP,
    )

    log.info("Output dir: %s", out_dir)  # Record the absolute path of the output directory

    # -----------------------------------------------------------------------
    # Section 1: load dataset + create a sealed donor-level held-out split
    # -----------------------------------------------------------------------
    if not config.DATA_PATH.exists():  # Check that the configured input file actually exists before proceeding
        raise FileNotFoundError(
            f"Input TSV not found: {config.DATA_PATH}. "
            f"Place the file in the input folder and re-run."
        )

    log.info("-" * 70)
    log.info("Section 1: load + sealed donor split")
    log.info("-" * 70)

    df, X, y, donors, gene_cols = helpers.load_dataset(  # Load the full input table and split it into metadata, gene matrix, labels, donor IDs, and gene names.
        config.DATA_PATH,  # Use the dataset path defined in configurations.py
        logger=log,  # Send progress messages to the pipeline log
    )

    donors_split = helpers.make_sealed_donor_split(
        df,
        X,  # Gene-expression matrix
        y,  # Binary class labels (1 = SLE, 0 = Healthy)
        donors,
        logger=log,
    )

    helpers.save_json(  # Save a small summary of the split for later checking and reporting
        {
            "n_samples_dev": int(len(donors_split["y_dev"])),  # Count how many samples ended up in the development set
            "n_samples_hold": int(len(donors_split["y_hold"])),  # Count how many samples ended up in the held-out set
            "n_donors_dev": int(donors_split["donors_dev"].nunique()),  # Count unique donors in development
            "n_donors_hold": int(donors_split["donors_hold"].nunique()),  # Count unique donors in held-out
            "n_sle_dev": int((donors_split["y_dev"] == 1).sum()),  # Count SLE-labelled samples in development
            "n_sle_hold": int((donors_split["y_hold"] == 1).sum()),  # Count SLE-labelled samples in held-out
            "holdout_donors": sorted(donors_split["hold_donors"]),  # Save the held-out donor IDs in sorted order for auditability
            "n_gene_columns": int(len(gene_cols)),  # Save how many gene columns were loaded
        },
        out_dir / "01_split_summary.json",  # Write the split summary into the output folder
        logger=log,  # Record the save action in the log
    )

    log.info("saved -> 01_split_summary.json")  # Confirm that the split summary file was written successfully

    # ------------------------------------------------------------------------------
    # Section 2: Outer cross-validation on the development set only.
    # ------------------------------------------------------------------------------
    log.info("-" * 70)
    log.info("Section 2: outer cross validation")
    log.info("-" * 70)

    results_df, panels_per_fold = helpers.run_outer_benchmark(
        donors_split["X_dev"],
        donors_split["y_dev"],
        donors_split["donors_dev"],
        logger=log,
    )

    results_df.to_csv(
        out_dir / "02_benchmark_results.csv",
        index=False,
    )

    log.info(
        "saved -> 02_benchmark_results.csv (%d rows)",
        len(results_df),
    )

    summary = helpers.summarise_benchmark(  # Collapse fold-level results into mean, sd, and count per FS/classifier combination.
        results_df,
        logger=log,
    )

    summary.to_csv(
        out_dir / "02_benchmark_summary.csv",
        index=False,
    )

    log.info("saved -> 02_benchmark_summary.csv")

    plots.plot_benchmark_boxplot(results_df, out_dir, log)  # Create a boxplot of benchmark AUC distributions across folds
    plots.plot_benchmark_barplot(summary, out_dir, log)  # Create a barplot of mean benchmark performance
    plots.plot_benchmark_heatmap(summary, out_dir, log)  # Create a heatmap comparing FS/classifier combinations

    # ----------------------------------------------
    # Section 3: final signature panel + bootstrap stability
    # ----------------------------------------------
    log.info("-" * 70)
    log.info("Section 3: final panel + bootstrap stability")
    log.info("-" * 70)

    final = helpers.finalise_panel(
        results_df,
        donors_split["X_dev"],
        donors_split["y_dev"],
        donors_split["donors_dev"],
        logger=log,
    )

    # Save full Pearson correlation matrix for the ranked stable genes
    corr_full = donors_split["X_dev"][final["high_freq_ranked_genes"]].corr(method="pearson")
    corr_full.to_csv(
        config.OUTPUT_DIR / "03_final_panel_pearson_full_matrix.csv",
        index=True,
    )
    log.info("saved -> 03_final_panel_pearson_full_matrix.csv")

    # Save only pairs above the threshold for quick inspection
    high_corr_rows = []
    for i, g1 in enumerate(final["high_freq_ranked_genes"]):
        for g2 in final["high_freq_ranked_genes"][i + 1:]:
            r = float(corr_full.loc[g1, g2])
            if abs(r) >= config.FINAL_PEARSON_THRESHOLD:
                high_corr_rows.append(
                    {
                        "gene_1": g1,
                        "gene_2": g2,
                        "pearson_r": r,
                        "abs_r": abs(r),
                    }
                )

    if high_corr_rows:
        pd.DataFrame(high_corr_rows).to_csv(
            config.OUTPUT_DIR / "03_final_panel_pearson_high_pairs.csv",
            index=False,
        )
        log.info(
            "saved -> 03_final_panel_pearson_high_pairs.csv (%d pairs |r| >= %.2f)",
            len(high_corr_rows),
            config.FINAL_PEARSON_THRESHOLD,
        )

    final["panel_table"].to_csv(
        out_dir / "03_final_panel.csv",
        index=False,
    )

    log.info(
        "saved -> 03_final_panel.csv (%d genes)",
        len(final["panel_table"]),
    )

    stability_series = _pick(
        final,
        "stability",
        "stability_freq",
        default=pd.Series(dtype=float),
    )
    if isinstance(stability_series, pd.Series) and not stability_series.empty:
        stability_series.rename("bootstrap_freq").to_csv(
            out_dir / "03_full_stability.csv",
            header=True,
        )

        log.info(
            "saved -> 03_full_stability.csv (%d unique genes)",
            len(stability_series),
        )

    pd.DataFrame(
        {
            "k": np.arange(2, 2 + len(final["auc_curve"])),
            "mean_inner_cv_auc": final["auc_curve"],
        }
    ).to_csv(
        out_dir / "03_elbow_curve.csv",
        index=False,
    )

    log.info("saved -> 03_elbow_curve.csv")

    helpers.save_json(
        {
            "winner": final["winner"],
            "ranking": final["ranking"],
            "panel": final["panel"],
            "high_freq_ranked_genes": final.get("high_freq_ranked_genes", []),
            "elbow_k": final["elbow_k"],
            "final_stability_threshold": config.FINAL_STABILITY_THRESHOLD,
            "final_pearson_threshold": config.FINAL_PEARSON_THRESHOLD,
        },
        out_dir / "03_final_panel_meta.json",
        logger=log,
    )

    # After you have `final` (with keys: ranking, high_freq_ranked_genes, panel, stability)

    stability_series = final["stability"]  # pd.Series: index = gene, value = bootstrap freq

    high_freq_genes = list(final["high_freq_ranked_genes"])

    high_freq_df = pd.DataFrame({
        "gene": high_freq_genes,
        # FS-only rank: position of gene in the full Boruta ranking
        "fs_rank": [final["ranking"].index(g) + 1 for g in high_freq_genes],
    })

    # Add the stability frequency (from the full stability bootstrap)
    high_freq_df["bootstrap_selection_freq"] = (
        high_freq_df["gene"].map(stability_series).fillna(0.0)
    )

    # Mark which of these 29 genes are in the final 4-gene panel
    high_freq_df["in_final_panel"] = high_freq_df["gene"].isin(final["panel"])

    # Save only this CSV
    high_freq_df.to_csv(
        config.OUTPUT_DIR / "03_high_freq_ranked_genes.csv",
        index=False,
    )
    log.info("saved -> 03_high_freq_ranked_genes.csv")   

    log.info("saved -> 03_final_panel_meta.json")

    corr_df = donors_split["X_dev"][final["panel"]].corr(method="pearson")
    corr_df.to_csv(out_dir / "03_final_panel_correlation_matrix.csv", index=True)
    log.info("saved -> 03_final_panel_correlation_matrix.csv")

    plots.plot_elbow_curve(
        final["auc_curve"],
        final["elbow_k"],
        out_dir,
        log,
    )

    plots.plot_stability_barplot(
        final["panel_table"],
        out_dir,
        log,
    )

    plots.plot_final_panel_correlation_heatmap(
        donors_split["X_dev"],
        final["panel"],
        out_dir,
        log,
    )

    pd.DataFrame({
    "gene": final["high_freq_ranked_genes"]
        }).to_csv(
            out_dir / "03_ranked_stable_genes_pre_pearson.csv",
            index=False,
        )
    log.info("saved -> 03_ranked_stable_genes_pre_pearson.csv")

    if isinstance(stability_series, pd.Series) and not stability_series.empty:
        plots.plot_full_stability_top20(
            stability_series,
            out_dir,
            log,
        )


    # ---------------------------------------------
    # Section 4: "TEST" held-out donor subset evaluation
    # ---------------------------------------------
    log.info("-" * 70)
    log.info("Section 4: held-out evaluation")
    log.info("-" * 70)


    from sklearn.metrics import roc_curve


    # Train the final model on the full development set.
    fitted_model = helpers.train_final_model(
        X_dev=donors_split["X_dev"],
        y_dev=donors_split["y_dev"],
        donors_dev=donors_split["donors_dev"],
        final_signature=final["panel"],
        classifier_name=final["winner"]["classifier"],
        logger=log,
        seed=config.SEED,
    )


    # Evaluate the trained final model on the sealed held-out set.
    holdout_metrics = helpers.evaluate_on_holdout(
        fitted=fitted_model,
        X_hold=donors_split["X_hold"],
        y_hold=donors_split["y_hold"],
        logger=log,
    )


    # Add ROC curve coordinates for downstream export and plotting.
    fpr, tpr, _ = roc_curve(
        np.asarray(holdout_metrics["y_true"]),
        np.asarray(holdout_metrics["y_proba"]),
    )
    holdout_metrics["fpr"] = fpr.tolist()
    holdout_metrics["tpr"] = tpr.tolist()


    # Held-out AUC bootstrap confidence interval.
    ci_info = helpers.bootstrap_auc_ci(
        ytrue=np.asarray(holdout_metrics["y_true"]),
        yproba=np.asarray(holdout_metrics["y_proba"]),
        nbootstrap=config.N_HOLDOUT_BOOTSTRAP,
        seed=config.SEED,
        logger=log,
    )


    # Copy AUC confidence interval values into holdout_metrics so everything
    # lives in one dictionary for logging, JSON export, and plotting.
    holdout_metrics["auc_ci_lower"] = ci_info["ci_lo"]
    holdout_metrics["auc_ci_upper"] = ci_info["ci_hi"]
    holdout_metrics["n_auc_bootstrap"] = ci_info["n_valid_bootstrap"]

    # Log a user-friendly held-out summary block.
    log.info("PART A: Held-Out Summary")
    log.info("  Metric                 Value")
    log.info("  ----------------------------------------------")
    log.info(
        "  AUC (95%% CI)           %.4f [%.4f-%.4f]",
        holdout_metrics["auc"],
        holdout_metrics["auc_ci_lower"],
        holdout_metrics["auc_ci_upper"],
    )
    log.info("  Sensitivity            %.4f", holdout_metrics["sensitivity"])
    log.info("  Specificity            %.4f", holdout_metrics["specificity"])
    log.info("  Accuracy               %.4f", holdout_metrics["accuracy"])
    log.info("  F1 (held-out ratio)    %.4f", holdout_metrics["f1"])
    log.info(
        "  TP / FP / FN / TN       %d/%d/%d/%d",
        holdout_metrics["tp"],
        holdout_metrics["fp"],
        holdout_metrics["fn"],
        holdout_metrics["tn"],
    )

    n_pos = holdout_metrics["tp"] + holdout_metrics["fn"]
    n_neg = holdout_metrics["tn"] + holdout_metrics["fp"]
    n_total = holdout_metrics["n_samples"]
    prev_pos = n_pos / n_total if n_total > 0 else 0.0

    log.info(
        "  Prevalence (held-out)  SLE %.1f%% (%d/%d), Healthy %.1f%% (%d/%d)",
        100 * prev_pos,
        n_pos,
        n_total,
        100 * (1 - prev_pos),
        n_neg,
        n_total,
    )

    log.info(
        "  NOTE: AUC 95%% CI computed by %d-iteration bootstrap on held-out predictions.",
        holdout_metrics["n_auc_bootstrap"],
    )


    # Donor-stratified permutation test on the development set.
    permutation_results = helpers.permutation(
        Xdev=donors_split["X_dev"],
        ydev=donors_split["y_dev"],
        donorsdev=donors_split["donors_dev"],
        finalsignature=fitted_model["signature"],
        classifiername=fitted_model["classifier"],
        bestparams=fitted_model["best_params"],
        logger=log,
        npermutations=config.N_PERMUTATIONS,
        nfolds=config.N_INNER_FOLDS,
        seed=config.SEED,
    )


    # Logistic-regression comparator on the same final panel (optional baseline).
    lr_comparison = None

    # Only run LR comparator if the winning classifier is not already logreg.
    if final["winner"]["classifier"] != "logreg":
        log.info(
            "PART B: LR comparator | Winner is %s — training logistic regression baseline on the same panel.",
            final["winner"]["classifier"],
        )
        try:
            lr_fitted = helpers.train_final_model(
                X_dev=donors_split["X_dev"],
                y_dev=donors_split["y_dev"],
                donors_dev=donors_split["donors_dev"],
                final_signature=final["panel"],
                classifier_name="logreg",
                logger=log,
                seed=config.SEED,
            )

            lr_holdout = helpers.evaluate_on_holdout(
                fitted=lr_fitted,
                X_hold=donors_split["X_hold"],
                y_hold=donors_split["y_hold"],
                logger=log,
            )

            lr_fpr, lr_tpr, _ = roc_curve(
                np.asarray(lr_holdout["y_true"]),
                np.asarray(lr_holdout["y_proba"]),
            )
            lr_holdout["fpr"] = lr_fpr.tolist()
            lr_holdout["tpr"] = lr_tpr.tolist()

            lr_ci = helpers.bootstrap_auc_ci(
                ytrue=np.asarray(lr_holdout["y_true"]),
                yproba=np.asarray(lr_holdout["y_proba"]),
                nbootstrap=config.N_HOLDOUT_BOOTSTRAP,
                seed=config.SEED,
                logger=log,
            )

            # Store everything needed for plotting and saving
            lr_comparison = {
                "fitted": lr_fitted,
                "holdout": lr_holdout,
                "auc": lr_holdout["auc"],
                "ci_lo": lr_ci["ci_lo"],
                "ci_hi": lr_ci["ci_hi"],
            }

            log.info(
                "PART B: LR comparator | Held-out AUC = %.4f [%.4f–%.4f] | F1 = %.3f",
                lr_holdout["auc"],
                lr_ci["ci_lo"],
                lr_ci["ci_hi"],
                lr_holdout["accuracy"],
            )

        except Exception as exc:
            log.warning("LR comparator skipped: %s", exc)
    else:
        log.info(
            "PART B: LR comparator skipped because winning classifier is already logreg."
        )
        lr_comparison = None


    # SHAP summary and optional held-out SHAP values.
    shap_df = pd.DataFrame()


    if "HAS_SHAP" in globals() and HAS_SHAP:
        try:
            Xp = donors_split["X_hold"][fitted_model["signature"]]
            Xs = fitted_model["scaler"].transform(Xp)
            Xs_df = pd.DataFrame(Xs, columns=fitted_model["signature"])


            model_name = str(fitted_model["classifier"]).lower()


            if model_name in {"rf", "xgboost"}:
                explainer = shap.TreeExplainer(fitted_model["model"])
                shap_values_raw = explainer.shap_values(Xs_df)
            else:
                bg = shap.kmeans(Xs_df, min(10, len(Xs_df)))
                explainer = shap.KernelExplainer(
                    lambda x: fitted_model["model"].predict_proba(x)[:, 1],
                    bg,
                )
                shap_values_raw = explainer.shap_values(
                    Xs_df,
                    nsamples=100,
                    silent=True,
                )


            if isinstance(shap_values_raw, list):
                if len(shap_values_raw) > 1:
                    final_shap_values = np.asarray(shap_values_raw[1])
                else:
                    final_shap_values = np.asarray(shap_values_raw[0])
            else:
                final_shap_values = np.asarray(shap_values_raw)


            if final_shap_values.ndim == 3:
                final_shap_values = final_shap_values[:, :, 1]


            shap_df = pd.DataFrame(
                {
                    "gene": fitted_model["signature"],
                    "mean_abs_shap": np.mean(np.abs(final_shap_values), axis=0),
                }
            ).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)


            pd.DataFrame(
                final_shap_values,
                columns=fitted_model["signature"],
            ).to_csv(
                out_dir / "04_shap_values.csv",
                index=False,
            )
            log.info("saved -> 04_shap_values.csv")


        except Exception as exc:
            log.warning("SHAP summary/values skipped: %s", exc)


    # Predictions table (one row per held-out sample).
    prediction_df = pd.DataFrame({
        "y_true": holdout_metrics["y_true"],
        "y_prediction": holdout_metrics["y_pred"],
        "y_probability": holdout_metrics["y_proba"],
    })
    prediction_df.to_csv(
        out_dir / "04_predictions.csv",
        index=False,
    )
    log.info("saved -> 04_predictions.csv")


    # ROC curve points.
    pd.DataFrame({
        "fpr": holdout_metrics["fpr"],
        "tpr": holdout_metrics["tpr"],
    }).to_csv(
        out_dir / "04_roc_data.csv",
        index=False,
    )
    log.info("saved -> 04_roc_data.csv")


    n_pos = holdout_metrics["tp"] + holdout_metrics["fn"]
    n_neg = holdout_metrics["tn"] + holdout_metrics["fp"]
    n_total = holdout_metrics["n_samples"]

    helpers.save_json(
        {
            "auc": holdout_metrics["auc"],
            "auc_ci_lower": holdout_metrics["auc_ci_lower"],
            "auc_ci_upper": holdout_metrics["auc_ci_upper"],
            "n_auc_bootstrap": holdout_metrics["n_auc_bootstrap"],
            "sensitivity": holdout_metrics["sensitivity"],
            "specificity": holdout_metrics["specificity"],
            "accuracy": holdout_metrics["accuracy"],
            "f1": holdout_metrics["f1"],
            "tp": holdout_metrics["tp"],
            "fp": holdout_metrics["fp"],
            "fn": holdout_metrics["fn"],
            "tn": holdout_metrics["tn"],
            "n_positive": int(n_pos),
            "n_negative": int(n_neg),
            "n_samples": int(n_total),
            "best_params": fitted_model["best_params"],
            "inner_cv_auc": fitted_model["inner_cv_auc"],
            "panel": fitted_model["signature"],
            "classifier": final["winner"]["classifier"],
            "fs_method": final["winner"]["fs_method"],
            "final_stability_threshold": config.FINAL_STABILITY_THRESHOLD,
            "final_pearson_threshold": config.FINAL_PEARSON_THRESHOLD,
        },
        out_dir / "04_holdout_metrics.json",
        logger=log,
    )


    pd.DataFrame({
        "null_auc": permutation_results["null_aucs"],
    }).to_csv(
        out_dir / "04_permutation_null.csv",
        index=False,
    )


    helpers.save_json(
        {
            "observed_auc": permutation_results["observed_auc"],
            "p_value": permutation_results["p_value"],
            "n_permutations": permutation_results["n_permutations"],
            "null_mean": (
                float(np.mean(permutation_results["null_aucs"]))
                if len(permutation_results["null_aucs"]) else None
            ),
        },
        out_dir / "04_permutation_meta.json",
        logger=log,
    )
    log.info("saved -> 04_permutation_null.csv + 04_permutation_meta.json")


    if not shap_df.empty:
        shap_df.to_csv(
            out_dir / "04_shap_summary.csv",
            index=False,
        )
        log.info("saved -> 04_shap_summary.csv")


    # LR comparator overlay ROC if available.
    if lr_comparison is not None:
        lr_holdout = lr_comparison["holdout"]

        helpers.save_json(
            {
                "auc": lr_comparison["auc"],
                "ci_lo": lr_comparison["ci_lo"],
                "ci_hi": lr_comparison["ci_hi"],
                "best_params": lr_comparison["fitted"]["best_params"],
            },
            out_dir / "04_lr_comparator.json",
            logger=log,
        )
        log.info("saved -> 04_lr_comparator.json")

        plots.plot_roc_curve_with_lr(
            holdout_metrics, lr_holdout, ci_info, lr_comparison, out_dir, log
        )
    else:
        plots.plot_roc_curve(
            holdout_metrics, None, ci_info, out_dir, log
        )


    plots.plot_confusion_matrices(holdout_metrics, out_dir, log)
    plots.plot_permutation_histogram(permutation_results, out_dir, log)
    plots.plot_predicted_probability_violin(holdout_metrics, out_dir, log)


    if not shap_df.empty and "HAS_SHAP" in globals() and HAS_SHAP:
        plots.plot_shap_barplot(shap_df, out_dir, log)


        try:
            Xp = donors_split["X_hold"][fitted_model["signature"]]
            Xs = fitted_model["scaler"].transform(Xp)
            Xs_df = pd.DataFrame(Xs, columns=fitted_model["signature"])


            if "04_shap_values.csv":
                final_shap_values = pd.read_csv(out_dir / "04_shap_values.csv").values
                plots.plot_shap_beeswarm(final_shap_values, Xs_df, out_dir, log)


        except Exception as exc:
            log.warning("SHAP beeswarm skipped: %s", exc)


if __name__ == "__main__":
    main()