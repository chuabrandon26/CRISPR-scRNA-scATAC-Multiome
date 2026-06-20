"""CRISPR like single cell multiome analysis pipeline.

This script loads the 10x Genomics mouse brain multiome dataset, preprocesses
RNA and ATAC data, creates virtual perturbation labels, and runs a conceptual
CRISPR style analysis for gene expression and chromatin accessibility.
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

from scipy import sparse, stats
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans


DATA_DIR = Path("/mnt/c/Users/chuab/Desktop/CRISPR_scRNA_scATAC")
FIGURE_DIR = DATA_DIR / "images"
FIGURE_DIR.mkdir(exist_ok=True)

H5_PATH = DATA_DIR / "Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_filtered_feature_bc_matrix.h5"
PEAK_BED_PATH = DATA_DIR / "Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_atac_peaks.bed"
RANDOM_SEED = 7

NEUROINFLAMMATION_GENES = [
    "Apoe",
    "Trem2",
    "Tyrobp",
    "C1qa",
    "C1qb",
    "C1qc",
    "Lpl",
    "Cst7",
    "Itgax",
    "B2m",
]

SYNAPTIC_GENES = [
    "Snap25",
    "Syp",
    "Rbfox3",
    "Map2",
    "Grin1",
    "Camk2a",
    "Slc17a7",
    "Gad1",
]


def save_figure(filename: str) -> Path:
    """Save the active matplotlib figure in the project image folder."""
    path = FIGURE_DIR / filename
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()
    print(f"Saved {path}")
    return path


def to_1d(values) -> np.ndarray:
    """Return a one dimensional NumPy array from sparse or dense input."""
    if sparse.issparse(values):
        return np.asarray(values.toarray()).ravel()
    return np.asarray(values).ravel()


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    """Benjamini Hochberg adjusted p values for many tests."""
    p_values = np.asarray(p_values, dtype=float)
    adjusted = np.full_like(p_values, np.nan, dtype=float)
    finite = np.isfinite(p_values)
    if not finite.any():
        return adjusted

    p = p_values[finite]
    order = np.argsort(p)
    ranked = p[order]
    n = len(ranked)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    finite_indices = np.where(finite)[0]
    adjusted[finite_indices[order]] = q
    return adjusted


def load_multiome_h5(h5_path: Path = H5_PATH) -> ad.AnnData:
    """Load the 10x multiome H5 file and keep both RNA and ATAC features."""
    if not h5_path.exists():
        raise FileNotFoundError(f"Could not find {h5_path}")

    # gex_only=False is important because the file contains genes and peaks.
    try:
        adata = sc.read_10x_h5(h5_path, gex_only=False)
    except TypeError:
        adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    print(adata)
    print(adata.var["feature_types"].value_counts())
    return adata


def split_modalities(adata: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    """Split one multiome AnnData object into RNA and ATAC AnnData objects."""
    if "feature_types" not in adata.var:
        raise KeyError("The 10x object does not contain adata.var['feature_types'].")

    feature_types = adata.var["feature_types"].astype(str)
    rna_mask = feature_types.str.contains("Gene Expression", case=False, na=False)
    atac_mask = feature_types.str.contains("Peak|Chromatin|Accessibility", case=False, na=False)

    if not atac_mask.any():
        atac_mask = ~rna_mask

    adata_rna = adata[:, rna_mask].copy()
    adata_atac = adata[:, atac_mask].copy()
    adata_rna.var_names_make_unique()
    adata_atac.var_names_make_unique()
    print(f"RNA shape: {adata_rna.shape}")
    print(f"ATAC shape: {adata_atac.shape}")
    return adata_rna, adata_atac


def add_rna_qc(adata_rna: ad.AnnData) -> None:
    """Add common RNA quality control metrics to adata_rna.obs."""
    var_upper = adata_rna.var_names.str.upper()
    adata_rna.var["mt"] = var_upper.str.startswith("MT-") | var_upper.str.startswith("MT.")
    if adata_rna.var["mt"].any():
        sc.pp.calculate_qc_metrics(
            adata_rna,
            qc_vars=["mt"],
            inplace=True,
            percent_top=None,
            log1p=False,
        )
    else:
        sc.pp.calculate_qc_metrics(
            adata_rna,
            inplace=True,
            percent_top=None,
            log1p=False,
        )
        adata_rna.obs["pct_counts_mt"] = 0.0


def filter_rna(adata_rna: ad.AnnData) -> ad.AnnData:
    """Remove low quality RNA cells and genes with simple, transparent rules."""
    add_rna_qc(adata_rna)
    sc.pp.filter_cells(adata_rna, min_genes=200)
    sc.pp.filter_cells(adata_rna, min_counts=500)
    if "pct_counts_mt" in adata_rna.obs:
        adata_rna = adata_rna[adata_rna.obs["pct_counts_mt"] < 20].copy()
    sc.pp.filter_genes(adata_rna, min_cells=3)
    print(f"RNA shape after QC: {adata_rna.shape}")
    return adata_rna


def add_atac_qc(adata_atac: ad.AnnData) -> None:
    """Add ATAC quality control metrics based on the peak count matrix."""
    if sparse.issparse(adata_atac.X):
        x_csr = adata_atac.X.tocsr()
        fragments = np.asarray(x_csr.sum(axis=1)).ravel()
        peaks = x_csr.getnnz(axis=1)
    else:
        fragments = np.asarray(adata_atac.X.sum(axis=1)).ravel()
        peaks = np.asarray((adata_atac.X > 0).sum(axis=1)).ravel()

    adata_atac.obs["fragments_per_cell"] = fragments
    adata_atac.obs["peaks_per_cell"] = peaks


def filter_atac(adata_atac: ad.AnnData) -> ad.AnnData:
    """Remove low quality ATAC cells and rarely detected peaks."""
    add_atac_qc(adata_atac)
    fragments = adata_atac.obs["fragments_per_cell"].to_numpy()
    peaks = adata_atac.obs["peaks_per_cell"].to_numpy()

    min_fragments = max(500.0, float(np.quantile(fragments, 0.02)))
    min_peaks = max(250.0, float(np.quantile(peaks, 0.02)))
    max_fragments = float(np.quantile(fragments, 0.995))

    cell_mask = (
        (fragments >= min_fragments)
        & (peaks >= min_peaks)
        & (fragments <= max_fragments)
    )
    adata_atac = adata_atac[cell_mask].copy()

    if sparse.issparse(adata_atac.X):
        peak_cells = adata_atac.X.tocsr().getnnz(axis=0)
    else:
        peak_cells = np.asarray((adata_atac.X > 0).sum(axis=0)).ravel()

    adata_atac = adata_atac[:, peak_cells >= 10].copy()
    add_atac_qc(adata_atac)
    print(f"ATAC shape after QC: {adata_atac.shape}")
    return adata_atac


def align_cells(adata_rna: ad.AnnData, adata_atac: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    """Keep only cells that pass QC in both RNA and ATAC."""
    common_cells = adata_rna.obs_names.intersection(adata_atac.obs_names)
    if len(common_cells) == 0:
        raise ValueError("No shared cells remain after QC.")

    adata_rna = adata_rna[common_cells].copy()
    adata_atac = adata_atac[common_cells].copy()
    print(f"Shared cells after QC: {len(common_cells)}")
    return adata_rna, adata_atac


def run_leiden_or_fallback(adata_obj: ad.AnnData, key_added: str, resolution: float) -> None:
    """Run Leiden clustering, with a small fallback if Leiden is unavailable."""
    try:
        sc.tl.leiden(
            adata_obj,
            key_added=key_added,
            resolution=resolution,
            random_state=RANDOM_SEED,
        )
    except Exception as exc:
        warnings.warn(
            f"Leiden clustering failed with {exc}. "
            f"Using KMeans labels in {key_added} so the rest of the workflow can continue."
        )
        rep = "X_pca" if "X_pca" in adata_obj.obsm else "X_lsi"
        matrix = adata_obj.obsm[rep]
        n_clusters = min(10, max(2, adata_obj.n_obs // 250))
        labels = KMeans(n_clusters=n_clusters, random_state=RANDOM_SEED, n_init=10).fit_predict(matrix)
        adata_obj.obs[key_added] = pd.Categorical(labels.astype(str))


def preprocess_rna(adata_rna: ad.AnnData) -> ad.AnnData:
    """Normalize, log transform, run PCA, compute neighbors, UMAP, and clusters."""
    adata_rna.layers["counts"] = adata_rna.X.copy()
    sc.pp.normalize_total(adata_rna, target_sum=1e4)
    sc.pp.log1p(adata_rna)

    n_top_genes = min(3000, adata_rna.n_vars)
    sc.pp.highly_variable_genes(adata_rna, n_top_genes=n_top_genes, flavor="seurat")
    sc.tl.pca(
        adata_rna,
        n_comps=min(50, adata_rna.n_obs - 1, adata_rna.n_vars - 1),
        use_highly_variable=True,
        svd_solver="arpack",
        random_state=RANDOM_SEED,
    )
    sc.pp.neighbors(adata_rna, n_neighbors=15, n_pcs=min(30, adata_rna.obsm["X_pca"].shape[1]))
    sc.tl.umap(adata_rna, random_state=RANDOM_SEED)
    run_leiden_or_fallback(adata_rna, "leiden_rna", resolution=0.8)
    return adata_rna


def run_tfidf_lsi(adata_atac: ad.AnnData, n_components: int = 50) -> ad.AnnData:
    """Run TF-IDF and truncated SVD to make an LSI embedding for scATAC."""
    x = adata_atac.X.tocsr() if sparse.issparse(adata_atac.X) else sparse.csr_matrix(adata_atac.X)
    cell_sums = np.asarray(x.sum(axis=1)).ravel()
    cell_sums[cell_sums == 0] = 1

    tf = sparse.diags(1.0 / cell_sums).dot(x)
    peak_cells = np.asarray((x > 0).sum(axis=0)).ravel()
    idf = np.log1p(x.shape[0] / (1.0 + peak_cells))
    x_tfidf = tf.multiply(idf)

    max_components = min(n_components, x_tfidf.shape[0] - 1, x_tfidf.shape[1] - 1)
    if max_components < 2:
        raise ValueError("Not enough cells or peaks remain to compute LSI.")

    svd = TruncatedSVD(n_components=max_components, random_state=RANDOM_SEED)
    lsi = svd.fit_transform(x_tfidf)

    from scipy.sparse import issparse, csr_matrix
    adata_atac.layers["tfidf"] = csr_matrix(x_tfidf)
    adata_atac.obsm["X_lsi"] = lsi
    adata_atac.uns["lsi_variance_ratio"] = svd.explained_variance_ratio_
    return adata_atac


def preprocess_atac(adata_atac: ad.AnnData) -> ad.AnnData:
    """Run LSI, compute neighbors, UMAP, and clusters for ATAC data."""
    run_tfidf_lsi(adata_atac)
    use_dims = min(30, adata_atac.obsm["X_lsi"].shape[1])
    sc.pp.neighbors(adata_atac, n_neighbors=15, use_rep="X_lsi", n_pcs=use_dims)
    sc.tl.umap(adata_atac, random_state=RANDOM_SEED)
    run_leiden_or_fallback(adata_atac, "leiden_atac", resolution=0.8)
    return adata_atac


def plot_qc(adata_rna: ad.AnnData, adata_atac: ad.AnnData) -> None:
    """Save compact RNA and ATAC QC figures."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    sns.histplot(adata_rna.obs["total_counts"], bins=60, ax=axes[0])
    axes[0].set_title("RNA counts per cell")
    sns.histplot(adata_rna.obs["n_genes_by_counts"], bins=60, ax=axes[1])
    axes[1].set_title("RNA genes per cell")
    sns.histplot(adata_rna.obs["pct_counts_mt"], bins=60, ax=axes[2])
    axes[2].set_title("Mitochondrial percent")
    plt.tight_layout()
    save_figure("rna_qc_metrics.png")

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    sns.histplot(adata_atac.obs["fragments_per_cell"], bins=60, ax=axes[0])
    axes[0].set_title("ATAC fragments per cell")
    sns.histplot(adata_atac.obs["peaks_per_cell"], bins=60, ax=axes[1])
    axes[1].set_title("ATAC peaks per cell")
    plt.tight_layout()
    save_figure("atac_qc_metrics.png")


def plot_umap(adata_obj: ad.AnnData, color: str, filename: str, title: str) -> None:
    """Save a Scanpy UMAP plot."""
    sc.pl.umap(adata_obj, color=color, show=False, frameon=False, title=title)
    save_figure(filename)


def add_gene_score(adata_rna: ad.AnnData, genes: list[str], score_name: str) -> list[str]:
    """Score cells by the average signal of a small gene set."""
    present = [gene for gene in genes if gene in adata_rna.var_names]
    if len(present) >= 2:
        sc.tl.score_genes(adata_rna, gene_list=present, score_name=score_name, random_state=RANDOM_SEED)
    elif len(present) == 1:
        adata_rna.obs[score_name] = to_1d(adata_rna[:, present[0]].X)
    else:
        adata_rna.obs[score_name] = adata_rna.obsm["X_pca"][:, 0]
        warnings.warn(
            f"No genes from {score_name} were found. "
            "Using the first RNA principal component as a placeholder score."
        )
    return present


def add_virtual_perturbations(adata_rna: ad.AnnData, adata_atac: ad.AnnData) -> tuple[ad.AnnData, ad.AnnData]:
    """Create conceptual CRISPR guide labels from an RNA disease pathway score."""
    neuro_genes = add_gene_score(adata_rna, NEUROINFLAMMATION_GENES, "neuroinflammation_score")
    synaptic_genes = add_gene_score(adata_rna, SYNAPTIC_GENES, "synaptic_score")

    cutoff = float(adata_rna.obs["neuroinflammation_score"].median())
    labels = np.where(
        adata_rna.obs["neuroinflammation_score"] >= cutoff,
        "Perturbation A",
        "Perturbation B",
    )
    categories = ["Perturbation A", "Perturbation B"]
    adata_rna.obs["virtual_perturbation"] = pd.Categorical(labels, categories=categories)
    adata_atac.obs["virtual_perturbation"] = adata_rna.obs.loc[
        adata_atac.obs_names,
        "virtual_perturbation",
    ].astype("category")

    adata_rna.uns["virtual_perturbation_note"] = {
        "rule": "Perturbation A is high neuroinflammation score. Perturbation B is low neuroinflammation score.",
        "neuroinflammation_genes_found": neuro_genes,
        "synaptic_genes_found": synaptic_genes,
    }
    return adata_rna, adata_atac


def choose_cluster_with_groups(
    adata_obj: ad.AnnData,
    cluster_key: str,
    group_key: str = "virtual_perturbation",
    min_cells_per_group: int = 20,
) -> str:
    """Choose a cluster that has enough cells from both perturbation groups."""
    cluster_sizes = adata_obj.obs[cluster_key].value_counts()
    for cluster in cluster_sizes.index:
        counts = adata_obj.obs.loc[adata_obj.obs[cluster_key] == cluster, group_key].value_counts()
        if all(counts.get(group, 0) >= min_cells_per_group for group in ["Perturbation A", "Perturbation B"]):
            return str(cluster)

    return str(cluster_sizes.index[0])


def run_rna_differential_expression(adata_rna: ad.AnnData) -> pd.DataFrame:
    """Compare Perturbation A and B inside one RNA cluster."""
    cluster = choose_cluster_with_groups(adata_rna, "leiden_rna")
    adata_cluster = adata_rna[adata_rna.obs["leiden_rna"] == cluster].copy()
    print(f"RNA differential expression cluster: {cluster}")

    sc.tl.rank_genes_groups(
        adata_cluster,
        groupby="virtual_perturbation",
        groups=["Perturbation A"],
        reference="Perturbation B",
        method="wilcoxon",
    )
    de = sc.get.rank_genes_groups_df(adata_cluster, group="Perturbation A")
    de = de.rename(columns={"names": "gene"})
    de["cluster"] = cluster
    de_path = DATA_DIR / "rna_de_perturbationA_vs_B.csv"
    de.to_csv(de_path, index=False)
    print(f"Saved {de_path}")
    return de


def plot_rna_volcano(de: pd.DataFrame) -> None:
    """Save a volcano plot for the RNA differential expression results."""
    plot_df = de.replace([np.inf, -np.inf], np.nan).dropna(subset=["logfoldchanges", "pvals_adj"]).copy()
    plot_df["minus_log10_padj"] = -np.log10(np.clip(plot_df["pvals_adj"], 1e-300, 1))

    plt.figure(figsize=(6, 5))
    sns.scatterplot(
        data=plot_df,
        x="logfoldchanges",
        y="minus_log10_padj",
        s=12,
        alpha=0.5,
        edgecolor=None,
    )
    plt.axvline(0, color="black", linewidth=0.8)
    plt.xlabel("log fold change, Perturbation A versus B")
    plt.ylabel("-log10 adjusted p value")
    plt.title("RNA differential expression")

    top = plot_df.sort_values(["pvals_adj", "logfoldchanges"], ascending=[True, False]).head(10)
    for _, row in top.iterrows():
        plt.text(row["logfoldchanges"], row["minus_log10_padj"], str(row["gene"]), fontsize=7)

    save_figure("rna_de_perturbationA_vs_B_volcano.png")


def plot_rna_pathway_scores(adata_rna: ad.AnnData) -> None:
    """Save pathway score summaries for virtual perturbation groups."""
    score_df = adata_rna.obs[
        ["virtual_perturbation", "neuroinflammation_score", "synaptic_score"]
    ].copy()
    score_long = score_df.melt(
        id_vars="virtual_perturbation",
        var_name="score",
        value_name="value",
    )

    plt.figure(figsize=(7, 4))
    ax = sns.boxplot(
        data=score_long,
        x="score",
        y="value",
        hue="virtual_perturbation",
        fliersize=0,
    )
    sns.stripplot(
        data=score_long.sample(min(2000, len(score_long)), random_state=RANDOM_SEED),
        x="score",
        y="value",
        hue="virtual_perturbation",
        dodge=True,
        color="black",
        alpha=0.15,
        size=1,
    )
    plt.xlabel("")
    plt.ylabel("score")
    plt.title("RNA pathway scores by virtual perturbation")
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[:2], labels[:2], title="virtual perturbation", bbox_to_anchor=(1.02, 1), loc="upper left")
    save_figure("rna_pathway_scores_virtual_perturbations.png")


def load_peak_bed(bed_path: Path = PEAK_BED_PATH) -> pd.DataFrame:
    """Load the 10x peak BED file and create peak IDs that match Scanpy names."""
    peak_df = pd.read_csv(
        bed_path,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
    )
    peak_df["peak"] = (
        peak_df["chrom"].astype(str)
        + ":"
        + peak_df["start"].astype(str)
        + "-"
        + peak_df["end"].astype(str)
    )
    peak_df["peak_width"] = peak_df["end"] - peak_df["start"]
    return peak_df


def annotate_da_peaks_with_coordinates(da: pd.DataFrame) -> pd.DataFrame:
    """Attach BED coordinates to differential accessibility results."""
    peaks = load_peak_bed()
    annotated = da.merge(peaks, how="left", on="peak")

    annotated["nearest_gene"] = "not annotated"
    annotated["motif_annotation"] = "not run"
    annotated["annotation_note"] = (
        "Coordinates come from the 10x BED file. Nearest genes need a mouse TSS table. "
        "Motifs need a genome FASTA and motif database."
    )
    return annotated


def run_atac_differential_accessibility(
    adata_atac: ad.AnnData,
    max_peaks: int = 30000,
) -> pd.DataFrame:
    """Compare peak accessibility between Perturbation A and B in one ATAC cluster."""
    cluster = choose_cluster_with_groups(adata_atac, "leiden_atac")
    adata_cluster = adata_atac[adata_atac.obs["leiden_atac"] == cluster].copy()
    print(f"ATAC differential accessibility cluster: {cluster}")

    x = adata_cluster.X.tocsr() if sparse.issparse(adata_cluster.X) else sparse.csr_matrix(adata_cluster.X)
    x_bin = x.copy()
    x_bin.data = np.ones_like(x_bin.data)

    groups = adata_cluster.obs["virtual_perturbation"].astype(str).to_numpy()
    mask_a = groups == "Perturbation A"
    mask_b = groups == "Perturbation B"
    n_a = int(mask_a.sum())
    n_b = int(mask_b.sum())
    if n_a < 2 or n_b < 2:
        raise ValueError("Need at least two cells in each perturbation group for ATAC testing.")

    detection = np.asarray(x_bin.sum(axis=0)).ravel()
    candidate_mask = (detection >= 5) & (detection <= (adata_cluster.n_obs - 5))
    candidate_indices = np.where(candidate_mask)[0]
    if len(candidate_indices) > max_peaks:
        variability = detection * (adata_cluster.n_obs - detection)
        keep_order = np.argsort(variability[candidate_indices])[-max_peaks:]
        candidate_indices = candidate_indices[keep_order]

    open_a = np.asarray(x_bin[mask_a, :][:, candidate_indices].sum(axis=0)).ravel()
    open_b = np.asarray(x_bin[mask_b, :][:, candidate_indices].sum(axis=0)).ravel()
    rate_a = open_a / n_a
    rate_b = open_b / n_b

    pooled = (open_a + open_b) / (n_a + n_b)
    se = np.sqrt(pooled * (1 - pooled) * ((1 / n_a) + (1 / n_b)))
    with np.errstate(divide="ignore", invalid="ignore"):
        z_scores = (rate_a - rate_b) / se
    p_values = 2 * stats.norm.sf(np.abs(z_scores))
    p_values[~np.isfinite(p_values)] = 1.0

    da = pd.DataFrame(
        {
            "peak": adata_cluster.var_names[candidate_indices],
            "cluster": cluster,
            "open_cells_A": open_a,
            "open_cells_B": open_b,
            "cells_A": n_a,
            "cells_B": n_b,
            "accessibility_rate_A": rate_a,
            "accessibility_rate_B": rate_b,
            "rate_difference_A_minus_B": rate_a - rate_b,
            "log2_rate_ratio_A_vs_B": np.log2((rate_a + 0.01) / (rate_b + 0.01)),
            "z_score": z_scores,
            "p_value": p_values,
        }
    )
    da["p_value_adj"] = bh_adjust(da["p_value"].to_numpy())
    da = annotate_da_peaks_with_coordinates(da)
    da = da.sort_values("p_value_adj")

    da_path = DATA_DIR / "atac_diff_peaks_perturbationA_vs_B.csv"
    da.to_csv(da_path, index=False)
    print(f"Saved {da_path}")
    return da


def plot_atac_volcano(da: pd.DataFrame) -> None:
    """Save a volcano plot for differential accessibility results."""
    plot_df = da.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["log2_rate_ratio_A_vs_B", "p_value_adj"]
    ).copy()
    plot_df["minus_log10_padj"] = -np.log10(np.clip(plot_df["p_value_adj"], 1e-300, 1))

    plt.figure(figsize=(6, 5))
    sns.scatterplot(
        data=plot_df,
        x="log2_rate_ratio_A_vs_B",
        y="minus_log10_padj",
        s=12,
        alpha=0.5,
        edgecolor=None,
    )
    plt.axvline(0, color="black", linewidth=0.8)
    plt.xlabel("log2 accessibility rate ratio, Perturbation A versus B")
    plt.ylabel("-log10 adjusted p value")
    plt.title("ATAC differential accessibility")

    top = plot_df.head(10)
    for _, row in top.iterrows():
        label = row["peak"]
        plt.text(row["log2_rate_ratio_A_vs_B"], row["minus_log10_padj"], label, fontsize=6)

    save_figure("atac_diff_peaks_perturbationA_vs_B_volcano.png")


def add_peak_set_score(adata_atac: ad.AnnData, peaks: list[str], score_name: str) -> None:
    """Score cells by the fraction of selected peaks that are accessible."""
    present = [peak for peak in peaks if peak in adata_atac.var_names]
    if not present:
        adata_atac.obs[score_name] = 0.0
        return

    x = adata_atac[:, present].X
    if sparse.issparse(x):
        score = np.asarray((x > 0).mean(axis=1)).ravel()
    else:
        score = np.asarray((x > 0).mean(axis=1)).ravel()
    adata_atac.obs[score_name] = score


def add_atac_program_scores(adata_atac: ad.AnnData, da: pd.DataFrame, n_peaks: int = 150) -> None:
    """Create simple ATAC scores from peaks that are more open in each group."""
    da_valid = da.dropna(subset=["p_value_adj", "log2_rate_ratio_A_vs_B"]).copy()
    top_a = da_valid[da_valid["log2_rate_ratio_A_vs_B"] > 0].sort_values("p_value_adj").head(n_peaks)
    top_b = da_valid[da_valid["log2_rate_ratio_A_vs_B"] < 0].sort_values("p_value_adj").head(n_peaks)

    add_peak_set_score(adata_atac, top_a["peak"].tolist(), "atac_A_open_peak_score")
    add_peak_set_score(adata_atac, top_b["peak"].tolist(), "atac_B_open_peak_score")


def plot_multiome_heatmap(adata_rna: ad.AnnData, adata_atac: ad.AnnData) -> pd.DataFrame:
    """Save a heatmap that links RNA and ATAC perturbation responses."""
    summary = pd.DataFrame(index=adata_rna.obs_names)
    summary["virtual_perturbation"] = adata_rna.obs["virtual_perturbation"].astype(str)
    summary["RNA neuroinflammation score"] = adata_rna.obs["neuroinflammation_score"]
    summary["RNA synaptic score"] = adata_rna.obs["synaptic_score"]
    summary["ATAC A-open peak score"] = adata_atac.obs.loc[summary.index, "atac_A_open_peak_score"]
    summary["ATAC B-open peak score"] = adata_atac.obs.loc[summary.index, "atac_B_open_peak_score"]

    heatmap_data = summary.groupby("virtual_perturbation").mean()
    z_data = (heatmap_data - heatmap_data.mean(axis=0)) / heatmap_data.std(axis=0).replace(0, np.nan)
    z_data = z_data.fillna(0)

    plt.figure(figsize=(8, 3.5))
    sns.heatmap(z_data, cmap="vlag", center=0, annot=heatmap_data.round(3), fmt="", linewidths=0.5)
    plt.title("RNA and ATAC response scores by virtual perturbation")
    plt.xlabel("score")
    plt.ylabel("")
    save_figure("multiome_virtual_perturbation_pathway_heatmap.png")
    return heatmap_data


def main() -> None:
    """Run the full analysis."""
    sns.set_theme(style="whitegrid", context="notebook")
    sc.settings.verbosity = 2

    adata = load_multiome_h5()
    adata_rna, adata_atac = split_modalities(adata)

    adata_rna = filter_rna(adata_rna)
    adata_atac = filter_atac(adata_atac)
    adata_rna, adata_atac = align_cells(adata_rna, adata_atac)
    plot_qc(adata_rna, adata_atac)

    adata_rna = preprocess_rna(adata_rna)
    adata_atac = preprocess_atac(adata_atac)

    plot_umap(adata_rna, "leiden_rna", "rna_umap_clusters.png", "RNA clusters")
    plot_umap(adata_atac, "leiden_atac", "atac_umap_lsi_clusters.png", "ATAC LSI clusters")

    adata_rna, adata_atac = add_virtual_perturbations(adata_rna, adata_atac)

    plot_umap(
        adata_rna,
        "virtual_perturbation",
        "rna_umap_virtual_perturbations.png",
        "RNA virtual perturbations",
    )
    plot_umap(
        adata_atac,
        "virtual_perturbation",
        "atac_umap_virtual_perturbations.png",
        "ATAC virtual perturbations",
    )

    plot_rna_pathway_scores(adata_rna)
    de = run_rna_differential_expression(adata_rna)
    plot_rna_volcano(de)

    da = run_atac_differential_accessibility(adata_atac)
    plot_atac_volcano(da)

    add_atac_program_scores(adata_atac, da)
    heatmap_data = plot_multiome_heatmap(adata_rna, adata_atac)
    heatmap_path = DATA_DIR / "multiome_virtual_perturbation_summary.csv"
    heatmap_data.to_csv(heatmap_path)
    print(f"Saved {heatmap_path}")

    adata_rna.write(DATA_DIR / "processed_rna_virtual_crispr.h5ad")
    adata_atac.write(DATA_DIR / "processed_atac_virtual_crispr.h5ad")
    print("Pipeline finished.")


if __name__ == "__main__":
    main()


import subprocess
from pathlib import Path

def run_motif_scanning(atac_da, data_dir, group="A", top_n=500):
    """
    Write top differential peaks to BED and run HOMER findMotifsGenome.pl.
    group: 'A' for peaks more open in Perturbation A, 'B' for Perturbation B
    """
    data_dir = Path(data_dir)
    
    # Filter top peaks for the chosen group
    if group == "A":
        peaks = atac_da[atac_da["log2_rate_ratio_A_vs_B"] > 0].nsmallest(top_n, "p_value")
    else:
        peaks = atac_da[atac_da["log2_rate_ratio_A_vs_B"] < 0].nsmallest(top_n, "p_value")
    
    # Write BED file (HOMER needs chr, start, end, name, score, strand)
    bed_path = data_dir / f"perturbation_{group}_peaks.bed"
    peaks[["chrom", "start", "end", "peak"]].assign(
        score=0, strand="."
    ).to_csv(bed_path, sep="\t", header=False, index=False)
    
    # Run HOMER
    out_dir = data_dir / f"homer_motifs_perturbation_{group}"
    out_dir.mkdir(exist_ok=True)
    
    cmd = [
        "findMotifsGenome.pl",
        str(bed_path),
        "mm10",
        str(out_dir),
        "-size", "given",   # use actual peak coordinates
        "-mask",            # mask repeats
        "-p", "4"           # 4 CPU threads
    ]
    subprocess.run(cmd, check=True)
    print(f"HOMER results saved to {out_dir}")
    return out_dir