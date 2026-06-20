# CRISPR-like scRNA + scATAC Multiome Analysis Pipeline

A conceptual CRISPR-like single-cell multiome analysis pipeline applied to a 10x Genomics mouse brain Alzheimer's disease dataset. This project simulates the downstream analytical workflow of a Perturb-seq or CRISPR-sciATAC experiment, combining scRNA-seq and scATAC-seq modalities to identify transcription factor programs and gene regulatory networks disrupted in neuroinflammatory states.

---

## Project Overview

This pipeline demonstrates the core computational skills required for multiome perturbation analysis:

- **Dual-modality preprocessing**: joint quality control, filtering, and normalization of RNA and ATAC data from the same cells
- **Virtual perturbation labelling**: cells are stratified into two groups (Perturbation A: high neuroinflammation score; Perturbation B: high synaptic score) to simulate a CRISPR knockin/knockout contrast
- **Differential gene expression**: Wilcoxon rank-sum testing with Benjamini-Hochberg correction across 30,000+ cells
- **Differential chromatin accessibility**: rate-ratio testing across ATAC peaks between virtual perturbation groups
- **Transcription factor motif scanning**: HOMER `findMotifsGenome.pl` on differential peaks against the mm10 genome to identify candidate TF regulators
- **Cross-modality integration**: RNA pathway scores and ATAC peak scores jointly visualized to confirm regulatory concordance

---

## Dataset

| Property | Value |
|---|---|
| Source | [10x Genomics Mouse Brain Multiome (Alzheimer's AppNote)](https://www.10xgenomics.com/datasets/fresh-cortex-from-alzheimer-s-disease-mouse-model-multiplexed-samples-2-standard) |
| Species | *Mus musculus* (mm10) |
| Modalities | scRNA-seq + scATAC-seq (same cells) |
| Cells (post-QC) | ~30,000 |
| RNA features | ~32,000 genes |
| ATAC features | ~110,000 peaks |

---

## Key Results

### 1. RNA UMAP, Virtual Perturbation Labels

Cells separate into neuroinflammation-high (Perturbation A, blue) and synaptic-high (Perturbation B, orange) populations across multiple clusters, consistent with a glial vs. neuronal cell state axis.

![RNA UMAP Virtual Perturbations](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/rna_umap_virtual_perturbations.png?raw=true)

### 2. RNA Clusters

Ten Leiden clusters resolved across the full transcriptome UMAP. Cluster identity corresponds to distinct brain cell types including microglia/astrocytes (neuroinflammation-enriched) and excitatory/inhibitory neurons (synaptic-enriched).

![RNA Clusters](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/rna_umap_clusters.png?raw=true)

### 3. RNA Quality Control Metrics

Post-QC distributions show a clean library: RNA counts peak below 25,000 per cell, gene detection spans 500-8,000 genes per cell, and mitochondrial fraction is tightly controlled below 5% in most cells.

![RNA QC Metrics](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/rna_qc_metrics.png?raw=true)

### 4. RNA Pathway Scores by Virtual Perturbation

Perturbation A cells score significantly higher on the neuroinflammation gene set, while Perturbation B cells score higher on the synaptic gene set, confirming that the virtual perturbation labels capture a biologically meaningful separation.

![RNA Pathway Scores](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/rna_pathway_scores_virtual_perturbations.png?raw=true)

### 5. RNA Differential Expression Volcano

*Apoe* is the most significantly upregulated gene in Perturbation A (-log10 p > 240), followed by *AY036118*, *Cst3*, *Ckb1*, *Calm1*, and *B2m*, all established markers of activated microglia and disease-associated glial states in Alzheimer's models.

![RNA Differential Expression Volcano](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/rna_de_perturbationA_vs_B_volcano.png?raw=true)

### 6. ATAC UMAP, Virtual Perturbation Labels

The ATAC chromatin landscape separates Perturbation A and B cells into distinct spatial regions, confirming that the transcriptional differences are paralleled by differential chromatin accessibility.

![ATAC UMAP Virtual Perturbations](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/atac_umap_virtual_perturbations.png?raw=true)

### 7. ATAC LSI Clusters

Ten LSI-based chromatin clusters reveal clear cell-type-specific accessibility patterns, structurally consistent with the RNA cluster topology.

![ATAC LSI Clusters](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/atac_umap_lsi_clusters.png?raw=true)

### 8. ATAC Quality Control Metrics

Fragment counts are broadly distributed between 1,000-10,000 per cell and peaks per cell between 500-5,000, consistent with high-quality ATAC-seq libraries.

![ATAC QC Metrics](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/atac_qc_metrics.png?raw=true)

### 9. ATAC Differential Accessibility Volcano

Top differentially accessible peaks include `chr1:24612324-24613211` and `chr7:139900547-139901380` (more open in Perturbation A) and several peaks on chr5 and chr1 (more open in Perturbation B), spanning ~1.5 log2 accessibility ratio.

![ATAC Differential Accessibility Volcano](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/atac_diff_peaks_perturbationA_vs_B_volcano.png?raw=true)

### 10. RNA + ATAC Cross-Modality Response Heatmap

The heatmap confirms concordance between RNA and ATAC modalities: Perturbation A shows a positive neuroinflammation RNA score (+0.132) and negative synaptic RNA score (-0.147), while Perturbation B shows the inverse pattern.

![Cross-Modality Heatmap](https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome/blob/main/images/multiome_virtual_perturbation_pathway_heatmap.png?raw=true)

---

## Biological Interpretation

The RNA and ATAC results converge on a consistent story across three levels of evidence:

1. **Cell state (UMAP)**: Perturbation A and B cells separate in both RNA and ATAC space
2. **Gene expression (DE)**: Perturbation A upregulates *Apoe*, *B2m*, *Cst3*, markers of activated microglia in Alzheimer's disease
3. **Regulatory chromatin (HOMER)**: Perturbation A peaks are bound by AP-1 and BATF (DAM program drivers); Perturbation B peaks are bound by RFX and CTCF (neuronal regulatory elements)

In a real CRISPR-sciATAC experiment, this workflow would identify which CRISPR guide RNA targets shift cells from the AP-1-driven inflammatory state toward the RFX-driven neuronal state, which is a directly testable mechanistic hypothesis for Alzheimer's disease therapeutic research.

---

## Repository Structure

```text
CRISPR-scRNA-scATAC-Multiome/
├── crispr_multiome_pipeline.py
├── CRISPR_scRNA_scATAC_analysis.ipynb
├── images/
│   ├── rna_umap_virtual_perturbations.png
│   ├── rna_umap_clusters.png
│   ├── rna_qc_metrics.png
│   ├── rna_pathway_scores_virtual_perturbations.png
│   ├── rna_de_perturbationA_vs_B_volcano.png
│   ├── atac_umap_virtual_perturbations.png
│   ├── atac_umap_lsi_clusters.png
│   ├── atac_qc_metrics.png
│   ├── atac_diff_peaks_perturbationA_vs_B_volcano.png
│   └── multiome_virtual_perturbation_pathway_heatmap.png
└── README.md
```

> **Note:** The raw 10x Genomics H5 data file is not included due to size. Download it from the [10x Genomics website](https://www.10xgenomics.com/datasets/fresh-cortex-from-alzheimer-s-disease-mouse-model-multiplexed-samples-2-standard) and place it in the project root before running the pipeline.

---

## Installation and Usage

### Prerequisites

```bash
conda create -n ad_multiome python=3.11 -y
conda activate ad_multiome
conda install -c conda-forge scanpy anndata scipy scikit-learn seaborn matplotlib -y
conda install -c bioconda homer -y

# Download mm10 genome for HOMER motif scanning
perl ~/miniforge3/envs/ad_multiome/share/homer/configureHomer.pl -install mm10
```

### Running the Pipeline

```bash
git clone https://github.com/chuabrandon26/CRISPR-scRNA-scATAC-Multiome.git
cd CRISPR-scRNA-scATAC-Multiome
python crispr_multiome_pipeline.py
```

Or step through interactively:

```bash
jupyter notebook CRISPR_scRNA_scATAC_analysis.ipynb
```

---

## Key Dependencies

| Package | Version | Purpose |
|---|---|---|
| `scanpy` | >=1.10 | scRNA-seq preprocessing, PCA, UMAP, clustering |
| `anndata` | >=0.10 | Single-cell data container |
| `scipy` | >=1.11 | Statistical tests (Wilcoxon rank-sum) |
| `scikit-learn` | >=1.3 | Truncated SVD for LSI (ATAC) |
| `seaborn` / `matplotlib` | latest | Visualization |
| `HOMER` | 5.1 | TF motif enrichment on differential ATAC peaks |

---

## Skills Demonstrated

- **Single-cell genomics**: dual-modality preprocessing (scRNA + scATAC), QC filtering, normalization, dimensionality reduction (PCA, LSI/TF-IDF)
- **CRISPR perturbation analysis**: virtual perturbation labelling, differential expression, differential accessibility
- **Chromatin biology**: ATAC-seq peak analysis, TF motif enrichment with HOMER, mm10 genome annotation
- **Statistical methods**: Wilcoxon rank-sum, Benjamini-Hochberg FDR correction, rate-ratio testing
- **Python bioinformatics stack**: Scanpy, AnnData, SciPy, scikit-learn, Pandas, NumPy
- **Workflow engineering**: modular pipeline design, reproducible analysis, cross-modality integration

---

## License

MIT License. See [LICENSE](LICENSE) for details.

---

*Built as part of a personal portfolio to develop skills in CRISPR functional genomics and single-cell multiome analysis.*
