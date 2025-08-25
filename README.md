# IntersectOmics

**IntersectOmics** is a computational framework for analyzing multi-omics datasets with time series or multiple conditions. It identifies **concordant** biomolecules—those that behave similarly across different omics layers (e.g., transcriptomics, proteomics, metabolomics)—using correlation, graph theory, and community detection.

This tool is ideal for detecting groups of genes, proteins, and metabolites that respond in a coordinated manner across different biological layers and time points. While future updates may support **discordant** relationships (opposite behaviors across layers), this version focuses **only on concordant biomolecules**.

---

## Installation

After cloning the repository, install the package in editable mode:

```bash
pip install -e .
```

---

## Supported Dataset

IntersectOmics supports multi-omics datasets with any number of replicates and experimental conditions. Each omics layer (transcriptomics, proteomics, metabolomics, etc.) should be provided as a separate table.

### Input Data Structure

Each input dataset should be structured as follows:

- **Rows**: Unique biomolecule identifiers (e.g., gene/protein/metabolite names)
- **Columns**: Sample measurements, ideally grouped by condition and replicate

### Accepted Formats

- **Flat columns** (e.g., `Day1_Rep1`, `Day2_Rep3`, etc.)
- **Multi-index columns**, where hierarchical levels (e.g., timepoint, replicate) are stored in tuples (e.g., `("Day1", "Rep1")`)

### Requirements

- All omics layers must use consistent biomolecule identifiers to enable graph intersection.
- Replicates must be distinguishable by naming convention or metadata.
- Handle missing values appropriately before using the tool.

### Example Table

| Biomolecule | Day1_Rep1 | Day1_Rep2 | Day1_Rep3 | Day2_Rep1 | Day2_Rep2 | Day2_Rep3 |
|-------------|-----------|-----------|-----------|-----------|-----------|-----------|
| GeneA       | 1.2       | 1.3       | 1.1       | 1.5       | 1.4       | 1.6       |
| GeneB       | 0.5       | 0.4       | 0.6       | 0.3       | 0.2       | 0.4       |
| GeneC       | 2.1       | 2.2       | 2.0       | 1.9       | 1.8       | 2.0       |

### Real Dataset Example

We use a dataset from [this publication](https://www.sciencedirect.com/science/article/pii/S0048969723003558), which measured transcriptomics and proteomics in springtail earthworms over time after insecticide exposure.

<p align="center">
  <img width="400" alt="example_input_data" src="https://github.com/Melclic/intersectomics/assets/4260862/5265380e-c6e9-4969-babb-cbd9dc882832">
</p>

---

## Correlation with Replicates

Standard correlation methods often require averaging across replicates, discarding important information about variability. IntersectOmics instead fits a distribution to each condition and draws bootstrap samples to estimate correlations more robustly.

<p align="center">
  <img width="700" alt="example_bootstrap" src="https://github.com/Melclic/intersectomics/assets/4260862/68140855-ad3c-43a2-ba5d-ad0643e8169b">
</p>

### Bootstrap Workflow

1. Fit a normal distribution at each time point using replicate values.
2. Sample one value per time point from the fitted distribution.
3. Repeat this sampling process *n* times to compute a distribution of correlation values.
4. Average these correlations.
5. Combine p-values using the Pearson method (`scipy.stats.combine_pvalues`).

> Note: While a normal distribution is currently used, future versions may include data-type-specific distributions (e.g., Poisson for RNA-seq).

---

## Supported Correlation Metrics

- **Spearman** (default): Ideal for rank-based, monotonic trends, especially in time series.
- **Pearson**: *TODO*
- **Euclidean Distance**: *TODO*

---

## Graph Construction

For each omics layer, a graph is created:

- **Nodes**: Biomolecules
- **Edges**: Similarity score (correlation) between biomolecules
- **Weights**: Mean correlation across bootstrap iterations

Only **concordant** edges (positive correlations) are retained in the graph.

<p align="center">
  <img width="800" alt="protein_spearman_graph" src="https://github.com/Melclic/intersectomics/assets/4260862/29b31e32-a2e6-4feb-a2ef-cad67ac219a8">
</p>

---

### Why We Ignore Discordant Edges

The current goal is to detect concordant biomolecules—those showing similar directional changes across omics layers. Including discordant edges (e.g., negative correlations) can lead to misleading results and incorrect community grouping.

<p align="center">
  <img width="400" alt="anticorrelation_graph" src="https://github.com/Melclic/intersectomics/assets/4260862/a4f0e411-d01b-4f86-a36e-118755195180">
</p>

<p align="center">
  <img width="500" alt="anticorrelation_mistake" src="https://github.com/Melclic/intersectomics/assets/4260862/19f93788-4426-4bad-af7f-dc9c4d0b06fd">
</p>

---

## Graph Intersection

We create one graph per omics layer, then compute the **intersection**:

- **Nodes**: Must exist in all graphs (matched by name)
- **Edges**: Retained only if present in all omics-specific graphs

The result is a unified graph of biomolecules that are **concordantly** related in every omics layer.

<p align="center">
  <img width="500" alt="graph_intersection" src="https://mathworld.wolfram.com/images/eps-svg/GraphIntersection_800.svg">
</p>

---

## Community Detection

We apply community detection to the intersected graph to identify clusters of biomolecules with shared behavior.

- These clusters represent **concordant groups**
- Edges are weighted by correlation strength
- Each community may contain genes, proteins, and metabolites

<p align="center">
  <img width="800" alt="G_inter_example" src="https://github.com/Melclic/intersectomics/assets/4260862/27fab3fd-fd73-4ce9-9f7b-010d399ffc55">
</p>

---

## Final Results

The result is a set of **concordant communities**—biomolecules across different omics layers that behave similarly over time or conditions.

<p align="center">
  <img width="500" alt="result" src="https://github.com/Melclic/intersectomics/assets/4260862/1495f700-da58-4d75-8347-1d43ba1d10bd">
</p>

> **Note:** Detection of **discordant** communities (e.g., genes increasing while proteins decrease) is not currently supported, but is planned for future versions.

---

## Inspiration

IntersectOmics is inspired by the work of [Nikolay Oskolkov](https://github.com/NikolayOskolkov) on [UMAPDataIntegration](https://github.com/NikolayOskolkov/UMAPDataIntegration), with key extensions for:

- Time series data
- Replicate-aware correlation
- Graph intersection
- Concordant community analysis