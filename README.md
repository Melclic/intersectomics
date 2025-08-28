# IntersectOmics

**IntersectOmics** is a computational framework for analyzing multi-omics datasets with time series or multiple conditions. It identifies biomolecules (e.g., genes, proteins, metabolites) that behave in a **coordinated** manner—either **concordantly** (similar direction) or **discordantly** (opposite direction)—across different omics layers, using correlation, graph theory, and community detection.

This tool is ideal for uncovering complex biological patterns across layers and conditions. The interpretation of whether clusters represent concordant or discordant behavior is **left to the user**, based on downstream inspection or visualization.

> **TODO:** Add support for automatic detection and labeling of concordant vs. discordant biomolecule relationships across omics layers along with biological explanation as detailed (here)[https://substack.com/history/post/157552707]

---

## Installation

After cloning the repository, install the package in editable mode:

```bash
pip install -e .
```

---

## Supported Dataset

IntersectOmics supports multi-omics datasets with any number of replicates and experimental conditions. Each omics layer (transcriptomics, proteomics, metabolomics, etc.) should be provided as a separate table. To ensure consistency between the metadata and data, the program requires columns that are multiindex. 

### Input Data Structure

Each input dataset should be structured as follows:

- **Rows**: Unique biomolecule identifiers (e.g., gene/protein/metabolite names)
- **Columns**: Sample measurements, ideally grouped by condition and replicate as multiindex

You need at least two different omics layers.

### Requirements

- All omics layers must use consistent biomolecule identifiers to enable graph intersection.
- Replicates must be distinguishable by naming convention or metadata.
- Handle missing values appropriately before using the tool.

---

## Example Dataset

We use a dataset from [this publication](https://www.sciencedirect.com/science/article/pii/S0048969723003558), which measured transcriptomics and proteomics in springtail earthworms over time after insecticide exposure.

<p align="center">
  <img width="400" alt="example_input_data" src="https://github.com/Melclic/intersectomics/assets/4260862/5265380e-c6e9-4969-babb-cbd9dc882832">
</p>

---

## Correlation with Replicates

Rather than averaging replicates—which can lose valuable variance information—IntersectOmics fits a distribution for each condition and bootstraps correlation values by sampling from these distributions.

<p align="center">
  <img width="700" alt="example_bootstrap" src="https://github.com/Melclic/intersectomics/assets/4260862/68140855-ad3c-43a2-ba5d-ad0643e8169b">
</p>

### Bootstrap Workflow

1. Fit a normal distribution at each time point using replicate values.
2. Sample one value per time point from the fitted distribution and perform correlation computation
3. Repeat this sampling process *n* times to compute a distribution of correlation values.
4. Average these correlations.
5. Combine p-values using the Pearson method (`scipy.stats.combine_pvalues`).

> Note: While a normal distribution is currently used, future versions may support data-type-specific distributions (e.g., Poisson for RNA-seq).

---

## Supported Correlation Metrics

- **Spearman** (default): Rank-based, ideal for monotonic or curvilinear trends
- **Pearson**: *TODO*
- **Euclidean Distance**: *TODO*

---

## Graph Construction

A graph is constructed for each omics layer. Note that by default any correlations that have a significance <=0.05 is ignore, but that threshold can be modified by the user:

- **Nodes**: Biomolecules
- **Edges**: Pairwise similarity scores between biomolecules
- **Weights**: Averaged correlation scores (from bootstrapped sampling)

Edges can reflect both **positive (concordant)** and **negative (discordant)** relationships.

<p align="center">
  <img width="800" alt="protein_spearman_graph" src="https://github.com/Melclic/intersectomics/assets/4260862/29b31e32-a2e6-4feb-a2ef-cad67ac219a8">
</p>

---

## Graph Intersection

Once a graph is built for each omics layer, their **intersection** is computed:

- **Nodes**: Must exist in all graphs
- **Edges**: Retained only if present in all graphs

This ensures only biomolecule relationships that are consistent across all omics layers are preserved, regardless of the direction of the correlation.

<p align="center">
  <img width="500" alt="graph_intersection" src="https://mathworld.wolfram.com/images/eps-svg/GraphIntersection_800.svg">
</p>

---

## Community Detection

Community detection is performed on the intersected graph:

- Identifies clusters of biomolecules that are similarly related across layers
- Uses edge weights (correlation) as a measure of connectivity

> **Important:** Communities may include biomolecules that are concordantly or discordantly related between omics layers. **It is the user’s responsibility to inspect each community** (e.g., via visualization or trend comparison) to interpret the biological meaning.

<p align="center">
  <img width="800" alt="G_inter_example" src="https://github.com/Melclic/intersectomics/assets/4260862/27fab3fd-fd73-4ce9-9f7b-010d399ffc55">
</p>

---

## Final Results

Each community represents a group of biomolecules with shared trends across layers. These may be:

- **Concordant**: Biomolecules change in the same direction
- **Discordant**: Biomolecules show opposing trends

<p align="center">
  <img width="500" alt="concordant_results" src="https://github.com/user-attachments/assets/fbafb051-2357-4bf2-a868-53830137cb7b" />
</p>

<p align="center">
  <img width="500" alt="discordant_results" src="https://github.com/user-attachments/assets/2ab6db9d-d28c-4d97-839d-6c33e7b0a384" />
</p>

Use visualization tools to evaluate and annotate the nature of each cluster.

--- 


## Notebooks

A complete example demonstrating the full pipeline, including data loading, correlation, graph construction, and visualization, can be found at:

`notebooks/run_example.ipynb`

---

## Inspiration

IntersectOmics is inspired by [Nikolay Oskolkov](https://github.com/NikolayOskolkov)’s [UMAPDataIntegration](https://github.com/NikolayOskolkov/UMAPDataIntegration), with major additions including:

- Support for time series data
- Bootstrap-based correlation with replicates
- Cross-layer graph intersection
- Discovery of both concordant and discordant multi-omics communities
