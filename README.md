# IntersectOmics
Method to analyze multi-omics datasets time series (or multiple condition) and uncover 
similarly behaving groups of biomoelcules.

To be able to combine the different omics layers, we need to convert each data
to a non-parametric space. In other words, we need to find a way to remove the 
"memory" of the measuring process of each data type. To that end we perform a 
pairwise comparison between each pair of each data type and then construct a graph.
We use the metric of similarity as an edge value and then perform community analysis
using the weights to find cluster of biomolecules that have similar behaviors. 

## Multi-omics data

We will use a dataset from the [following paper](https://www.sciencedirect.com/science/article/pii/S0048969723003558), 
that measured transcriptomics and proteomics of springtail earth worm over several days after exposure to an insecticide.

The time series data has multiple omics layers and has three replicates
for each time point. Note that the data needs to be in the form of a table, where the index
of the table are the names of the genes/protein/metabolite (biomolecule) and the columns represents
the metadata associated with the sample.  

<img width="200" alt="example_input_data" src="https://github.com/Melclic/intersectomics/assets/4260862/5265380e-c6e9-4969-babb-cbd9dc882832">


### Correlation with replicates

When performing correlation analysis with replicates you are forced to take the
mean of the samples. By doing so, you loose information regarding the variability of the
sample. This may lead to false negative results when performing correlation analysis.
We propose a method that is more robust but more computationally intensive by bootstrapping
random variables extracted from fitted distributions for each replicate.

To this end, this package enables the user to perform bootstrap analysis on two
vectors that contain known replicates. The process may be summarized as follows:
1) Given two vectors with replicates, fit a normal distribution at each time points
2) Loop n times and each time, sample a single value at each timepoint
3) Calculate the  correlation between the two
4) Take the mean of the correlations and combine the p-values correcting for multiple
tesing

The reported p-value is combined using the pearson method. See the scipy documentation for `combined_pvalues`.

TODO: As of now, we use a normal distribution, but more appropriate distributions
should be used depending on the data type. For example, RNA-seq should use Poisson
distribution instead. To that end, we should run R packages to process the data and
extract better distribution parameters

Below are the supported correlation types:

#### Spearman

The default. This works very well when parameters have curvilinear relationship.
In our example dataset that is time series, an increase could mean a decrease in 
another. Spearman correlation is the most appropriate method.

#### Pearson

TODO

#### Euclidian

TODO

#### 

### Turning the results to a graph

The graph represents the pairwise similaritly between each biomolecule for each 
omics layer. For example, in the example we have a transcriptomics graph and a 
protemics graph. Each node is either a gene or a protein and each edge represents 
the correlation score.

<p align="center">
  <img width="200" alt="protein_spearman_graph" src="https://github.com/Melclic/intersectomics/assets/4260862/29b31e32-a2e6-4feb-a2ef-cad67ac219a8">
</p>

#### Ignoring the anti-correlation

The goal of the analysis is to find collection of genes that behave similarly. 
If anticorrelations and correlations are used to contruct the graph, you would get
local subgrpahs that are highly connected between biomolecules, but you would not be 
able to distinguish between those that behave similarly and those that do not. Below is 
a small example of three biomolecules that are anticorrelated to each other. By nature
of the 

<p align="center">
  <img width="200" alt="anticorrelation_graph" src="https://github.com/Melclic/intersectomics/assets/4260862/a4f0e411-d01b-4f86-a36e-118755195180">
</p>

<p align="center">
  <img width="200" alt="anticorrelation_mistake" src="https://github.com/Melclic/intersectomics/assets/4260862/19f93788-4426-4bad-af7f-dc9c4d0b06fd">
</p>

### Graph intersection

Now that we have multiple graphs for each omics layer, we combine them by taking
the interection between each graph. This means that we keep only an edge if it exists
in all of the graphs. Note that the nodes also need to be present in each of the graphs
and must have the same names. If not, the result will be orphan nodes that will be removed 
from the resulting graph.

<p align="center">
  <img width="200" alt="graph_intersection" src="https://mathworld.wolfram.com/images/eps-svg/GraphIntersection_800.svg">
</p>

### Community Analysis

Now that we have a consensus graph, we need can analyze the results and extract 
groups of omics layers that behave similarly. Note that each layer may not behave 
the same, but each group would

To that end we use community analysis to detect groups of nodes that are well 
connected. We use the correlation metric of our choice as a numerical value of 
closeness between the two. 

<p align="center">
  <img width="500" alt="G_inter_example" src="https://github.com/Melclic/intersectomics/assets/4260862/27fab3fd-fd73-4ce9-9f7b-010d399ffc55">
</p>

### Result

The result are collection of proteins, genes, and metabolites that are grouped together because 
they behave the same. Note that in the example below, the protein and genes behave the same over 
time, but this is not always the case. Here is an example of a single community in the above graph.

<p align="center">
  <img width="200" alt="result" src="https://github.com/Melclic/intersectomics/assets/4260862/1495f700-da58-4d75-8347-1d43ba1d10bd">
</p>

## Inspiration

The idea to convert multiple data types to a non-parametric space
and perform an intersection study has been inspired by [Nikolay Oskolkov](https://github.com/NikolayOskolkov)
and the following [github](https://github.com/NikolayOskolkov/UMAPDataIntegration).
The added features are limitations I have found when implementing the method with
time series data with replicates.
