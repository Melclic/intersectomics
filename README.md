# IntersectOmics
Method to study multi-omics datasets, measured with replicates and find those 
that have similar patterns of expression over time or different measuring 
types.

## Multi-omics data

We will use a dataset from the [following paper](https://www.sciencedirect.com/science/article/pii/S0048969723003558), 
that measured transcriptomics and proteomics of springtail earth worm over several days after exposure to an insecticide.

The time series data has multiple omics layers and has three replicates
for each time point. Note that the data needs to be in the form of a table, where the index
of the table are the names of the genes/protein/metabolite and the columns represents
the metadata associated with the sample. 

The data 

### Correlation with replicates

When performing correlation analysis with replicates you are forced to take the
mean of the samples. Doing so looses information regarding the variability of the
sample and may lead to false negative results when performing correlation analysis.
We propose a method that is more robust but more computationally intensive.

#show an example

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
another.

#### Pearson

TODO

#### Euclidian

TODO

#### 

### Turning the results to a graph

#### Ignoring the anti-correlation

The goal of the analysis is to find collection of genes that behave similarly. 
When constructing the graph, if you use anticorrelations then you would not get 
the results you expect. Below is an example of the plot of anticorrelations graph.
The first plot is anticorrelated to the second. The second is anticorrelated to the 
third, but the first is correlated to the third. This means as a graph it would not be 
useful.


![anticorrelation_mistake](https://github.com/Melclic/intersectomics/assets/4260862/19f93788-4426-4bad-af7f-dc9c4d0b06fd)


To be able to combine the different omics layers, we need to convert each data
to a non-parametric space. In other words, we need to find a way to remove the 
"memory" of the measuring process of each data type. To that end we perform a 
pairwise comparison between each pair of each data type and then construct a graph.
We use the metric of similarity as an edge value.

### Graph intersection

Now that we have multiple graphs for each omics layer, we combine them by taking
the interection between the them. This means that we keep only an edge if it exists
in all of the graphs. Note that the nodes also need to be present in each of the graphs.
If not, the result will be orphan nodes that will be removed from the resulting graph.

### Community Analysis

Now that we have a consensus graph, we need can analyze the results and extract 
groups of omics layers that behave similarly. Note that each layer may not behave 
the same, but each group would

To that end we use community analysis to detect groups of nodes that are well 
connected. We use the correlation metric of our choice as a numerical value of 
closeness between the two. 

## Inspiration

The idea of this method is to convert multiple data types to a non-parametric space
and perform an intersection study has been inspired by [Nikolay Oskolkov](https://github.com/NikolayOskolkov)
and the following [github](https://github.com/NikolayOskolkov/UMAPDataIntegration).
The added features are limitations I have found when implementing the method with
different types of data.
