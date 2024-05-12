# IntersectOmics
Method to study multi-omics datasets, measured with replicates and find those 
that have similar patterns of behaviors.

## Multi-omics data

We will use a dataset from the [following paper](https://www.sciencedirect.com/science/article/pii/S0048969723003558), 
that measured transcriptomics and proteomics of an earth worm  over several days after exposure to an insecticide.

#paste an example picture of the data

This time series data contains all the elements we need. It is multi-omics and has replicates
for each time point. Note that the data needs to be in the form of a table, where the index
of the table are the names of the genes/protein/metabolite and the columns represents
the metadata associated with the sample. We have a helper function that can combine
the two:

python```
#shoe the example code
```

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

The reported p-value is combined using the pearson method. See the scipy documentation for combined_pvalues.

TODO: As of now, we use a normal distribution, but more appropriate distributions
should be used depending on the data type. For example, RNA-seq should use Poisson
distribution instead. To that end, we should run R packages to process the data and
extract better distribution parameters

Below are the supported correlationtyps

#### Spearman

The default 

#### Pearson

#### Euclidian

#### 

### Turning the results to a graph

#### Ignoring the anti-correlation

The goal of the analysis is to find collection of genes that behave similarly,

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
and perform an intersection study between the two. 

##
