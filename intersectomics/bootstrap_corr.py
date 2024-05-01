import pandas as pd
import numpy as np
from typing import  Dict, Any, Tuple
from scipy.stats import spearmanr, norm, combine_pvalues


def _calc_norm_dict(
    sample_series: pd.Series
) -> Dict[any, norm]:
    """
    Calculate the normal distributions for each unique index value in a Series.
    
    The function resets the index of the series, computes the mean and standard deviation
    for each unique original index, and returns a dictionary mapping each index to a
    normal distribution (`scipy.stats.norm`) with the corresponding mean and standard deviation.
    
    Args:
        series (Series): A pandas Series with a named index and numerical values.
        
    Returns:
        Dict[any, norm]: A dictionary where each key is an index from the original Series and each value
                         is a `scipy.stats.norm` object parameterized by the mean and standard deviation
                         of the values associated with that index.
    """
    df = sample_series.reset_index()
    # Calculate the mean and standard deviation for each replicate
    means = df.groupby(by=sample_series.index.name)[sample_series.name].mean().to_dict()
    stds = df.groupby(by=sample_series.index.name)[sample_series.name].std().to_dict()
    # Return the normal distributions given the replicates
    return {i: norm(loc=means[i], scale=stds[i]) for i in means}


def bootstrap_spearman_correlation(
    sample_series_1: pd.Series,
    sample_series_2: pd.Series,
    n_iterations: int = 10,
) -> Tuple[float, Tuple[float, float], float]:
    """
    Perform a bootstrap Spearman correlation analysis between two genes over time.
    Usumes that the replicate distributions are normally distributed.

    TODO: need to add native sample distributions for each omics level

    Combined pvalues are performed using the pearson method (see scipy method)

    Args:
        sample_series_1 (pd.Series): Series of the first sample. The index will be used to group replicates
        sample_series_2 (pd.Series): Series of the second sample. The index will be used to group replicates
        n_iterations (int): Number of bootstrap iterations to perform.

    Returns:
        mean_correlation (float): The mean Spearman correlation coefficient across all bootstrapped samples.
        confidence_interval (Tuple[float, float]): The 2.5th and 97.5th percentiles of the bootstrapped correlation coefficients.
        combined_p_value (float): Combined p-value from the bootstrap correlations using Pearson's method.

    """
    normal_dict_1 = _calc_norm_dict(sample_series_1)
    normal_dict_2 = _calc_norm_dict(sample_series_2)
    
    # Bootstrapping
    correlations = []
    pvalues = []
    for _ in range(n_iterations):
        a = [normal_dict_1[i].rvs(size=1)[0] for i in normal_dict_1]
        b = [normal_dict_2[i].rvs(size=1)[0] for i in normal_dict_2]
        corr, pvalue = spearmanr(a, b)
        correlations.append(corr)
        pvalues.append(pvalue)

    # Calculate mean correlation and confidence interval
    mean_correlation = np.mean(correlations)
    confidence_interval = np.percentile(correlations, [2.5, 97.5])

    # Combine p-values using Pearson's method
    _, combined_p_value = combine_pvalues(pvalues, method='pearson')

    return mean_correlation, confidence_interval, combined_p_value
