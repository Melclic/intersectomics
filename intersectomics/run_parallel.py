import pandas as pd
from tqdm import tqdm

from typing import Iterator, Dict, Any, List, Tuple
from itertools import combinations
from multiprocessing import Pool

from intersectomics.utils import combinations_count
from intersectomics.bootstrap_corr import bootstrap_spearman_correlation


def _pair_interator(
    combo_iterator: Iterator[Tuple[str, str]],
    dict_of_series: Dict[str, pd.Series],
    n_iterations: int = 10,
) -> Iterator[List[Tuple[str, str, pd.Series, pd.Series, int]]]:
    """Yields an iterator used for multiprocessing the spearman correlation

    Returns a tuple that is input to the _compute_pair_spearman_correlation fucntion

    Args:
        combo_iterator (Iterator[Tuple[str, str]]): An iterator that yields combinations of items, each a tuple of two elements.
        dict_of_series (Dict[str, pd.Series]): A dictionary mapping keys from the combinations to values that need to be processed alongside.
        n_iterations (int): A constant that represents the number of iterations, to be included in the output for each combination.

    Yields:
        Iterator[List[Tuple[str, str, any, any, int]]]: An iterator over chunks, where each chunk is a list of tuples. Each tuple
        contains the combination pair, the corresponding values from dict_of_series for these keys, and the number of iterations.
    """
    for combination in combo_iterator:
        yield (combination[0], combination[1], dict_of_series[combination[0]], dict_of_series[combination[1]], n_iterations)


def _compute_pair_spearman_correlation(
    args: Tuple[Tuple[str, str], pd.DataFrame, str, str, str, int]
) -> Tuple[str, str, float, float]:
    """Helper function that passes a multiprocessing tuple to the bootstrap_spearman_correlation
    function

    Calculates the correlation and p-values using Spearman correlation for a given pair of gene/protein/metabolite 
    given an expression table. Takes a tuple of arguments necessary for computing the Spearman correlation
    between two specified gene/protei/metabolite over a series of time points within a provided DataFrame.

    Args:
        args (Tuple[Tuple[str, str], pd.DataFrame, str, str, str, int]):
            - s1 (str): name of the first gene/protein/metabolite
            - s2 (str): name of the second gene/protein/metabolite
            - series_1 (pd.Series): First Series with expression levels as values and replicate as index
            - series_2 (pd.Series): Second Series with expression levels as values and replicate as index
            - n_iterations (int): Number of bootstrap iterations to perform.

    Returns:
        Tuple[str, str, float, float]:
            - Name of the first gene/protein/metabolite
            - Name of the second gene/protein/metabolite
            - Spearman correlation coefficient between the two genes.
            - Combined p-value of the correlation.

    """
    s1, s2, series_1, series_2, n_iterations = args

    corr, c_i, combined_pvalue = bootstrap_spearman_correlation(
        series_1,
        series_2,
        n_iterations=n_iterations
    )
    return s1, s2, corr, combined_pvalue


def bootstrap_spearman_corr_parallel(
    df_input: pd.DataFrame,
    replicate_column_name: str,
    n_processes: int = 4,
    n_iterations: int = 10,
    chunk_size: int = 10,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]]]:
    """Performs bootstrap Spearman correlation analysis for all unique gene pairs in a dataframe, using parallel processing.

    This function uses multiprocessing to parallelize the bootstrap Spearman 
    correlation analysis for pairs of samples with replicates across multiple 
    processes. It computes Spearman correlation coefficients and combined p-values 
    for each pair, returning results in dictionaries.

    Args:
        df_input (pd.DataFrame): The dataframe containing omics data in a wide format with 
            the rows as genes/proteins/metabolites and the columns as samples.
            The column may be a MultiIndex.
        replicate_column_name (str): Column name of replicates that will be grouped together. 
        n_processes (int): Number of processes to use for parallel computation.
        n_iterations (int): Number of bootstrap iterations to perform for each gene pair.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]]:
            - DataFrame of correlation coefficients for each pair.
            - DataFrame of combined p-values for the correlations of each pair.

    Raises:
        IndexError: If any specified column names are not present in the dataframe.
    """
    if not replicate_column_name in df_input.columns.names:
        raise IndexError('Could not find one of the specified columns in the dataframe. Valid options are: ' + str(df_input.columns.names))
        
    # Make sure that the grouped axis are sorted and that the multiindex is sinlge
    df = df_input.copy(deep=True)
    df = df.sort_index(axis=1, level=replicate_column_name)
    df.columns = df.columns.get_level_values(replicate_column_name)

    iter_pairs = combinations(df.index.unique(), 2)
    iter_length = combinations_count(len(df.index.unique()), 2)
    dict_of_series = {index: row for index, row in df.iterrows()}

    correlations = {}
    pvalues = {}

    with Pool(processes=n_processes) as pool:
        results = tqdm(
            pool.imap_unordered(
                _compute_pair_spearman_correlation, 
                _pair_interator(
                    combo_iterator=iter_pairs, 
                    dict_of_series=dict_of_series, 
                    n_iterations=n_iterations, 
                ),
                chunksize=chunk_size,
            ),
            total=iter_length,
        )

        for s1, s2, corr, combined_pvalue in results:
            if s1 not in correlations:
                correlations[s1] = {}
            correlations[s1][s2] = corr

            if s1 not in pvalues:
                pvalues[s1] = {}
            pvalues[s1][s2] = combined_pvalue

    return pd.DataFrame(correlations), pd.DataFrame(pvalues)
