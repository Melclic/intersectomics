from math import factorial
import pandas as pd
import networkx as nx
from typing import Union


def remove_G_pair_isolates(G_in: Union[nx.Graph, nx.DiGraph]) -> Union[nx.Graph, nx.DiGraph]:
    """
    Remove nodes that are connected only to one another and isolate nodes from a graph.

    Args:
        G_in (Union[nx.Graph, nx.DiGraph]): Input graph from which nodes will be removed.

    Returns:
        Union[nx.Graph, nx.DiGraph]: A copy of the input graph with specified nodes removed.

    Raises:
        None
    """
    G = G_in.copy()
    to_remove = []

    for node in G.nodes():
        if node not in to_remove:
            neighbors = list(G.neighbors(node))
            if len(neighbors) == 0:
                to_remove.append(node)
            elif len(neighbors) == 1:
                neighbors_of_neighbor = list(G.neighbors(neighbors[0]))
                if len(neighbors_of_neighbor) == 1 and neighbors_of_neighbor[0] == node:
                    to_remove.append(node)
                    to_remove.append(neighbors[0])

    G.remove_nodes_from(to_remove)
    G.remove_nodes_from(list(nx.isolates(G)))  # Remove isolated nodes

    return G


def combinations_count(n: int, r: int) -> int:
    """
    Calculate the number of combinations for choosing r elements from a set of n elements.

    Args:
        n (int): Total number of elements.
        r (int): Number of elements to choose.

    Returns:
        int: Number of combinations.
    """
    if r > n:
        return 0  # If the number of chosen elements is greater than available, return 0
    return factorial(n) // (factorial(r) * factorial(n - r))


def add_cols_multi_index(
    df_data: pd.DataFrame, 
    df_metadata: pd.DataFrame, 
    metadata_index_column: str
) -> pd.DataFrame:
    """
    Integrates metadata into the DataFrame's columns as a MultiIndex. This function assumes
    that the DataFrame's columns should match the entries in a specified column of the metadata DataFrame.

    Args:
        df_data (pd.DataFrame): The original DataFrame whose columns are to be enriched with metadata.
        df_metadata (pd.DataFrame): The metadata DataFrame containing additional information that will be 
                                    used to form the MultiIndex.
        metadata_index_column (str): The column name in `df_metadata` that corresponds to the columns of `df_data`.
                                    This column must perfectly match the columns of `df_data` in sorted order.

    Returns:
        pd.DataFrame: A new DataFrame with columns transformed into a MultiIndex using metadata from `df_metadata`.

    Raises:
        IndexError: If `metadata_index_column` is not in `df_metadata` or if the columns of `df_data` do not match
                    the sorted values of `metadata_index_column` in `df_metadata`.
    """
    if metadata_index_column not in df_metadata:
        raise IndexError(f'The {metadata_index_column} is not in the columns of the metadata')

    #TODO: check if the metadata matches the shape of input dataframe

    tmp_df = df_data.copy(deep=True)
    tmp_metdata_df = df_metadata.copy(deep=True)
    
    tmp_df = tmp_df[sorted(tmp_df.columns)]
    tmp_metdata_df = tmp_metdata_df.sort_values(by=metadata_index_column)

    if tmp_df.columns.to_list() != tmp_metdata_df[metadata_index_column].to_list():
        raise IndexError('The columns index and the metadata_index_column do not match')

    tmp_df.columns = pd.MultiIndex.from_frame(tmp_metdata_df)
    return tmp_df
