import networkx as nx
import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Union


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


def make_corr_graph(
    corr: pd.DataFrame, 
    pvalues: pd.DataFrame,
    pvalue_cutoff: float = 0.05,
    corr_cutoff: float = 0.6,
) -> nx.Graph:
    """
    Create a correlation graph based on significant correlations.

    Args:
        corr (pd.DataFrame): DataFrame matrix containing correlation values for each pair of biomolecule.
        pvalues (pd.DataFrame): DataFrame matrix containing p-values corresponding to the correlations for each pair of biomolecule.
        pvalue_cutoff (float, optional): P-value threshold to consider a correlation significant. Defaults to 0.05.
        corr_cutoff (float, optional): Absolute correlation value threshold to include an edge in the graph. Defaults to 0.6.

    Returns:
        nx.Graph: A NetworkX graph where nodes represent biomolecules and edges represent significant correlations.

    Raises:
        None
    """
    signi_corrs = corr[pvalues <= pvalue_cutoff].melt(ignore_index=False).dropna().reset_index().rename(columns={'index': 'biomolecule1', 'variable': 'biomolecule2'})
    signi_corrs['abs_value'] = np.abs(signi_corrs['value'])
    signi_corrs = signi_corrs.sort_values(by='abs_value')
    signi_corrs = signi_corrs[signi_corrs['abs_value'] >= corr_cutoff]
    signi_corrs_dict = signi_corrs.T.to_dict()
    
    G = nx.Graph()
    G.add_nodes_from(np.unique(list(signi_corrs['biomolecule1']) + list(signi_corrs['biomolecule2'])))

    for i in signi_corrs_dict:
        # We only use the positive correlations, see README for a detailed
        # explanation why negative correlations are not useful when making 
        # graphs
        if signi_corrs_dict[i]['value'] >= 0:
            G.add_edge(
                signi_corrs_dict[i]['biomolecule1'],
                signi_corrs_dict[i]['biomolecule2'],
                weight=signi_corrs_dict[i]['abs_value'],
            )

    #remove unique
    G.remove_nodes_from(list(nx.isolates(G)))
    return G


def make_intersection_graph(
    corrs: List[pd.DataFrame],
    pvalues: List[pd.DataFrame],
    pvalue_cutoff: float = 0.05,
    corr_cutoff: float = 0.6,
) -> Tuple[nx.Graph, cm.ScalarMappable]:
    """
    Create an intersection graph from multiple correlation matrices and perform community detection.

    Args:
        corrs (List[pd.DataFrame]): List of dataframes containing correlation values.
        pvalues (List[pd.DataFrame]): List of dataframes containing p-values corresponding to the correlations.
        pvalue_cutoff (float, optional): P-value threshold to consider a correlation significant. Defaults to 0.05.
        corr_cutoff (float, optional): Absolute correlation value threshold to include an edge in the graph. Defaults to 0.6.

    Returns:
        Tuple[nx.Graph, dict]: A tuple containing the intersection graph and a colormap for community detection.
    """
    # Make the graphs
    graphs = []
    for corr, pvalue in zip(corrs, pvalues):
        graphs.append(
            make_corr_graph(
                corr=corr, 
                pvalues=pvalue,
                pvalue_cutoff=pvalue_cutoff,
                corr_cutoff=corr_cutoff,
            )
        )

    # Graph intersection
    G_inter = nx.intersection_all(graphs)
    G_inter = remove_G_pair_isolates(G_inter)

    # Perform community analysis
    communities = nx.community.louvain_communities(G_inter)
    comm_trans: Dict[int, List[str]] = {}
    trans_comm: Dict[str, int] = {}
    count = 0
    for community in communities:
        comm_trans[count] = list(community)
        for biomolecule in community:
            trans_comm[biomolecule] = count
        count += 1
    
    return G_inter, trans_comm

