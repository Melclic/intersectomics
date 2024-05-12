from typing import List, Optional, Union
import pandas as pd
import plotly.express as px


def _melt_omics_layer(
    df: pd.DataFrame,
    col_replicate_name: str,
    biomolecules: List[str] = None,
    name: Optional[str] = None,
    normalize: Optional[bool] = True
) -> pd.DataFrame:
    """
    Transform a dataframe of omics data by melting it for visualization, optionally normalizing the data.

    Args:
        df (pd.DataFrame): Dataframe containing omics data with time as one of the column levels.
        col_replicate_name (str): The multiindex column level that contains the x-axis values and are replicates
        biomolecules (List[str], optional): List of biomolecules (e.g., proteins, transcripts, metabolites) to include in the transformation.
        name (Optional[str], optional): Name to assign to the 'Data Type' column in the resulting dataframe. Defaults to None.
        normalize (bool, optional): Whether to normalize the data. Defaults to True.

    Raises:
        KeyError: If any biomolecule in `biomolecules` is not found in `df`.
        KeyError: If `col_replicate_name` is not in the column names of `df`.

    Returns:
        pd.DataFrame: Melted dataframe ready for visualization, with optional normalization.
    """
    if col_replicate_name not in df.columns.names:
        raise KeyError('There must be {col_replicate_name} in the multiindex columns')
    if biomolecules:
        try:
            tmp_df = df.loc[biomolecules]
        except KeyError:
            raise KeyError(f'The biomolecule input does not exist in {name}')
    else:
        tmp_df = df
     
    if normalize:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
        tmp_df = pd.DataFrame(
            scaler.fit_transform(tmp_df.T), 
            index=tmp_df.columns,
            columns=tmp_df.index,
        )
        tmp_df = tmp_df.T 
    
    tmp_df = tmp_df.melt(ignore_index=False)
    tmp_df = tmp_df[[col_replicate_name, 'value']]
    if name:
        tmp_df['Data Type'] = [name] * tmp_df.shape[0]
    tmp_df = tmp_df.reset_index()
    tmp_df = tmp_df.rename(
            columns={
                'index': 'Biomolecule Name', 
                col_replicate_name: col_replicate_name.title()
                }
            )
    return tmp_df
    


def plot_time_single_omics_layer(
    df: pd.DataFrame,
    col_replicate_name: str,
    biomolecules: List[str] = None,
    num_col: int = 1,
    plot_width: int = 800,
    plot_height: int = 1000,
    trendline: str = "lowess",
    normalize: bool = True,
    title: str = None,
    y_title: str = None,
    x_title: str = None,
) -> px.scatter:
    """
    Plot time series data for multiple omics biomolecules from a single dataframe.
    The index are unique identifier of the biomolecule and the columns is a 
    multiindex with sample metadata.

    Args:
        df (pd.DataFrame): Dataframe containing time series data.
        col_replicate_name (str): The multiindex column level that contains the x-axis values and are replicates
        biomolecules (List[str]): List of biomolecules to plot.
        num_col (int, optional): Number of columns in the facet grid. Defaults to 1.
        plot_width (int, optional): Width of the plot. Defaults to 800.
        plot_height (int, optional): Height of the plot. Defaults to 1000.
        trendline (str, optional): Type of trendline to add to the plot. Defaults to "lowess".
        normalize (bool, optional): Whether to normalize the data. Defaults to True.
        title (str, optional): Title of the plot. Defaults to ''.
        title (Optional[str]): The title for the plot
        y_title (str, optional): Title for the y-axis. 
            Defaults to 'Normalized Values' if `normalize` is True or 'Values' if `normalize` is False.
        x_title (Optional[str]): The title of the x-axis. Defaults to `col_replicate_name`

    Raises:
        KeyError: If any biomolecule in `biomolecules` is not found in `df`.
        KeyError: If 'time' is not in the column names of `df`.

    Returns:
        px.scatter: A Plotly Express scatter plot with the specified parameters.
    """

    tmp_df = _melt_omics_layer(
            df=df,
            col_replicate_name=col_replicate_name,
            biomolecules=biomolecules,
            normalize=normalize)

    if not x_title:
        x_title = col_replicate_name.title()
    tmp_df = tmp_df.rename(columns={col_replicate_name: x_title})

    if not y_title:
        if normalize:
            y_title = 'Normalized Values'
        else:
            y_title = 'Values' 
    tmp_df = tmp_df.rename(columns={'value': y_title})

    # Plot
    fig = px.scatter(
        tmp_df, 
        x=x_title, 
        y=y_title, 
        facet_col='Biomolecule Name', 
        facet_col_wrap=num_col, 
        log_y=False,
        width=plot_width,
        height=plot_height,
        trendline=trendline,
        title=title,
    )
    fig.update_yaxes(matches=None)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
    
    return fig

def plot_time_multiple_omics_layers(
    dfs: List[pd.DataFrame],
    col_replicate_name: str,
    biomolecules: Optional[List[str]] = None,
    dfs_names: Optional[List[str]] = None,
    num_col: int = 1,
    plot_width: int = 800,
    plot_height: int = 1000,
    trendline: str = "lowess",
    normalize: bool = True,
    title: str = None,
    y_title: str = None,
    x_title: str = None,
) -> px.scatter:
    """
    Plot time series data for multiple biomolecules from multiple dataframes.

    Args:
        dfs (List[pd.DataFrame]): List of dataframes containing time series data. 
            The index are the genes/proteins/metabolites and the columns are the 
            samples with the column as metadata
        col_replicate_name (str): The multiindex column level that contains the x-axis values and are replicates
        biomolecules (Optional[List[str]], optional): List of biomolecules to plot. 
            Defaults to None, translates to all rows of dataframe.
        dfs_names (Optional[List[str]], optional): List of names corresponding to each dataframe. Defaults to None.
        num_col (int, optional): Number of columns in the facet grid. Defaults to 1.
        plot_width (int, optional): Width of the plot. Defaults to 800.
        plot_height (int, optional): Height of the plot. Defaults to 1000.
        trendline (str, optional): Type of trendline to add to the plot. Defaults to "lowess".
        normalize (bool, optional): Whether to normalize the data. Defaults to True.
        title (Optional[str]): The title for the plot
        y_title (str, optional): Title for the y-axis. 
            Defaults to 'Normalized Values' if `normalize` is True or 'Values' if `normalize` is False.
        x_title (Optional[str]): The title of the x-axis. Defaults to `col_replicate_name`

    Returns:
        px.scatter: A Plotly Express scatter plot with the specified parameters.
    """
    if dfs_names is None:
        dfs_names = [f'df_{i}' for i in range(len(dfs))]
    
    if len(dfs) != len(dfs_names):
        raise TypeError('The number of dataframes does not match the number of names')
    
    res_melted_dfs = []
    
    for df, name in zip(dfs, dfs_names):
        tmp_df = _melt_omics_layer(
                df,
                col_replicate_name=col_replicate_name,
                biomolecules=biomolecules,
                name=name,
                normalize=normalize)
        res_melted_dfs.append(tmp_df)
    
    # Merge all dataframes
    res_melted_dfs = pd.concat(res_melted_dfs)

    if not x_title:
        x_title = col_replicate_name.title()
    res_melted_dfs = res_melted_dfs.rename(columns={col_replicate_name: x_title})

    if not y_title:
        if normalize:
            y_title = 'Normalized Values'
        else:
            y_title = 'Values' 
    res_melted_dfs = res_melted_dfs.rename(columns={'value': y_title})
     
    # Plot
    fig = px.scatter(
        res_melted_dfs, 
        x=x_title, 
        y=y_title, 
        title=title,
        facet_col='Biomolecule Name', 
        facet_col_wrap=num_col, 
        color='Data Type',
        width=plot_width,
        height=plot_height,
        trendline=trendline,
    )
    fig.update_yaxes(matches=None)
    fig.for_each_yaxis(lambda yaxis: yaxis.update(showticklabels=True))
    
    return fig
