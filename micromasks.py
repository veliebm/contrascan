#!/usr/bin/env python3
"""
Isolate a specific cluster from a collection of clusters.
"""
from os import PathLike
import pandas
import subprocess


def main(clusters_image: PathLike, clusters_1d_file: PathLike, out_prefix: PathLike, get_positive_cluster: bool=True) -> None:
    """
    Isolate a specific cluster from a collection of clusters.
    """
    clusters_info = read_clusters_1D(clusters_1d_file)
    clusters_info["Strength"] = clusters_info["Volume"] * clusters_info["Mean"]
    cluster = _get_strongest_cluster(clusters_info, get_positive_cluster)

    args = f"""
    3dcalc
    -a {clusters_image}
    -expr equals(a,{cluster})
    -prefix {out_prefix}
    """.split()

    subprocess.run(args)


def _get_strongest_cluster(clusters_info: pandas.DataFrame, get_positive_cluster: bool) -> int:
    """
    Returns the index of the cluster we're interested in.
    """
    column = clusters_info[["Strength"]]
    index = column.idxmax()
    if get_positive_cluster == False:
        index = column.idxmin()
    
    return int(index)


def read_clusters_1D(in_path: PathLike) -> pandas.DataFrame:
    """
    Read specifically the clusters 1D file. Gives you headers!!! Adjusts your index!!!!!! Yay!
    """
    header = "Volume CM_RL CM_AP CM_IS minRL maxRL minAP maxAP minIS maxIS Mean SEM Max_Int MI_RL MI_AP MI_IS".split()
    table = read_1D(in_path)
    table.columns=header
    table.index += 1

    return table


def read_1D(in_path: PathLike) -> pandas.DataFrame:
    """
    Read an AFNI 1D file.
    """
    return pandas.read_table(in_path, sep=" +", header=None, index_col=False, comment="#", engine="python")
