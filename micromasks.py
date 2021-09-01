#!/usr/bin/env python3
"""
Isolate a specific cluster from a collection of clusters.
"""
from os import PathLike
import pandas
import subprocess
import argparse


def main(clusters_image: PathLike, clusters_1d_file: PathLike, out_prefix: PathLike, get_negative_cluster: bool) -> None:
    """
    Isolate a specific cluster from a collection of clusters.
    """
    clusters_info = read_clusters_1D(clusters_1d_file)
    clusters_info["Strength"] = clusters_info["Volume"] * clusters_info["Mean"]
    cluster = _get_strongest_cluster(clusters_info, get_negative_cluster)

    args = f"""
    3dcalc
    -a {clusters_image}
    -expr equals(a,{cluster})
    -prefix {out_prefix}
    """.split()

    subprocess.run(args)


def _get_strongest_cluster(clusters_info: pandas.DataFrame, get_negative_cluster: bool) -> int:
    """
    Returns the index of the cluster we're interested in.
    """
    column = clusters_info[["Strength"]]
    index = column.idxmax()
    if get_negative_cluster == True:
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


if __name__ == "__main__":
    """
    This code only runs if this file is called as a script.
    """
    parser = argparse.ArgumentParser(description="Given an image clusterized by AFNI, turn the largest cluster in that image into a mask")

    parser.add_argument("--clusters_image", required=True, help="The clusterized image")
    parser.add_argument("--clusters_1d_file", required=True, help="The 1d file describing the clusters in the image")
    parser.add_argument("--out_prefix", required=True, help="Where to output the mask we make")
    parser.add_argument("--get_negative_cluster", action="store_true", help="If this is true, get the largest NEGATIVE cluster rather than the largest POSITIVE cluster")

    parsed_args = parser.parse_args()
    parsed_args_dict = vars(parsed_args)
    main(**parsed_args_dict)
