"""
ttest the results of using 3dmaskave to get the average correlation within an ROI.

Created 10/18/2021 by Benjamin Velie
veliebm@gmail.com
"""
import scipy.stats
import pandas
from os import PathLike
from typing import List
from pathlib import Path


def main(from_text_files: List[PathLike], to_table: PathLike, to_statistics: PathLike) -> None:
    """
    ttest the results of using 3dmaskave to get the average correlation within an ROI.
    """
    # Make out dirs.
    for path in (to_table, to_statistics):
        Path(path).parent.mkdir(exist_ok=True, parents=True)

    # Extract correlation average from each text file.
    correlations = {source_file: extract_correlation(source_file) for source_file in from_text_files}
    
    # Concatenate the correlations into a summary pandas table.
    labelled_correlations = dict()
    labelled_correlations["filename"] = correlations.keys()
    labelled_correlations["correlation"] = correlations.values()
    correlations_table = pandas.DataFrame(labelled_correlations)

    # Calculate t-test.
    results = scipy.stats.ttest_1samp(a=correlations_table["correlation"], popmean=0, alternative='two-sided')

    # Save t-test results.
    results_dict = {"pvalue": results.pvalue, "statistic": results.statistic}
    results_table = pandas.DataFrame(results_dict, index=[0])
    results_table.to_csv(to_statistics)

    # Save summary table.
    correlations_table.to_csv(to_table)


def extract_correlation(path_to_file: PathLike) -> float:
    """
    Extract a correlation from a text file.

    Reads the first line of the text file and uses whatever number is finds there.
    """
    with open(path_to_file, "r") as io:
        lines = io.readlines()
        average = float(lines[0].split()[0])
        return average


def _test_module() -> None:
    """
    Test this module to make sure it works.
    """
    main(
        **{'from_text_files': ['./processed/correlation_averages/sub-104_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-106_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-107_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-108_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-109_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-110_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-111_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-112_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-113_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-115_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-116_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-117_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-120_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-121_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-122_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-123_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-124_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt', './processed/correlation_averages/sub-125_startvolume-1_source-occipital_data-SNRs_correlation_average_alpha.txt'], 'to_table': './processed/correlation_averages_ttests/variable-SNRs_startvolume-1_mask-occipital_analysis-alpha_outfile-table.csv', 'to_statistics': './processed/correlation_averages_ttests/variable-SNRs_startvolume-1_mask-occipital_analysis-alpha_outfile-results.csv'}
    )
