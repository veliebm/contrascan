"""
Get distributions of our maxes and mins for our permutations.

Created 12/8/2021 by Benjamin Velie.
veliebm@gmail.com
"""
from os import PathLike
from typing import List
import nibabel
import numpy
import pandas
from pathlib import Path


def main(in_permutations: List[PathLike], out_maxes: PathLike, out_mins: PathLike) -> None:
    """
    Get distributions of our maxes and mins for our permutations.

    Args:
        in_permutations (List[PathLike]): List of paths to permutation images.
        out_maxes (PathLike): Table of max values from those images.
        out_mins (PathLike): Table of min values from those images.
    """
    # Read images.
    permutation_images = [nibabel.load(path) for path in in_permutations]

    # Get data of images as list of numpy arrays.
    try:
        permutation_arrays = [image.dataobj[:, :, :, 0] for image in permutation_images]
    except IndexError:
        print("At least one image isn't 4d. Trying to read as 3d images.")
        permutation_arrays = [image.dataobj[:, :, :] for image in permutation_images]

    # Combine list of arrays into one big array. Actual data located in index 0 of array.
    all_permutations_array = numpy.stack(permutation_arrays, 3)

    # Flatten array into 2 dimensions: 1 is the index of the image, 0 is the data of that image.
    num_voxels_in_each_image = all_permutations_array.shape[0] * all_permutations_array.shape[1] * all_permutations_array.shape[2]
    array_length = all_permutations_array.shape[3]
    flattened_array = numpy.reshape(all_permutations_array, (num_voxels_in_each_image, array_length))

    # Calculate maxes and mins for each image in big array.
    maxes_array, mins_array = flattened_array.max(axis=0), flattened_array.min(axis=0)

    # Store arrays into a DataFrame.
    maxes_dataframe, mins_dataframe = (pandas.DataFrame(data=quantile_array.T, index=in_permutations, columns=["value"]) for quantile_array in (maxes_array, mins_array))

    for quantile_dataframe in (maxes_dataframe, mins_dataframe):

        # Rank data.
        num_values = len(quantile_dataframe["value"])
        quantile_dataframe["percentile"] = (quantile_dataframe["value"].rank() - 1) / num_values * 100

        # Sort data.
        quantile_dataframe.sort_values("percentile", ascending=False, inplace=True)

    # Save results to disk.
    make_parent_dir(out_maxes)
    maxes_dataframe.to_csv(out_maxes)
    mins_dataframe.to_csv(out_mins)


def _test_module():
    """
    Test this module.
    """
    kwargs = {'out_maxes': './processed/maxes_and_mins_distributions/startvolume-5_variable-amplitude_baselined-false_ssvep_maxes_table.csv', 'out_mins': './processed/maxes_and_mins_distributions/startvolume-5_variable-amplitude_baselined-false_ssvep_mins_table.csv', 'in_permutations': ['./processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_1+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_2+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_3+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_4+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_5+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_6+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_7+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_8+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_9+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_10+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_11+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_12+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_13+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_14+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_15+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_16+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_17+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_18+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_19+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_20+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_21+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_22+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_23+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_24+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_25+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_26+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_27+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_28+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_29+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_30+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_31+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_32+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_33+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_34+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_35+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_36+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_37+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_38+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_39+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_40+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_41+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_42+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_43+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_44+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_45+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_46+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_47+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_48+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_49+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_50+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_51+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_52+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_53+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_54+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_55+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_56+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_57+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_58+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_59+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_60+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_61+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_62+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_63+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_64+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_65+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_66+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_67+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_68+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_69+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_70+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_71+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_72+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_73+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_74+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_75+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_76+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_77+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_78+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_79+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_80+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_81+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_82+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_83+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_84+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_85+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_86+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_87+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_88+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_89+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_90+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_91+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_92+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_93+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_94+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_95+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_96+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_97+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_98+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_99+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_100+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_101+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_102+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_103+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_104+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_105+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_106+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_107+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_108+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_109+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_110+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_111+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_112+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_113+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_114+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_115+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_116+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_117+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_118+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_119+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_120+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_121+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_122+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_123+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_124+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_125+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_126+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_127+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_128+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_129+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_130+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_131+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_132+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_133+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_134+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_135+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_136+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_137+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_138+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_139+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_140+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_141+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_142+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_143+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_144+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_145+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_146+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_147+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_148+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_149+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_150+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_151+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_152+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_153+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_154+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_155+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_156+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_157+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_158+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_159+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_160+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_161+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_162+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_163+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_164+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_165+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_166+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_167+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_168+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_169+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_170+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_171+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_172+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_173+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_174+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_175+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_176+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_177+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_178+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_179+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_180+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_181+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_182+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_183+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_184+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_185+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_186+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_187+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_188+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_189+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_190+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_191+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_192+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_193+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_194+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_195+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_196+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_197+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_198+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_199+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/ssvep/startvolume-5_variable-amplitude_ssvep_200+tlrc.HEAD']}
    main(**kwargs)


def make_parent_dir(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)
