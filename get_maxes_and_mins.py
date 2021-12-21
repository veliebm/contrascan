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
    permutation_arrays = [image.dataobj[:, :, :, 0] for image in permutation_images]

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
    kwargs = {'in_permutations': ['./processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_1+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_2+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_3+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_4+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_5+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_6+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_7+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_8+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_9+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_10+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_11+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_12+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_13+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_14+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_15+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_16+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_17+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_18+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_19+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_20+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_21+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_22+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_23+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_24+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_25+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_26+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_27+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_28+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_29+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_30+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_31+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_32+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_33+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_34+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_35+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_36+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_37+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_38+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_39+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_40+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_41+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_42+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_43+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_44+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_45+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_46+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_47+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_48+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_49+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_50+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_51+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_52+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_53+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_54+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_55+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_56+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_57+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_58+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_59+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_60+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_61+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_62+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_63+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_64+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_65+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_66+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_67+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_68+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_69+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_70+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_71+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_72+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_73+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_74+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_75+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_76+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_77+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_78+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_79+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_80+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_81+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_82+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_83+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_84+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_85+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_86+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_87+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_88+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_89+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_90+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_91+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_92+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_93+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_94+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_95+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_96+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_97+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_98+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_99+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_100+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_101+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_102+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_103+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_104+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_105+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_106+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_107+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_108+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_109+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_110+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_111+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_112+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_113+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_114+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_115+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_116+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_117+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_118+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_119+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_120+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_121+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_122+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_123+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_124+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_125+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_126+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_127+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_128+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_129+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_130+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_131+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_132+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_133+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_134+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_135+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_136+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_137+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_138+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_139+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_140+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_141+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_142+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_143+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_144+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_145+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_146+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_147+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_148+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_149+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_150+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_151+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_152+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_153+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_154+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_155+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_156+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_157+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_158+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_159+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_160+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_161+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_162+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_163+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_164+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_165+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_166+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_167+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_168+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_169+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_170+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_171+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_172+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_173+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_174+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_175+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_176+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_177+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_178+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_179+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_180+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_181+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_182+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_183+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_184+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_185+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_186+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_187+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_188+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_189+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_190+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_191+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_192+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_193+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_194+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_195+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_196+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_197+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_198+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_199+tlrc.HEAD', './processed/correlation_whole_brain_ttest/permutation_testing/alpha/startvolume-4_variable-amplitude_alpha_200+tlrc.HEAD'], 'out_maxes': './processed/maxes_mins_quantiles/startvolume-4_variable-amplitude_alpha_maxes.csv', 'out_mins': './processed/maxes_mins_quantiles/startvolume-4_variable-amplitude_alpha_mins.csv'}
    main(**kwargs)


def make_parent_dir(path: PathLike) -> None:
    """
    Create the parent directory for a path.
    """
    Path(path).parent.mkdir(exist_ok=True, parents=True)