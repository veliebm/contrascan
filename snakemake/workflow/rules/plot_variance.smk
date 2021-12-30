"""
Plot mean variance images.
"""

rule strip_skull_of_variance:
    """
    Mask out the skull from our variance images.
    """
    input:
        mask="../data/misc/kastner_cortex_masks/MNI152_T1_2.5mm_full_mask.nii.gz",
        image="../processed/variance_whole_brain/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz"
    output:
        image="results/plot_variance/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.nii.gz",
    log:
        "logs/plot_variance/skull_stripped/startvolume-{startvolume}_variable-{variable}_baselined-{baselined}_{analysis}.log"
    conda:
        "../envs/afni.yaml"
    shell:
        "3dcalc -float -a {input.image} -b {input.mask} -expr 'a*step(b)' -prefix {output.image} 2> {log}"
