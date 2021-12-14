


rule plot_occipital:
    input:
        underlay="../data/misc/kastner_cortex_masks/MNI152_T1_1mm_masked.nii.gz",
        overlay="../processed/correlation_whole_brain_ttest/startvolume-{startvolume}_data-values_correlations_alpha_ttest+tlrc.HEAD"
    output:
        plot="results/plots/startvolume-{startvolume}.png"
    params:
        threshold=0.1,
        coordinates=config["occipital coordinates"],
        overlay_subbrick=0
    conda:
        "../envs/neuroimaging.yaml"
    script:
        "../scripts/plot_fmri.py"
