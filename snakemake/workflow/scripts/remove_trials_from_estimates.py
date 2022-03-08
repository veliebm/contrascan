"""
Remove trials from our quadratic estimates.

Created 9/30/2022 By Benjamin Velie.
"""
import numpy


def main() -> None:
    """
    Entrypoint of the script.
    """
    # Load data.
    estimate = numpy.loadtxt(snakemake.input.estimate)
    trials = numpy.loadtxt(snakemake.input.volumes_containing_trials)

    # Trim trials to be same length as estimates.
    trimmed_trials = trials[:len(estimate)]
    
    # Remove trials from data.
    masked_estimate = numpy.ma.array(estimate, mask=trimmed_trials)
    non_trial_values = masked_estimate[~masked_estimate.mask]

    # Save data.
    numpy.savetxt(snakemake.output.trialless_estimate, non_trial_values)


if __name__ == "__main__":
    main()
