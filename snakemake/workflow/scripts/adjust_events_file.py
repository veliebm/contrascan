"""
Adjust a BIDS event file to make its onsets start at a new time.

Created 3/2/22 by Benjamin Velie.
"""
import pandas


def main():
    """
    Entrypoint of the script.
    """
    events = pandas.read_csv(snakemake.input.events, sep="\t")
    events["onset"] += snakemake.params.adjust_by
    events.to_csv(snakemake.output.events, sep="\t", index=False)


if __name__ == "__main__":
    main()
