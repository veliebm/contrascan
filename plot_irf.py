"""
Make a nice plot of a representative occipital IRF.

Created 1/6/22 by Benjamin Velie.
"""

import matplotlib.pyplot as plt


def main():
    """
    Entrypoint of program.
    """
    calcarine_timeseries = [0, -.07243, -.088155, -.171879, -.143361, .037872, .066945, .032208, .013069, 0]
    occipital_timeseries = [0, -.051867, .238449, .579762, .771096, .711357, .295275, -.034379, 0.182821, 0]
    time_vector = [time*2 for time in range(len(calcarine_timeseries))]

    plot_irf(x_values=time_vector, calcarine_values=calcarine_timeseries, occipital_values=occipital_timeseries)


def plot_irf(x_values, calcarine_values, occipital_values) -> None:
    """
    Plot the IRFs
    """
    plt.rcParams.update({'font.size': 40})
    plt.rcParams.update({'axes.linewidth': 5})
    fig, ax = plt.subplots()
    plt.axhline(y=0, color='black', linestyle='--', linewidth=5)
    ax.plot(x_values, calcarine_values, linewidth=10, color='#377eb8')
    ax.plot(x_values, calcarine_values, color="#377eb8", markevery=[4], marker="o", ls="", label="points", ms=20)
    ax.plot(x_values, occipital_values, linewidth=10, color='#ff7f00')
    ax.plot(x_values, occipital_values, color="#ff7f00", markevery=[4], marker="o", ls="", label="points", ms=20)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('BOLD Δ%')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(width=5, length=10)
    fig.set_size_inches(10, 9)
    plt.xticks(x_values[::2])
    plt.yticks([-.2, 0, .8])
    plt.tight_layout()
    plt.savefig('processed/plots/representative_IRFs.png')


if __name__ == "__main__":
    main()
