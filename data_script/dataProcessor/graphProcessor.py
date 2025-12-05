import pandas as pd
from matplotlib import rcParams
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import warnings

# Used to set some plotting options.
dashes = ["-", "--", "-.", ":"]
markers = ["o", "s", "^", "p", "h", "<", ">", "x", "D"]
colors = [
    "tab:blue",
    "tab:red",
    "tab:green",
    "tab:orange",
    "tab:pink",
    "tab:brown",
    "tab:olive",
    "tab:cyan",
    "tab:pink",
]


def createGraph(output, table, instances, algorithms, graphType, timeLimit):
    warnings.simplefilter(action="ignore", category=FutureWarning)

    xLabel = {"time": "Time (s)", "gap": "Root Gap (\\%)", "node": "Node Count"}

    # Set pandas dataframe.
    df = pd.DataFrame(
        columns=["Instance", "Algorithm", "Time (s)", "Root Gap (\\%)", "Node Count"]
    )

    i = 0
    n_instances = 0
    for instance in instances:
        n_instances += 1
        for alg in algorithms:
            if alg not in table[instance]:
                continue

            df.loc[i] = [
                instance,
                alg,
                table[instance][alg]["Time (s)"],
                table[instance][alg]["Gap"],
                table[instance][alg]["Node Count"],
            ]
            i += 1

    print("Considered instances: " + str(n_instances) + " / " + str(len(instances)))

    # Set some plot parameters.
    pgf_with_custom_preamble = {
        "font.family": "serif",  # use serif/main font for text elements
        "text.usetex": True,  # don't setup fonts from rc parameters
        "figure.figsize": (19, 15),
        "pgf.preamble": "\n".join(
            [
                r"\\usepackage{units}",  # load additional packages
                r"\\usepackage{metalogo}",
                r"\\usepackage{unicode-math}",  # unicode math setup
                r"\setmathfont{xits-math.otf}",
                r"\setmainfont{Arial}",  # serif font via preamble
            ]
        ),
    }

    rcParams.update(pgf_with_custom_preamble)
    plt.rc("xtick", labelsize=40)
    plt.rc("ytick", labelsize=40)

    ax = plt.subplot(1, 1, 1)

    # Create x ticks.
    xTicks = list(range(1, timeLimit))
    if graphType == "gap":
        xTicks = [x * 0.1 for x in range(0, 501)]
    if graphType == "node":
        xTicks = [x for x in range(1, 30000)]
        ax.set_xscale("log")

    # Plot lines.
    for alg in algorithms:
        ratioList = []

        for x in xTicks:
            auxDf = df.loc[
                (df["Algorithm"] == alg) & (df[xLabel[graphType]] <= x - 1e-6)
            ]
            auxDf2 = df.loc[
                (df["Algorithm"] == alg)
                & (df[xLabel[graphType]] >= (x - 1) + 1e-6)
                & (df[xLabel[graphType]] <= x - 1e-6)
            ]
            if graphType == "node":
                auxDf = auxDf.loc[(df["Time (s)"] <= timeLimit - 1e-6)]

            ratioList.append(float(auxDf.shape[0]) / n_instances)

        ax.step(
            xTicks,
            ratioList,
            dashes[algorithms.index(alg) % len(dashes)],
            marker=markers[algorithms.index(alg) % len(markers)],
            markersize=20,
            markevery=0.1,
            fillstyle="full",
            color=colors[algorithms.index(alg) % len(colors)],
            label="\\textsc{" + alg.lower() + "}",
            linewidth=3,
            where="pre",
        )

    # Set some other parameters and save figure.
    rcParams.update({"font.size": 50})
    ax.yaxis.grid(color="gray", linestyle="dashed")
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    plt.rcParams.update(pgf_with_custom_preamble)
    plt.rc("xtick", labelsize=40)
    plt.rc("ytick", labelsize=40)
    # ax.set_xscale("log")

    ax.legend(loc="best", ncol=1)
    plt.xlabel(xLabel[graphType], fontsize=50)
    plt.ylabel("Instances ratio", fontsize=50)
    plt.ylim((-0.1, 1.1))
    plt.savefig(output, bbox_inches="tight", pad_inches=0.2)


def createTimeAndGapGraph(output, table, instances, algorithms, timeLimit):
    warnings.simplefilter(action="ignore", category=FutureWarning)

    # set pandas dataframe
    df = pd.DataFrame(
        columns=[
            "Instance",
            "Algorithm",
            "Time (s)",
            "Final Gap",
            "Node Count",
        ]
    )

    i = 0
    n_instances = len(instances)
    for instance in instances:
        for alg in algorithms:
            df.loc[i] = [
                instance,
                alg,
                table[instance][alg]["Time (s)"],
                table[instance][alg]["Final Gap"],
                table[instance][alg]["Node Count"],
            ]
            i += 1

    pgf_with_custom_preamble = {
        "font.family": "serif",  # use serif/main font for text elements
        "text.usetex": True,  # don't setup fonts from rc parameters
        "figure.figsize": (35, 15),
        "pgf.preamble": "\n".join(
            [
                r"\\usepackage{units}",  # load additional packages
                r"\\usepackage{metalogo}",
                r"\\usepackage{unicode-math}",  # unicode math setup
                r"\setmathfont{xits-math.otf}",
                r"\setmainfont{Arial}",  # serif font via preamble
            ]
        ),
    }

    plt.rcParams.update(pgf_with_custom_preamble)
    plt.rc("xtick", labelsize=40)
    plt.rc("ytick", labelsize=40)
    timeTicks = list(range(1, timeLimit))
    gapTicks = [x * 0.1 for x in range(0, 151)]

    fig, axs = plt.subplots(1, 2, gridspec_kw={"width_ratios": [7, 3]})
    axs[0].set_ylim((-0.1, 1.1))

    axs[0].grid(axis="y", linestyle="--", alpha=0.7)
    axs[0].spines["right"].set_visible(False)

    axs[0].set_xlabel("Time (s)", fontsize=50)
    axs[0].set_ylabel("Instances Ratio", fontsize=50)
    axs[0].set_xlim(1, timeLimit)  # or max_time
    axs[0].grid(axis="x", linestyle="--", alpha=0.7)

    axs[1].set_ylim((-0.1, 1.1))

    axs[1].grid(axis="y", linestyle="--", alpha=0.7)
    # axs[1].yaxis.set_visible(False)
    axs[1].tick_params(axis="y", which="both", left=False, right=True, labelsize=False)

    axs[1].spines["left"].set_visible(False)

    axs[1].set_xlabel("Final Gap (\\%)", fontsize=50)
    axs[1].set_xlim(0.1, 15.0)
    axs[1].grid(axis="x", linestyle="--", alpha=0.7)

    for alg in algorithms:
        ratioListTime = []
        ratioListGap = []

        for x in timeTicks:
            auxDf = df.loc[(df["Algorithm"] == alg) & (df["Time (s)"] <= x)]
            ratioListTime.append(float(auxDf.shape[0]) / n_instances)

        for x in gapTicks:
            auxDf = df.loc[(df["Algorithm"] == alg) & (df["Final Gap"] <= x)]
            ratioListGap.append(float(auxDf.shape[0]) / n_instances)

        axs[0].step(
            timeTicks,
            ratioListTime,
            dashes[algorithms.index(alg) % len(dashes)],
            marker=markers[algorithms.index(alg) % len(markers)],
            markersize=20,
            markevery=0.1,
            fillstyle="full",
            color=colors[algorithms.index(alg) % len(colors)],
            label=alg,
            linewidth=3,
            where="pre",
        )

        axs[1].step(
            gapTicks,
            ratioListGap,
            dashes[algorithms.index(alg) % len(dashes)],
            marker=markers[algorithms.index(alg) % len(markers)],
            markersize=20,
            markevery=(0.1, 0.1),
            fillstyle="full",
            color=colors[algorithms.index(alg) % len(colors)],
            label=alg,
            linewidth=3,
            where="pre",
        )

    axs[0].axvline(
        x=timeLimit, color="black", linewidth=4, linestyle="--"
    )  # vertical line
    rcParams.update({"font.size": 50})
    axs[1].yaxis.grid(color="gray", linestyle="dashed")
    axs[1].yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    axs[1].xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    axs[1].legend(loc="lower right", ncol=1)
    # axs[1].legend(bbox_to_anchor=(3.00, 1), loc="upper right", ncol=1)
    plt.subplots_adjust(wspace=0)
    plt.savefig(output, bbox_inches="tight", pad_inches=0.2)
