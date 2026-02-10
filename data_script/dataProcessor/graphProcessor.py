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


def createPaper1Graph(
    ccvrp_table,
    jabali_table,
    ccvrp_instances,
    jabali_instances,
    algorithms,
    timelimit,
    output,
):
    warnings.simplefilter(action="ignore", category=FutureWarning)

    xLabel = {"time": "Time (s)", "gap": "Root Gap (\\%)", "node": "Node Count"}

    # Set pandas dataframes.
    df = pd.DataFrame(
        columns=[
            "InstanceType",
            "Instance",
            "Algorithm",
            "Time (s)",
            "Root Gap (\\%)",
            "Node Count",
        ]
    )

    i = 0
    for instanceType in ["CCVRP", "Jabali"]:
        instances = ccvrp_instances if instanceType == "CCVRP" else jabali_instances
        table = ccvrp_table if instanceType == "CCVRP" else jabali_table

        for instance in instances:
            if instance not in table:
                continue

            for alg in algorithms:
                if alg not in table[instance]:
                    continue

                df.loc[i] = [
                    instanceType,
                    instance,
                    alg,
                    table[instance][alg]["Time (s)"],
                    table[instance][alg]["Gap"],
                    table[instance][alg]["Node Count"],
                ]
                i += 1

    # Number of instances for each type
    n_ccvrp = len(ccvrp_instances)
    n_jabali = len(jabali_instances)

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

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(38, 15))

    # Define x ticks for both graph types
    graphTypes = ["time", "gap"]
    axes = [ax1, ax2]
    xTicksDict = {
        "time": list(range(1, timelimit)),
        "gap": [x * 0.1 for x in range(0, 501)],
    }

    # Plot lines for each graph type
    for graphType, ax in zip(graphTypes, axes):
        xTicks = xTicksDict[graphType]

        for alg in algorithms:
            ratioList = []
            for x in xTicks:
                # Count instances from CCVRP
                ccvrp_count = df.loc[
                    (df["Algorithm"] == alg)
                    & (df["InstanceType"] == "CCVRP")
                    & (df[xLabel[graphType]] <= x - 1e-6)
                ].shape[0]

                # Count instances from Jabali
                jabali_count = df.loc[
                    (df["Algorithm"] == alg)
                    & (df["InstanceType"] == "Jabali")
                    & (df[xLabel[graphType]] <= x - 1e-6)
                ].shape[0]

                # Compute weighted average.
                ccvrp_ratio = ccvrp_count / n_ccvrp if n_ccvrp > 0 else 0
                jabali_ratio = jabali_count / n_jabali if n_jabali > 0 else 0
                weighted_ratio = 0.5 * ccvrp_ratio + 0.5 * jabali_ratio

                ratioList.append(weighted_ratio)

            ax.step(
                xTicks,
                ratioList,
                dashes[algorithms.index(alg) % len(dashes)],
                marker=markers[algorithms.index(alg) % len(markers)],
                markersize=30,
                markevery=0.1,
                fillstyle="full",
                color=colors[algorithms.index(alg) % len(colors)],
                label="\\textsc{" + alg.lower() + "}",
                linewidth=5,
                where="pre",
            )

        # Set parameters for each subplot
        ax.yaxis.grid(color="gray", linestyle="dashed")
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.set_xlabel(xLabel[graphType], fontsize=50)
        ax.set_ylabel("Instances ratio", fontsize=50)
        ax.set_ylim((-0.1, 1.1))

    # Create a single legend below the plots
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(algorithms),
        fontsize=50,
        bbox_to_anchor=(0.5, -0.02),
    )

    # Adjust layout to make room for the legend below
    plt.tight_layout(rect=[0, 0.08, 1, 1], w_pad=3.0)
    plt.subplots_adjust(wspace=0.2)

    plt.savefig(output, bbox_inches="tight", pad_inches=0.2)
    plt.close()


def createPaper2Graph(
    ccvrp_table,
    ccvrp_table_basic,
    jabali_table,
    jabali_table_basic,
    ccvrp_instances,
    jabali_instances,
    algorithms,
    timelimit,
    output,
):
    warnings.simplefilter(action="ignore", category=FutureWarning)

    xLabel = {"time": "Time (s)", "gap": "Root Gap (\\%)", "node": "Node Count"}

    # Set pandas dataframes.
    df = pd.DataFrame(
        columns=[
            "InstanceType",
            "Instance",
            "Algorithm",
            "Time (s)",
            "Root Gap (\\%)",
            "Node Count",
            "Basic",
        ]
    )

    i = 0
    for instanceType, isBasic in [
        ("CCVRP", False),
        ("CCVRP", True),
        ("Jabali", False),
        ("Jabali", True),
    ]:
        instances = ccvrp_instances if instanceType == "CCVRP" else jabali_instances

        if instanceType == "CCVRP" and not isBasic:
            table = ccvrp_table
        elif instanceType == "CCVRP" and isBasic:
            table = ccvrp_table_basic
        elif instanceType == "Jabali" and not isBasic:
            table = jabali_table
        else:
            table = jabali_table_basic

        for instance in instances:
            if instance not in table:
                continue

            for alg in algorithms:
                if alg not in table[instance]:
                    continue

                df.loc[i] = [
                    instanceType,
                    instance,
                    alg,
                    table[instance][alg]["Time (s)"],
                    table[instance][alg]["Gap"],
                    table[instance][alg]["Node Count"],
                    isBasic,
                ]
                i += 1

    # Number of instances for each type
    n_ccvrp = len(ccvrp_instances)
    n_jabali = len(jabali_instances)

    # Set some plot parameters.
    pgf_with_custom_preamble = {
        "font.family": "serif",
        "text.usetex": True,
        "figure.figsize": (38, 15),
        "pgf.preamble": "\n".join(
            [
                r"\\usepackage{units}",
                r"\\usepackage{metalogo}",
                r"\\usepackage{unicode-math}",
                r"\setmathfont{xits-math.otf}",
                r"\setmainfont{Arial}",
            ]
        ),
    }

    rcParams.update(pgf_with_custom_preamble)
    plt.rc("xtick", labelsize=40)
    plt.rc("ytick", labelsize=40)

    # Create figure with two subplots (Normal vs Basic)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(38, 15))

    # Define x ticks for time
    xTicks = list(range(1, timelimit))

    # Plot for both Normal (left) and Basic (right) instances
    for ax, isBasic, title in zip(
        [ax1, ax2],
        [False, True],
        [
            r"$\mathcal{X}_{\textsc{cvrp}}$",
            r"$\mathcal{X}_{\textsc{sub}}$",
        ],
    ):
        for alg in algorithms:
            ratioList = []
            for x in xTicks:
                # Count instances from CCVRP
                ccvrp_count = df.loc[
                    (df["Algorithm"] == alg)
                    & (df["InstanceType"] == "CCVRP")
                    & (df["Basic"] == isBasic)
                    & (df[xLabel["time"]] <= x - 1e-6)
                ].shape[0]

                # Count instances from Jabali
                jabali_count = df.loc[
                    (df["Algorithm"] == alg)
                    & (df["InstanceType"] == "Jabali")
                    & (df["Basic"] == isBasic)
                    & (df[xLabel["time"]] <= x - 1e-6)
                ].shape[0]

                # Compute weighted average
                ccvrp_ratio = ccvrp_count / n_ccvrp if n_ccvrp > 0 else 0
                jabali_ratio = jabali_count / n_jabali if n_jabali > 0 else 0
                weighted_ratio = 0.5 * ccvrp_ratio + 0.5 * jabali_ratio

                ratioList.append(weighted_ratio)

            ax.step(
                xTicks,
                ratioList,
                dashes[algorithms.index(alg) % len(dashes)],
                marker=markers[algorithms.index(alg) % len(markers)],
                markersize=30,
                markevery=0.1,
                fillstyle="full",
                color=colors[algorithms.index(alg) % len(colors)],
                label="\\textsc{" + alg.lower() + "}",
                linewidth=5,
                where="pre",
            )

        # Set parameters for each subplot
        ax.yaxis.grid(color="gray", linestyle="dashed")
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.set_xlabel(xLabel["time"], fontsize=50)
        ax.set_ylabel("Instances ratio", fontsize=50)
        ax.set_ylim((-0.1, 1.1))
        ax.set_title(title, fontsize=50, pad=25)

    # Create a single legend below the plots
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=len(algorithms),
        fontsize=50,
        bbox_to_anchor=(0.5, -0.02),
    )

    # Adjust layout to make room for the legend below
    plt.tight_layout(rect=[0, 0.08, 1, 1], w_pad=3.0)
    plt.subplots_adjust(wspace=0.2)

    plt.savefig(output, bbox_inches="tight", pad_inches=0.2)
    plt.close()
