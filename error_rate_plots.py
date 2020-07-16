import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from constants import (
    BASES,
    CHANGES,
    LANES,
    PEOPLE,
    age_lane_map,
    base_changes_map,
    base_color_map,
    change_color_map,
)
from lookup import gene_df, gene_seqs_map, seq_df, seq_data_df

# for gene_number in gene_df.index:
for gene_number in [0]:
    print(gene_number)
    gene = gene_df.loc[gene_number, :]

    sequences = gene_seqs_map[gene_number]

    df = seq_data_df(sequences.index[0], group_by="position")
    for sequence in sequences.index[1:]:
        df = pd.concat(
            [df, seq_data_df(sequence, group_by="position")]
        ).drop_duplicates(keep=False)

    # positions = np.arange(gene["start"], gene["end"] + 1)

    plot_title = "Gene {} ".format(gene_number)
    for sequence in sequences.index:
        plot_title += seq_df.at[sequence, "strand"]

    # maps used for plotting
    base_fig_map = {}
    base_axs_map = {}
    change_error_rates_map = {}
    change_positions_map = {}

    for base in BASES:
        base_fig_map[base], base_axs_map[base] = plt.subplots(2, 2, figsize=(15, 8))
        for sub in base_changes_map[base]:
            print(sub)
            df_change = df.loc[df["sub"] == sub]
            error_rates = df_change["num subs"] / df_change["num consensus molecules"]
            change_error_rates_map[sub] = error_rates
            change_positions_map[sub] = df_change["position"]

    for base in BASES:
        axs = base_axs_map[base]
        axs = axs.flatten()
        for j, ax in enumerate(axs[:3]):
            for i, sub in enumerate(base_changes_map[base]):
                if i == j:
                    color = change_color_map[sub]
                    ax.set(title=sub)
                else:
                    color = "lightgrey"
                ax.plot(
                    change_positions_map[sub],
                    change_error_rates_map[sub],
                    label=sub,
                    linestyle="None",
                    marker="o",
                    color=color,
                    alpha=0.5,
                )
            ax.ticklabel_format(useOffset=False, style="plain")
            ax.set(xlabel="position", ylabel="error rate", yscale="log")
        for sub in base_changes_map[base]:
            axs[-1].plot(
                change_positions_map[sub],
                change_error_rates_map[sub],
                label=sub,
                linestyle="None",
                marker="o",
                color=change_color_map[sub],
                alpha=0.5,
            )
        axs[-1].legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
        axs[-1].set(title=base, xlabel="position", ylabel="error rate", yscale="log")

        for ax in axs:
            axtwin = ax.twinx()
            axtwin.plot(
                df["position"],
                df["num consensus molecules"],
                label="consensus",
                linestyle="None",
                marker="+",
                markersize="0.5",
                color="black",
                alpha=0.25,
            )
            axtwin.set(ylabel="number of consensus molecules")

        fig = base_fig_map[base]
        fig.suptitle("{} {}".format(plot_title, base), size=16, y=0.52)
        fig.subplots_adjust(top=0.8)
        fig.tight_layout()
        fig.savefig("plots\\{}_{}_error_rate.png".format(plot_title, base))
    plt.show()
    plt.close("all")
