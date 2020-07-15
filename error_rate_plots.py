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
for gene_number in [244]:
    print(gene_number)
    change_error_rates_map = {}
    gene = gene_df.loc[gene_number, :]

    sequences = gene_seqs_map[gene_number].index

    df = seq_data_df(sequences[0])
    for sequence in sequences[1:]:
        df = pd.concat([df, seq_data_df(sequence)]).drop_duplicates(keep=False)

    positions = np.arange(gene["start"], gene["end"])

    plot_title = "Gene {} ".format(gene_number)
    for sequence in sequences:
        plot_title += seq_df.at[sequence, "strand"]

    base_fig_map = {}
    base_axs_map = {}

    for base in BASES:
        base_fig_map[base], base_axs_map[base] = plt.subplots(2, 2, figsize=(15, 8))
        for sub in base_changes_map[base]:
            print(sub)
            df_change = df.loc[df["sub"] == sub]
            # df.to_csv("data_files\\spam.csv")
            error_rates = np.zeros(positions.size)
            for i, position in enumerate(positions):
                df_position = df_change.loc[df_change["position"] == position]
                if len(df_position.index) != 0:
                    num_changes = np.sum(df_position["num subs"])
                    total_consensus = np.sum(df_position["num consensus molecules"])
                    error_rates[i] = num_changes / total_consensus

            change_error_rates_map[sub] = error_rates

    consensuses = np.zeros(positions.size)
    for i, position in enumerate(positions):
        df_position = df.loc[df["position"] == position]
        consensuses[i] = np.sum(df_position["num consensus molecules"])

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
                    positions,
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
                positions,
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
                positions, consensuses, label="consensus", color="black", alpha=0.25
            )
            axtwin.set(ylabel="number of consensus molecules")

        fig = base_fig_map[base]
        fig.suptitle("{} {}".format(plot_title, base), size=16, y=0.52)
        fig.subplots_adjust(top=0.8)
        fig.tight_layout()
        fig.savefig("plots\\{}_{}_error_rate.png".format(plot_title, base))
    plt.close("all")
