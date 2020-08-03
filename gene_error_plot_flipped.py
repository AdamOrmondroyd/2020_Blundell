"""
# 2020 Blundell lab internship

Generates plots of error rates for given genes using the flipped data
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from constants import BASES, base_subs_map, sub_color_map
from lookup import gene_df, gene_seqs_map, seq_data_df


def gene_error_plot_flipped(gene_number, downsample=True):
    """
    Saves plots of errors for a given gene, separating + and - strand data.

    downsample = True will plot the downsampled number of errors.
    downsample = False will plot the error rate for the full number of errors.
    """

    print(gene_number)
    gene = gene_df.loc[gene_number, :]

    plot_title = "Gene {} flipped ".format(gene_number)

    seq_df = gene_seqs_map[gene_number]

    df = pd.DataFrame()
    for i in seq_df.index:
        df = pd.concat([df, seq_data_df(i, group_by="position")]).drop_duplicates(
            keep="first"
        )

    # Make up plot title
    for strand in seq_df["strand"]:
        plot_title += strand

    for base in BASES:
        fig, axs = plt.subplots(2, 2, figsize=(15, 8))

        axs = axs.flatten()
        for j, ax in enumerate(axs):
            for i, sub in enumerate(base_subs_map[base]):
                if (3 == j) or (i == j):
                    color = sub_color_map[sub]
                    if 3 == j:
                        ax.set(title=base)
                    else:
                        ax.set(title=sub)
                else:
                    color = "lightgrey"

                sub_df = df.loc[df["sub"] == sub]

                if downsample:

                    ax.plot(
                        sub_df["position"],
                        sub_df["downsample"],
                        label=sub + "+",
                        linestyle="None",
                        marker="^",
                        color=color,
                        alpha=0.5,
                    )
                else:
                    ax.plot(
                        sub_df["position"],
                        sub_df["num subs"] / sub_df["num consensus molecules"],
                        label=sub + "+",
                        linestyle="None",
                        marker="^",
                        color=color,
                        alpha=0.5,
                    )

            ax.ticklabel_format(useOffset=False, style="plain")

            if downsample:
                ax.set(xlabel="position", ylabel="downsampled errors", yscale="log")
            else:
                ax.set(xlabel="position", ylabel="error rate", yscale="log")

            axtwin = ax.twinx()
            axtwin.plot(
                df["position"],
                df["num consensus molecules"],
                label="consensus",
                color="black",
                alpha=0.25,
            )
            axtwin.set(ylabel="number of consensus molecules")
            for i, seq in seq_df.iterrows():
                axtwin.plot(
                    [seq["start"], seq["end"]],
                    [0, 0],
                    label="seq {}".format(i),
                    marker="|",
                )
        axs[-1].legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")

        fig.suptitle("{} {}".format(plot_title, base), size=16, y=0.52)
        fig.subplots_adjust(top=0.8)
        fig.tight_layout()
        if downsample:
            fig.savefig(
                "plots\\downsampled_errors\\{}_{}_downsampled.png".format(
                    plot_title, base
                )
            )
        else:
            fig.savefig(
                "plots\\error_rates\\{}_{}_error_rate.png".format(plot_title, base)
            )

    plt.close("all")
