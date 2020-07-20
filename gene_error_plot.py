"""
Generates plots of error rates for given genes.
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from constants import BASES, CHANGES, base_subs_map, sub_color_map
from lookup import gene_df, gene_seqs_map, seq_data_df


def gene_error_plot(gene_number, downsample=True):
    """
    Saves plots of errors for a given gene, separating + and - strand data.

    downsample = True will plot the downsampled number of errors.
    downsample = False will plot the error rate for the full number of errors.
    """

    print(gene_number)
    gene = gene_df.loc[gene_number, :]

    plot_title = "Gene {} ".format(gene_number)

    seq_df = gene_seqs_map[gene_number]

    consensus_df = pd.DataFrame()
    for i in seq_df.index:
        consensus_df = pd.concat(
            [consensus_df, seq_data_df(i, group_by="position")]
        ).drop_duplicates(keep="first")

    # Bring in data from pos and neg strands of the gene
    pos = False
    if os.path.isfile("data_files\\genes\\gene_{}_pos.csv".format(gene_number)):
        pos = True
        pos_df = pd.read_csv("data_files\\genes\\gene_{}_pos.csv".format(gene_number))
    neg = False
    if os.path.isfile("data_files\\genes\\gene_{}_neg.csv".format(gene_number)):
        neg = True
        neg_df = pd.read_csv("data_files\\genes\\gene_{}_neg.csv".format(gene_number))

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

                if downsample:
                    if pos:
                        pos_sub_df = pos_df.loc[pos_df["sub"] == sub]
                        ax.plot(
                            pos_sub_df["position"],
                            pos_sub_df["downsample"],
                            label=sub + "+",
                            linestyle="None",
                            marker="^",
                            color=color,
                            alpha=0.5,
                        )
                    if neg:
                        neg_sub_df = neg_df.loc[neg_df["sub"] == sub]
                        ax.plot(
                            neg_sub_df["position"],
                            neg_sub_df["downsample"],
                            label=sub + "-",
                            linestyle="None",
                            marker="v",
                            color=color,
                            alpha=0.5,
                        )
                else:
                    if pos:
                        pos_sub_df = pos_df.loc[pos_df["sub"] == sub]
                        ax.plot(
                            pos_sub_df["position"],
                            pos_sub_df["num subs"]
                            / pos_sub_df["num consensus molecules"],
                            label=sub + "+",
                            linestyle="None",
                            marker="^",
                            color=color,
                            alpha=0.5,
                        )
                    if neg:
                        neg_sub_df = neg_df.loc[neg_df["sub"] == sub]
                        ax.plot(
                            neg_sub_df["position"],
                            neg_sub_df["num subs"]
                            / neg_sub_df["num consensus molecules"],
                            label=sub + "-",
                            linestyle="None",
                            marker="v",
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
                consensus_df["position"],
                consensus_df["num consensus molecules"],
                label="consensus",
                color="black",
                alpha=0.25,
            )
            axtwin.set(ylabel="number of consensus molecules")
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
