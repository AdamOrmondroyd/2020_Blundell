"""
Generates plots of downsampled errors for given genes.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from constants import BASES, CHANGES, base_subs_map, sub_color_map
from lookup import gene_df, gene_seqs_map, seq_data_df

for gene_number in np.arange(0, 10):
    # for gene_number in [6]:
    print(gene_number)
    gene = gene_df.loc[gene_number, :]

    plot_title = "Gene {} ".format(gene_number)

    seq_df = gene_seqs_map[gene_number]

    consensus_df = pd.DataFrame()
    for i in seq_df.index:
        consensus_df = pd.concat(
            [consensus_df, seq_data_df(i, group_by="position")]
        ).drop_duplicates(keep="first")

    pos_seq_df = seq_df.loc[seq_df["strand"] == "+"]
    neg_seq_df = seq_df.loc[seq_df["strand"] == "-"]

    # Bring together the data for each sequence, separated by pos and neg strand, grouped by position
    pos = False
    if len(pos_seq_df.index):
        pos = True
        pos_df = pd.DataFrame()
        for i in pos_seq_df.index:
            pos_df = pd.concat(
                [pos_df, seq_data_df(i, group_by="position")]
            ).drop_duplicates(keep=False)

    negative = False
    if len(neg_seq_df.index):
        neg = True
        neg_df = pd.DataFrame()
        for i in neg_seq_df.index:
            neg_df = pd.concat(
                [neg_df, seq_data_df(i, group_by="position")]
            ).drop_duplicates(keep=False)

    # Drop rows that appear in the other strand
    if pos and neg:
        pos_cond = ~pos_df["position"].isin(neg_df["position"])
        neg_cond = ~neg_df["position"].isin(pos_df["position"])
        pos_df = pos_df.loc[pos_cond, :]
        neg_df = neg_df.loc[neg_cond, :]

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
            ax.ticklabel_format(useOffset=False, style="plain")
            ax.set(xlabel="position", ylabel="downsampled errors", yscale="log")

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
        fig.savefig(
            "plots\\downsampled_errors\\{}_{}_downsampled.png".format(plot_title, base)
        )

    plt.close("all")
