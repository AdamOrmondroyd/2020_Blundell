import numpy as np
import pandas as pd
from copy import deepcopy
import matplotlib.pyplot as plt
from constants import BASES, CHANGES, base_subs_map, sub_color_map
from lookup import gene_df, gene_seqs_map, seq_data_df

for gene_number in [1]:
    print(gene_number)
    gene = gene_df.loc[gene_number, :]

    plot_title = "Gene {} ".format(gene_number)

    seq_df = gene_seqs_map[gene_number]

    df = pd.DataFrame()

    # Bring together the data for each sequence, grouped by position and dropping overlaps between sequences
    for i in seq_df.index:
        df = pd.concat([df, seq_data_df(i, group_by="position")]).drop_duplicates(
            keep=False
        )

    # Make up plot title
    for strand in seq_df["strand"]:
        plot_title += strand

    for base in BASES:
        fig, axs = plt.subplots(2, 2, figsize=(15, 8))

        axs = axs.flatten()
        for j, ax in enumerate(axs):
            for i, sub in enumerate(base_subs_map[base]):
                if 4 == i or i == j:
                    color = sub_color_map[sub]
                    if i == j:
                        ax.set(title=sub)
                else:
                    color = "lightgrey"

                df_sub = df.loc[df["sub"] == sub]

                ax.plot(
                    df_sub["position"],
                    df_sub["num subs"] / df_sub["num consensus molecules"],
                    label="sub",
                    linestyle="None",
                    marker="o",
                    color=color,
                    alpha=0.5,
                )
            ax.ticklabel_format(useOffset=False, style="plain")
            ax.set(xlabel="position", ylabel="error rate", yscale="log")

            axtwin = ax.twinx()
            axtwin.plot(
                df["position"],
                df["num consensus molecules"],
                label="consensus",
                linestyle="None",
                marker="+",
                markersize="1",
                color="black",
                alpha=0.25,
            )
            axtwin.set(ylabel="number of consensus molecules")
        axs[-1].legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")

        fig.suptitle("{} {}".format(plot_title, base), size=16, y=0.52)
        fig.subplots_adjust(top=0.8)
        fig.tight_layout()

    plt.show()
    plt.close()
