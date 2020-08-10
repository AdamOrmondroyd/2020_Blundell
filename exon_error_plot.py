"""
# 2020 Blundell lab internship

Generates plots of error rates for given exons
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from constants import BASES, base_variants_map, variant_color_map
from lookup import exon_df, exon_tiles_map, tile_data_df


def exon_error_plot(exon_number, downsample=True, trim_and_flip=True, save=True):
    """
    Saves plots of errors for a given exon, separating + and - strand data.

    downsample = True will plot the downsampled number of errors.
    downsample = False will plot the error rate for the full number of errors.
    """

    print(exon_number)
    exon = exon_df.loc[exon_number, :]
    tile_df = exon_tiles_map[exon_number]
    df = pd.DataFrame()
    plot_title = "Gene {} ".format(exon_number)

    if trim_and_flip:
        plot_title += "flipped "
        for i in tile_df.index:
            df = pd.concat([df, tile_data_df(i, group_by="position")]).drop_duplicates(
                keep="first"
            )
    else:

        for i in tile_df.index:
            df = pd.concat(
                [df, tile_data_df(i, group_by="position", trim_and_flip=False)]
            ).drop_duplicates(keep="first")

    # Make up plot title
    for strand in tile_df["strand"]:
        plot_title += strand

    for base in BASES:
        fig, axs = plt.subplots(2, 2, figsize=(15, 8))

        axs = axs.flatten()
        for j, ax in enumerate(axs):
            for i, variant in enumerate(base_variants_map[base]):
                if (3 == j) or (i == j):
                    color = variant_color_map[variant]
                    if 3 == j:
                        ax.set(title=base)
                    else:
                        ax.set(title=variant)
                else:
                    color = "lightgrey"

                variant_df = df.loc[df["variant"] == variant]
                if downsample:
                    ax.plot(
                        variant_df["position"],
                        variant_df["downsample"],
                        label=variant + "+",
                        linestyle="None",
                        marker="^",
                        color=color,
                        alpha=0.5,
                    )
                else:
                    ax.plot(
                        variant_df["position"],
                        variant_df["num variants"]
                        / variant_df["num consensus molecules"],
                        label=variant + "+",
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
            for i, tile in tile_df.iterrows():
                axtwin.plot(
                    [tile["start"], tile["end"]],
                    [0, 0],
                    label="tile {}".format(i),
                    marker="|",
                )
        axs[-1].legend(bbox_to_anchor=(1.05, 1.0), loc="upper left", frameon=False)

        fig.suptitle("{} {}".format(plot_title, base), size=16, y=0.52)
        fig.subplots_adjust(top=0.8)
        fig.tight_layout()
        if save:
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
    if ~save:
        plt.show()

    plt.close("all")
