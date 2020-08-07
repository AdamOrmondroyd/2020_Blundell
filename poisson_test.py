import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lookup import exon_tiles_map, tile_df, tile_data_df
from constants import SUBS


def mean_var(tile_number, variant, chromosome):
    """
    Returns  of mean and var for a given tile and variant.
    """
    df = tile_data_df(tile_number, trim_and_flip=False)
    df = df.loc[df["variant"] == variant]

    mean = np.mean(df["downsample"])
    variance = np.var(df["downsample"])
    return mean, variance


def plot_all_mean_var(save=True):
    for variant in SUBS:
        print(variant)
        fig, ax = plt.subplots(3, figsize=(8, 8))

        for chromosome in pd.unique(tile_df["chromosome"]):
            chr_tile_df = tile_df.loc[
                (tile_df["chromosome"] == chromosome) & (tile_df["strand"] == "+")
            ]
            means = np.zeros(len(chr_tile_df.index))
            variances = np.zeros(len(chr_tile_df.index))
            for i, index in enumerate(chr_tile_df.index):
                print(i)
                means[i], variances[i] = mean_var(index, variant, chromosome)

            ax[0].plot(
                chr_tile_df.index,
                means,
                label="means",
                marker="${}$".format("+"),
                linestyle="None",
            )

            ax[1].plot(
                chr_tile_df.index,
                variances,
                label="variances",
                marker="${}$".format("+"),
                linestyle="None",
            )

            ax[2].plot(
                chr_tile_df.index,
                variances / means,
                label="means",
                marker="${}$".format("+"),
                linestyle="None",
            )
        for chromosome in pd.unique(tile_df["chromosome"]):
            chr_tile_df = tile_df.loc[
                (tile_df["chromosome"] == chromosome) & (tile_df["strand"] == "-")
            ]
            means = np.zeros(len(chr_tile_df.index))
            variances = np.zeros(len(chr_tile_df.index))
            for i, index in enumerate(chr_tile_df.index):
                print(i)
                means[i], variances[i] = mean_var(index, variant, chromosome)

            ax[0].plot(
                chr_tile_df.index,
                means,
                label="means",
                marker="${}$".format("-"),
                linestyle="None",
            )

            ax[1].plot(
                chr_tile_df.index,
                variances,
                label="variances",
                marker="${}$".format("-"),
                linestyle="None",
            )

            ax[2].plot(
                chr_tile_df.index,
                variances / means,
                label="means",
                marker="${}$".format("-"),
                linestyle="None",
            )
        ax[0].set(
            title="{} means".format(variant), xlabel="tile", ylabel="mean", yscale="log"
        )
        ax[1].set(title="variances", xlabel="tile", ylabel="variance", yscale="log")
        ax[2].set(title="ratios", xlabel="tile", ylabel="ratios", yscale="log")
        fig.tight_layout()
        if save:
            fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(variant))
        else:
            plt.show()
        plt.close("all")


def plot_exon_mean_var(exon_number, save=True):
    df = exon_tiles_map[exon_number]
    means = np.zeros(len(df.index))
    variances = np.zeros(len(df.index))

    for variant in SUBS:
        print(variant)
        fig, ax = plt.subplots(3, figsize=(8, 8))
        for i in df.index:
            print(i)
            means[i], variances[i] = mean_var(i, variant)
        ax[0].plot(df.index, means, label="means", marker="+", linestyle="None")
        ax[0].set(title="means", xlabel="tile", ylabel="mean", yscale="log")

        ax[1].plot(df.index, variances, label="variances", marker="+", linestyle="None")
        ax[1].set(title="variances", xlabel="tile", ylabel="variance", yscale="log")

        ax[2].plot(
            df.index, variances / means, label="means", marker="+", linestyle="None"
        )
        ax[2].set(title="ratios", xlabel="tile", ylabel="ratios", yscale="log")
        fig.tight_layout()
        if save:
            fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(variant))
        else:
            plt.show()
        plt.close("all")


def plot_tile_variant_hist(tile_number):
    """
    Plots a histogram of the number of downsampled variantstitutions for a given tile
    """
    plot_title = "histogram"

    df = tile_data_df(tile_number)

    for variant in SUBS:
        print(variant)
        fig, ax = plt.subplots()

        change_df = df.loc[variant == df["variant"]]

        bins = np.arange(-0.5, np.max(change_df["num variants"]) + 0.5)
        print(bins)
        ax.hist(change_df["num variants"], bins=bins)
        ax.set(title="{}_{}".format(plot_title, variant), yscale="log")
        plt.show()
