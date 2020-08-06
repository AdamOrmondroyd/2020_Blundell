import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lookup import gene_exons_map, exon_df, exon_data_df
from constants import SUBS


def mean_var(exon_number, sub, chromosome):
    """
    Returns  of mean and var for a given exon and sub.
    """
    df = exon_data_df(exon_number, trim_and_flip=True)
    df = df.loc[df["sub"] == sub]

    mean = np.mean(df["downsample"])
    variance = np.var(df["downsample"])
    return mean, variance


def plot_all_mean_var():
    for sub in SUBS:
        print(sub)
        fig, ax = plt.subplots(3, figsize=(8, 8))

        for chromosome in pd.unique(exon_df["chromosome"]):
            chr_exon_df = exon_df.loc[
                (exon_df["chromosome"] == chromosome) & (exon_df["strand"] == "+")
            ]
            means = np.zeros(len(chr_exon_df.index))
            variances = np.zeros(len(chr_exon_df.index))
            for i, index in enumerate(chr_exon_df.index):
                print(i)
                means[i], variances[i] = mean_var(index, sub, chromosome)

            ax[0].plot(
                chr_exon_df.index,
                means,
                label="means",
                marker="${}$".format("+"),
                linestyle="None",
            )

            ax[1].plot(
                chr_exon_df.index,
                variances,
                label="variances",
                marker="${}$".format("+"),
                linestyle="None",
            )

            ax[2].plot(
                chr_exon_df.index,
                variances / means,
                label="means",
                marker="${}$".format("+"),
                linestyle="None",
            )
        for chromosome in pd.unique(exon_df["chromosome"]):
            chr_exon_df = exon_df.loc[
                (exon_df["chromosome"] == chromosome) & (exon_df["strand"] == "-")
            ]
            means = np.zeros(len(chr_exon_df.index))
            variances = np.zeros(len(chr_exon_df.index))
            for i, index in enumerate(chr_exon_df.index):
                print(i)
                means[i], variances[i] = mean_var(index, sub, chromosome)

            ax[0].plot(
                chr_exon_df.index,
                means,
                label="means",
                marker="${}$".format("-"),
                linestyle="None",
            )

            ax[1].plot(
                chr_exon_df.index,
                variances,
                label="variances",
                marker="${}$".format("-"),
                linestyle="None",
            )

            ax[2].plot(
                chr_exon_df.index,
                variances / means,
                label="means",
                marker="${}$".format("-"),
                linestyle="None",
            )
        ax[0].set(title="means", xlabel="exon", ylabel="mean", yscale="log")
        ax[1].set(title="variances", xlabel="exon", ylabel="variance", yscale="log")
        ax[2].set(title="ratios", xlabel="exon", ylabel="ratios", yscale="log")
        fig.tight_layout()
        # fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(sub))
        plt.show()
        plt.close("all")


def plot_gene_mean_var(gene_number):
    df = gene_exons_map[gene_number]
    means = np.zeros(len(df.index))
    variances = np.zeros(len(df.index))

    for sub in SUBS:
        print(sub)
        fig, ax = plt.subplots(3, figsize=(8, 8))
        for i in df.index:
            print(i)
            means[i], variances[i] = mean_var(i, sub)
        ax[0].plot(df.index, means, label="means", marker="+", linestyle="None")
        ax[0].set(title="means", xlabel="exon", ylabel="mean", yscale="log")

        ax[1].plot(df.index, variances, label="variances", marker="+", linestyle="None")
        ax[1].set(title="variances", xlabel="exon", ylabel="variance", yscale="log")

        ax[2].plot(
            df.index, variances / means, label="means", marker="+", linestyle="None"
        )
        ax[2].set(title="ratios", xlabel="exon", ylabel="ratios", yscale="log")
        fig.tight_layout()
        # fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(sub))
        plt.show()
        plt.close("all")
