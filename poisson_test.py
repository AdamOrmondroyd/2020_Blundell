import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lookup import exon_seqs_map, seq_df, seq_data_df
from constants import SUBS


def mean_var(seq_number, sub, chromosome):
    """
    Returns  of mean and var for a given seq and sub.
    """
    df = seq_data_df(seq_number, trim_and_flip=False)
    df = df.loc[df["sub"] == sub]

    mean = np.mean(df["downsample"])
    variance = np.var(df["downsample"])
    return mean, variance


def plot_all_mean_var(save=True):
    for sub in SUBS:
        print(sub)
        fig, ax = plt.subplots(3, figsize=(8, 8))

        for chromosome in pd.unique(seq_df["chromosome"]):
            chr_seq_df = seq_df.loc[
                (seq_df["chromosome"] == chromosome) & (seq_df["strand"] == "+")
            ]
            means = np.zeros(len(chr_seq_df.index))
            variances = np.zeros(len(chr_seq_df.index))
            for i, index in enumerate(chr_seq_df.index):
                print(i)
                means[i], variances[i] = mean_var(index, sub, chromosome)

            ax[0].plot(
                chr_seq_df.index,
                means,
                label="means",
                marker="${}$".format("+"),
                linestyle="None",
            )

            ax[1].plot(
                chr_seq_df.index,
                variances,
                label="variances",
                marker="${}$".format("+"),
                linestyle="None",
            )

            ax[2].plot(
                chr_seq_df.index,
                variances / means,
                label="means",
                marker="${}$".format("+"),
                linestyle="None",
            )
        for chromosome in pd.unique(seq_df["chromosome"]):
            chr_seq_df = seq_df.loc[
                (seq_df["chromosome"] == chromosome) & (seq_df["strand"] == "-")
            ]
            means = np.zeros(len(chr_seq_df.index))
            variances = np.zeros(len(chr_seq_df.index))
            for i, index in enumerate(chr_seq_df.index):
                print(i)
                means[i], variances[i] = mean_var(index, sub, chromosome)

            ax[0].plot(
                chr_seq_df.index,
                means,
                label="means",
                marker="${}$".format("-"),
                linestyle="None",
            )

            ax[1].plot(
                chr_seq_df.index,
                variances,
                label="variances",
                marker="${}$".format("-"),
                linestyle="None",
            )

            ax[2].plot(
                chr_seq_df.index,
                variances / means,
                label="means",
                marker="${}$".format("-"),
                linestyle="None",
            )
        ax[0].set(
            title="{} means".format(sub), xlabel="seq", ylabel="mean", yscale="log"
        )
        ax[1].set(title="variances", xlabel="seq", ylabel="variance", yscale="log")
        ax[2].set(title="ratios", xlabel="seq", ylabel="ratios", yscale="log")
        fig.tight_layout()
        if save:
            fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(sub))
        else:
            plt.show()
        plt.close("all")


def plot_exon_mean_var(exon_number, save=True):
    df = exon_seqs_map[exon_number]
    means = np.zeros(len(df.index))
    variances = np.zeros(len(df.index))

    for sub in SUBS:
        print(sub)
        fig, ax = plt.subplots(3, figsize=(8, 8))
        for i in df.index:
            print(i)
            means[i], variances[i] = mean_var(i, sub)
        ax[0].plot(df.index, means, label="means", marker="+", linestyle="None")
        ax[0].set(title="means", xlabel="seq", ylabel="mean", yscale="log")

        ax[1].plot(df.index, variances, label="variances", marker="+", linestyle="None")
        ax[1].set(title="variances", xlabel="seq", ylabel="variance", yscale="log")

        ax[2].plot(
            df.index, variances / means, label="means", marker="+", linestyle="None"
        )
        ax[2].set(title="ratios", xlabel="seq", ylabel="ratios", yscale="log")
        fig.tight_layout()
        if save:
            fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(sub))
        else:
            plt.show()
        plt.close("all")


def plot_seq_sub_hist(seq_number):
    """
    Plots a histogram of the number of downsampled substitutions for a given seq
    """
    plot_title = "histogram"

    df = seq_data_df(seq_number)

    for sub in SUBS:
        print(sub)
        fig, ax = plt.subplots()

        change_df = df.loc[sub == df["sub"]]

        bins = np.arange(-0.5, np.max(change_df["num subs"]) + 0.5)
        print(bins)
        ax.hist(change_df["num subs"], bins=bins)
        ax.set(title="{}_{}".format(plot_title, sub), yscale="log")
        plt.show()
