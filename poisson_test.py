import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lookup import gene_seqs_map, seq_df, seq_data_df
from constants import SUBS


def mean_var(sequence_number, sub, chromosome):
    """
    Returns  of mean and var for a given sequence and sub.
    """
    df = seq_data_df(sequence_number, trimmed_and_flipped=True)
    df = df.loc[df["sub"] == sub]

    mean = np.mean(df["downsample"])
    variance = np.var(df["downsample"])
    return mean, variance


def plot_all_mean_var():
    for sub in SUBS:
        print(sub)
        fig, ax = plt.subplots(3, figsize=(8, 8))

        for chromosome in pd.unique(seq_df["chromosome"]):
            chr_seq_df = seq_df.loc[seq_df["chromosome"] == chromosome]
            means = np.zeros(len(chr_seq_df.index))
            variances = np.zeros(len(chr_seq_df.index))
            for i, index in enumerate(chr_seq_df.index):
                print(i)
                means[i], variances[i] = mean_var(index, sub, chromosome)

            ax[0].plot(
                chr_seq_df.index,
                means,
                label="means",
                marker="${}$".format(chromosome),
                linestyle="None",
            )

            ax[1].plot(
                chr_seq_df.index,
                variances,
                label="variances",
                marker="${}$".format(chromosome),
                linestyle="None",
            )

            ax[2].plot(
                chr_seq_df.index,
                variances / means,
                label="means",
                marker="${}$".format(chromosome),
                linestyle="None",
            )
        ax[0].set(title="means", xlabel="sequence", ylabel="mean", yscale="log")
        ax[1].set(title="variances", xlabel="sequence", ylabel="variance", yscale="log")
        ax[2].set(title="ratios", xlabel="sequence", ylabel="ratios", yscale="log")
        fig.tight_layout()
        fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(sub))
        # plt.show()
        plt.close("all")


def plot_gene_mean_var(gene_number):
    df = gene_seqs_map[gene_number]
    means = np.zeros(len(df.index))
    variances = np.zeros(len(df.index))

    for sub in SUBS:
        print(sub)
        fig, ax = plt.subplots(3, figsize=(8, 8))
        for i in df.index:
            print(i)
            means[i], variances[i] = mean_var(i, sub)
        ax[0].plot(df.index, means, label="means", marker="+", linestyle="None")
        ax[0].set(title="means", xlabel="sequence", ylabel="mean", yscale="log")

        ax[1].plot(df.index, variances, label="variances", marker="+", linestyle="None")
        ax[1].set(title="variances", xlabel="sequence", ylabel="variance", yscale="log")

        ax[2].plot(
            df.index, variances / means, label="means", marker="+", linestyle="None"
        )
        ax[2].set(title="ratios", xlabel="sequence", ylabel="ratios", yscale="log")
        fig.tight_layout()
        # fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(sub))
        plt.show()
        plt.close("all")