import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from constants import row_age_map
from lookup import chromosome_tiles_map, exon_tiles_map, tile_df, tile_data_df, juicy_df
from constants import VARIANTS
from scipy.stats import betabinom, poisson
from scipy.optimize import curve_fit, fsolve, newton
from collections.abc import Iterable


def mean_var(tile_numbers, variant, trim_and_flip=True, age="all"):
    """Returns mean and var for a given tile and variant."""
    df = pd.DataFrame()
    if not isinstance(tile_numbers, Iterable):
        tile_numbers = [tile_numbers]
    for tile_number in tile_numbers:
        df = df.append(tile_data_df(tile_number, trim_and_flip=trim_and_flip))

    if True == df.empty:
        return 0, 0
    df = df.loc[df["variant"] == variant]
    if True == df.empty:
        return 0, 0
    if age != "all":
        df = df.loc[row_age_map(df["sample ID"]) == age]

    mean = np.mean(df["downsample"])
    variance = np.var(df["downsample"])
    return mean, variance


def position_mean_var(tile_number, position, trim_and_flip=True):
    """Returns mean and var for a given position."""
    df = tile_data_df(tile_number, trim_and_flip=trim_and_flip)

    df = df.loc[df["position"] == position]
    return np.mean(df["downsample"]), np.var(df["downsample"])


def plot_tile_mean_var(tile_number, save=False, age="all"):
    """Plots all the means and variances for all positions on a given tile."""
    for variant in VARIANTS:
        print(variant)
        fig, axs = plt.subplots(3, figsize=(8, 8))

        df = tile_data_df(tile_number)
        df = df.loc[df["variant"] == variant]
        if age != "all":
            df = df.loc[row_age_map(df["sample ID"]) == age]

        positions = pd.unique(df["position"])
        length = len(positions)
        means = np.zeros(length)
        variances = np.zeros(length)

        for i, position in enumerate(positions):
            means[i], variances[i] = position_mean_var(tile_number, position)

        non_zero_means = means[np.nonzero(means)]
        non_zero_variances = variances[np.nonzero(means)]
        non_zero_positions = positions[np.nonzero(means)]

        Ds = non_zero_variances / non_zero_means

        marker = "+"

        axs[0].plot(
            positions, means, label="means", marker=marker, linestyle="None",
        )

        axs[1].plot(
            positions, variances, label="variances", marker=marker, linestyle="None",
        )

        axs[2].plot(
            non_zero_positions,
            Ds,
            label="index of dispersion",
            marker=marker,
            linestyle="None",
        )

        axs[0].set(
            title="{} means".format(variant), xlabel="position", ylabel="mean",
        )
        axs[1].set(title="variances", xlabel="position", ylabel="variance")
        axs[2].set(
            title="Index of dispersion", xlabel="position", ylabel="D = Var/mean",
        )
        for ax in axs:
            ax.ticklabel_format(useOffset=False, style="plain")
            ax.set(yscale="log")

        fig.tight_layout()
        if save:
            file_name = "plots\\means_and_variances\\{}_mean_var_positions".format(
                variant
            )
            if show_strands:
                file_name += "_strands"
            if trim_and_flip:
                file_name += "_t&f"
            if group_chromosomes:
                file_name += "_group_chromosomes"
            fig.savefig(file_name + ".png")
            fig.savefig(file_name + ".svg", dpi=1200)
        else:
            plt.show()
        plt.close("all")


def plot_all_mean_var(
    show_strands=True,
    save=True,
    trim_and_flip=True,
    group_chromosomes=False,
    age="all",
    chromosome="all",
):
    """Plots the mean, variance and the index of dispersion for all tiles."""
    for variant in VARIANTS:
        print(variant)
        fig, axs = plt.subplots(3, figsize=(8, 8))

        for strand in ["+", "-"]:

            for j, chromosome in enumerate(pd.unique(tile_df["chromosome"])):
                print(j)
                chr_tile_df = tile_df.loc[
                    (tile_df["chromosome"] == chromosome)
                    & (tile_df["strand"] == strand)
                ]

                if group_chromosomes:
                    means, variances = mean_var(
                        chr_tile_df.index, variant, trim_and_flip, age
                    )
                    xs = j
                else:
                    means = np.zeros(len(chr_tile_df.index))
                    variances = np.zeros(len(chr_tile_df.index))
                    for i, index in enumerate(chr_tile_df.index):
                        print(i)
                        means[i], variances[i] = mean_var(
                            index, variant, trim_and_flip, age
                        )
                    xs = chr_tile_df.index  # janky trick to get x axis to work
                if show_strands:
                    marker = "${}$".format(strand)
                else:
                    marker = "${}$".format(chromosome)

                axs[0].plot(
                    xs, means, label="means", marker=marker, linestyle="None",
                )

                axs[1].plot(
                    xs, variances, label="variances", marker=marker, linestyle="None",
                )
                if not group_chromosomes:
                    xs = xs[np.nonzero(means)]
                    variances = variances[np.nonzero(means)]
                    means = means[np.nonzero(means)]

                axs[2].plot(
                    xs,
                    variances / means,
                    label="index of dispersion",
                    marker=marker,
                    linestyle="None",
                )

        axs[0].set(
            title="{} means".format(variant), xlabel="tile", ylabel="mean", yscale="log"
        )
        axs[1].set(title="variances", xlabel="tile", ylabel="variance", yscale="log")
        axs[2].set(
            title="Index of dispersion",
            xlabel="tile",
            ylabel="D = Var/mean",
            yscale="log",
        )
        if age != "all":
            axs[0].set(title="{} means, age {}".format(variant, age))
        fig.tight_layout()
        if save:
            if age != "all":
                file_name = "plots\\means_and_variances\\split_ages\\{}_mean_var".format(
                    variant
                )
            else:
                file_name = "plots\\means_and_variances\\{}_mean_var".format(variant)
            if show_strands:
                file_name += "_strands"
            if trim_and_flip:
                file_name += "_t&f"
            if group_chromosomes:
                file_name += "_group_chromosomes"
            if age != "all":
                file_name += "_age_{}".format(age)
            fig.savefig(file_name + ".png")
            fig.savefig(file_name + ".svg", dpi=1200)
        else:
            plt.show()
        plt.close("all")


def plot_exon_mean_var(exon_number, save=True, trim_and_flip=True):
    df = exon_tiles_map[exon_number]
    means = np.zeros(len(df.index))
    variances = np.zeros(len(df.index))

    for variant in VARIANTS:
        print(variant)
        fig, axs = plt.subplots(3, figsize=(8, 8))
        for i in df.index:
            print(i)
            means[i], variances[i] = mean_var(i, variant, trim_and_flip)
        axs[0].plot(df.index, means, label="means", marker="+", linestyle="None")
        axs[0].set(title="means", xlabel="tile", ylabel="mean", yscale="log")

        axs[1].plot(
            df.index, variances, label="variances", marker="+", linestyle="None"
        )
        axs[1].set(title="variances", xlabel="tile", ylabel="variance", yscale="log")

        axs[2].plot(
            df.index, variances / means, label="means", marker="+", linestyle="None"
        )
        axs[2].set(title="ratios", xlabel="tile", ylabel="ratios", yscale="log")
        fig.tight_layout()
        if save:
            fig.savefig("plots\\means_and_variances\\{}_mean_var.png".format(variant))
        else:
            plt.show()
        plt.close("all")


def plot_chromosome_variant_hist(
    chromosome, fit=None, save=True, strand=None, bins_to_fit=-1
):
    """
    Plots a histogram of the number of downsampled variants for a given tile

    fit = None/"Poisson"/"beta-binomial"

    strand = None/"+"/"-"
    """
    tile_numbers = chromosome_tiles_map[str(chromosome)].index
    df = pd.DataFrame()
    for tile_number in tile_numbers:
        if strand is None or tile_df.at[tile_number, "strand"] == strand:
            df = df.append(tile_data_df(tile_number))
            print(tile_number)

    for variant in VARIANTS:
        print(variant)
        fig, ax = plt.subplots()
        change_df = df.loc[variant == df["variant"]]
        mean = np.mean(change_df["downsample"])
        variance = np.var(change_df["downsample"])
        D = variance / mean  # index of dispersion
        print("mean: {}, variance: {}, variance/mean: {}".format(mean, variance, D))
        n = 6348  # 50th percentile of smaller data file
        N = len(change_df.index) * n  # To adjust normalisation of distributions

        maximum = np.amax(change_df["downsample"])
        print(maximum)

        bins = np.arange(-0.5, maximum + 1.5)
        xs = np.arange(maximum + 1)

        hs, hs_bin_edges = np.histogram(change_df["downsample"], bins)
        print(hs)

        ax.hist(
            change_df["downsample"], bins=bins, color="c", linestyle="-", edgecolor="k",
        )

        if fit == "Poisson":
            ys = poisson.pmf(xs, mean) * N
            ax.plot(xs, ys, color="k", marker="+")

        if fit == "beta-binomial":

            def b(a):
                return a * (n / mean - 1.0)

            def f(x, a):
                return N * betabinom.pmf(x, n, a, b(a))

            a, pcov = curve_fit(f, xs[:bins_to_fit], hs[:bins_to_fit])
            fit_mean = betabinom.mean(n, a, b(a))
            print("fit mean: {}".format(fit_mean))
            fit_var = betabinom.var(n, a, b(a))
            print("fit variance: {}".format(fit_var))

            ax.plot(xs, f(xs, a), color="k", marker="+")

        plot_title = variant
        if strand is not None:
            plot_title += " " + strand
        ax.set(
            title=plot_title,
            xlabel="number of variants",
            ylabel="frequency",
            yscale="log",
        )
        ax.text(0.8, 0.9, "D = {:.2f}".format(D), transform=ax.transAxes)
        if save:
            file_name = "plots\\variant_histograms\\variant_hist_chr{}_{}".format(
                chromosome, variant
            )
            if fit is not None:
                file_name += "_" + fit
            if strand is not None:
                file_name += "_" + strand
            fig.savefig(file_name + ".png")
            fig.savefig(file_name + ".svg", dpi=1200)
        else:
            plt.show()


def plot_juicy_hist(fit=None, bins_to_fit=-1):
    """Plots histograms of the juiciest positions and variants."""

    for i, juicy_row in juicy_df.iterrows():
        juicy_data_df = tile_data_df(juicy_row["tile"])
        juicy_data_df = juicy_data_df.loc[
            (juicy_data_df["position"] == juicy_row["position"])
            & (juicy_data_df["variant"] == juicy_row["variant"])
        ]

        df = juicy_data_df
        plot_title = "chr {}, position {}, {}".format(
            juicy_row["chromosome"], juicy_row["position"], juicy_row["variant"]
        )

        n = 6348
        N = len(df.index)  # To adjust normalisation of distributions
        print("N = {}".format(N))
        mean = np.mean(df["downsample"])
        print(mean)

        fig, ax = plt.subplots()

        maximum = np.amax(df["downsample"])
        bins = np.arange(-0.5, maximum + 1.5)
        hs, hs_bin_edges = np.histogram(df["downsample"], bins)

        ax.hist(
            df["downsample"],
            bins=bins[: bins_to_fit + 1],
            color="green",
            linestyle="-",
            edgecolor="k",
        )
        ax.hist(
            df["downsample"],
            bins=bins[bins_to_fit:],
            color="c",
            linestyle="-",
            edgecolor="k",
        )
        xs = np.arange(maximum + 1)

        if fit == "Poisson fix mean" or fit == "all":
            ys = poisson.pmf(xs, mean) * N
            ax.plot(
                xs, ys, marker="+", color="xkcd:piss yellow", label="Poisson fit mean"
            )

        if fit == "Poisson" or fit == "all":

            def f(x, mean):
                return poisson.pmf(x, mean) * N

            mean, pcov = curve_fit(f, xs[:bins_to_fit], hs[:bins_to_fit])
            ys = f(xs, mean)
            ax.plot(xs, ys, marker="+", color="xkcd:puke green", label="Poisson")

        if fit == "beta-binomial" or fit == "both beta-binomial" or fit == "all":

            def b(a):
                return a * (n / mean - 1.0)

            def f(x, a):
                return betabinom.pmf(x, n, a, b(a)) * N

            a, pcov = curve_fit(f, xs[:bins_to_fit], hs[:bins_to_fit])
            fit_mean = betabinom.mean(n, a, b(a))
            print("fit mean: {}".format(fit_mean))
            fit_var = betabinom.var(n, a, b(a))
            print("fit variance: {}".format(fit_var))

            ax.plot(xs, f(xs, a), marker="+", color="k", label="beta-binomial")

        if (
            fit == "beta-binomial fix mean"
            or fit == "both beta-binomial"
            or fit == "all"
        ):

            def f(x, a, b):
                return betabinom.pmf(x, n, a, b) * N

            (a, b), pcov = curve_fit(f, xs[:bins_to_fit], hs[:bins_to_fit])
            fit_mean = betabinom.mean(n, a, b)
            print("fit mean: {}".format(fit_mean))
            fit_var = betabinom.var(n, a, b)
            print("fit variance: {}".format(fit_var))

            ax.plot(
                xs, f(xs, a, b), marker="+", color="r", label="beta-binomial fixed mean"
            )

        ax.plot([xs[0], xs[-1]], [1.0, 1.0], label="1/N")

        ax.set(title=plot_title, yscale="log")
        ax.legend()

        def func(x):
            return betabinom.pmf(x, n, a, b) * N - 1

        for x in xs:
            if func(x) <= 0:
                x
                break

        print(xs)
        print(func(xs))
        print("first unlikely data: {}".format(x))
        plt.show()

