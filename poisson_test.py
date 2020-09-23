import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from constants import row_age_map, variant_contexts_map
from lookup import (
    chromosomes,
    chromosome_tiles_map,
    exon_tiles_map,
    tile_df,
    tile_data_df,
    juicy_df,
)
from constants import VARIANTS
from scipy.stats import betabinom, poisson
from scipy.optimize import curve_fit, fsolve, newton
from collections.abc import Iterable


def mean_var(tile_numbers, variant, trim_and_flip=True, age="all", context=None):
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
    if context is not None:
        df = df.loc[df["context"] == context]

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
    var_against_mean=False,
    show_exon_numbers=False,
    yscale="log",
    context=False,
):
    """Plots the mean, variance and the index of dispersion for all tiles."""
    for variant in VARIANTS:
        print(variant)
        if context:
            contexts_to_iterate = variant_contexts_map[variant]
        else:
            contexts_to_iterate = ["no context"]

        for context_string in contexts_to_iterate:
            print(context_string)

            if var_against_mean:
                fig, ax = plt.subplots()
            else:
                fig, axs = plt.subplots(3, figsize=(8, 8))

            if chromosome == "all":
                chrs_to_enumerate = chromosomes
            else:
                chrs_to_enumerate = chromosome

            for strand, color in zip(["+", "-"], ["blue", "orange"]):
                for j, chromosome in enumerate(chrs_to_enumerate):
                    # print(j)
                    chr_tile_df = tile_df.loc[
                        (tile_df["chromosome"] == chromosome)
                        & (tile_df["strand"] == strand)
                    ]
                    if show_exon_numbers:
                        marker_style = dict(
                            markersize=10, mfc="C0", mec=None, linewidth=0.5,
                        )
                        # marker_style.update(mec="None", markersize=100)
                        for exon_number in pd.unique(chr_tile_df["exon"]):
                            chr_exon_tile_df = chr_tile_df.loc[
                                chr_tile_df["exon"] == exon_number
                            ]

                            means = np.zeros(len(chr_exon_tile_df.index))
                            variances = np.zeros(len(chr_exon_tile_df.index))
                            for i, index in enumerate(chr_exon_tile_df.index):
                                # print(i)
                                means[i], variances[i] = mean_var(
                                    index, variant, trim_and_flip, age, context_string
                                )
                            xs = chr_exon_tile_df.index
                            marker = "${}$".format(exon_number)

                            if var_against_mean:
                                ax.plot(
                                    means,
                                    variances,
                                    label="variances against means",
                                    marker=marker,
                                    linestyle="None",
                                    color=color,
                                    **marker_style,
                                )
                            else:
                                axs[0].plot(
                                    xs,
                                    means,
                                    label="means",
                                    marker=marker,
                                    linestyle="None",
                                    color=color,
                                    **marker_style,
                                )

                                axs[1].plot(
                                    xs,
                                    variances,
                                    label="variances",
                                    marker=marker,
                                    linestyle="None",
                                    color=color,
                                    **marker_style,
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
                                    color=color,
                                    **marker_style,
                                )
                    else:
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

                        if var_against_mean:
                            ax.plot(
                                means,
                                variances,
                                label="variances against means",
                                marker=marker,
                                linestyle="None",
                            )
                        else:
                            axs[0].plot(
                                xs,
                                means,
                                label="means",
                                marker=marker,
                                linestyle="None",
                            )

                            axs[1].plot(
                                xs,
                                variances,
                                label="variances",
                                marker=marker,
                                linestyle="None",
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
                                color=color,
                            )
            if var_against_mean:
                axs0_title = "{} var against mean".format(variant)
            else:
                axs0_title = "{} means".format(variant)
            if context:
                axs0_title += " {}".format(context_string)
            if show_strands:
                axs0_title += " strands"
            if trim_and_flip:
                axs0_title += " t&f"
            if group_chromosomes:
                axs0_title += " group chromosomes"
            if age != "all":
                axs0_title += " age {}".format(age)
            if chromosome != "all":
                axs0_title += " {}".format(chromosome)
            if var_against_mean:
                axs0_title += " var against mean"
                ax.set(
                    title=axs0_title,
                    xlabel="mean",
                    ylabel="variance",
                    xscale="log",
                    yscale="log",
                )
            else:
                axs[0].set(
                    title=axs0_title, xlabel="tile", ylabel="mean", yscale=yscale
                )
                axs[1].set(
                    title="variances", xlabel="tile", ylabel="variance", yscale=yscale
                )
                axs[2].set(
                    title="Index of dispersion",
                    xlabel="tile",
                    ylabel="D = Var/mean",
                    yscale=yscale,
                )

            fig.tight_layout()
            if save:
                location = "plots\\means_and_variances\\"
                if age != "all":
                    location += "split_ages\\"
                if var_against_mean:
                    location += "var_against_mean\\"
                if context:
                    location += "context\\"

                file_name = "{}_mean_var".format(variant)
                if context:
                    file_name += "_{}".format(context_string)
                if show_strands:
                    file_name += "_strands"
                if trim_and_flip:
                    file_name += "_t&f"
                if group_chromosomes:
                    file_name += "_group_chromosomes"
                if age != "all":
                    file_name += "_age_{}".format(age)
                if chromosome != "all":
                    file_name += "_{}".format(chromosome)
                fig.savefig(location + file_name + ".png", dpi=600)
                fig.savefig(location + file_name + ".svg", dpi=1200)
            else:
                plt.show()
            plt.close("all")


def plot_all_mean_var_tile(
    save=True,
    trim_and_flip=True,
    age="all",
    chromosome="all",
    xscale="log",
    yscale="log",
    context=False,
):
    """Plots mean against variance for all tiles"""
    if context:
        contexts_to_iterate = CONTEXTS
    else:
        contexts_to_iterate = ["no context"]

    for context_string in contexts_to_iterate:
        for variant in VARIANTS:
            print(variant)

            fig, ax = plt.subplots()

            if chromosome == "all":
                chrs_to_use = chromosomes
            else:
                chrs_to_use = [chromosome]

            for strand, color in zip(["+", "-"], ["blue", "orange"]):
                strand_tile_df = tile_df.loc[
                    (tile_df["strand"] == strand)
                    & (tile_df["chromosome"].isin(chrs_to_use))
                ]

                means = []
                variances = []
                for tile_number in strand_tile_df.index:
                    data_df = tile_data_df(tile_number)
                    data_df = data_df.loc[data_df["variant"] == variant]
                    positions = data_df["position"].unique()
                    tile_means = np.zeros(len(positions))
                    tile_vars = np.zeros(len(positions))
                    for i, position in enumerate(positions):
                        df = data_df.loc[data_df["position"] == position]
                        tile_means[i] = np.mean(df["downsample"])
                        tile_vars[i] = np.var(df["downsample"])
                    means.append(tile_means)
                    variances.append(tile_vars)
                means = np.hstack(means)
                variances = np.hstack(variances)

                ax.plot(
                    means,
                    variances,
                    color=color,
                    marker="${}$".format(strand),
                    linestyle="None",
                )
                ax.set(xlabel="mean", ylabel="variance", xscale="log", yscale="log")
            plt.show()


def plot_len_tiles(show_strands=False, chromosome="all"):
    """Plots length of the tiles."""
    fig, ax = plt.subplots()

    plot_title = "tile lengths"
    if chromosome == "all":
        chrs_to_enumerate = chromosomes
    else:
        chrs_to_enumerate = chromosome
        plot_title += " {}".format(chromosome)

    for strand in ["+", "-"]:
        for chromosome in chrs_to_enumerate:
            chr_tile_df = tile_df.loc[
                (tile_df["chromosome"] == chromosome) & (tile_df["strand"] == strand)
            ]
            xs = chr_tile_df.index
            lengths = chr_tile_df["end"] - chr_tile_df["start"]

            if show_strands:
                marker = "${}$".format(strand)
            else:
                marker = "${}$".format(chromosome)

            ax.plot(xs, lengths, label="lengths", marker=marker, linestyle="None")
    ax.set(title=plot_title, xlabel="tile number", ylabel="length of tile")
    plt.show()


def plot_num_samples(show_strands=False, chromosome="all"):
    """Plots the number of people that survived the downsampling in each tile."""
    for variant in VARIANTS:
        print(variant)
        fig, ax = plt.subplots()

        plot_title = "num samples"
        if chromosome == "all":
            chrs_to_enumerate = chromosomes
        else:
            chrs_to_enumerate = chromosome
            plot_title += " {}".format(chromosome)

        for strand in ["+", "-"]:
            for chromosome in chrs_to_enumerate:
                chr_tile_df = tile_df.loc[
                    (tile_df["chromosome"] == chromosome)
                    & (tile_df["strand"] == strand)
                ]
                xs = chr_tile_df.index
                nums = np.zeros(len(xs))
                for x, i in zip(range(len(xs)), chr_tile_df.index):
                    print(x, i)
                    data_df = tile_data_df(i)
                    data_df = data_df.loc[data_df["variant"] == variant]
                    nums[x] = len(data_df.index)

                if show_strands:
                    marker = "${}$".format(strand)
                else:
                    marker = "${}$".format(chromosome)

                ax.plot(
                    xs, nums, label="number of samples", marker=marker, linestyle="None"
                )
        ax.set(
            title="{} ".format(variant) + plot_title,
            xlabel="tile number",
            ylabel="number of samples",
        )
        plt.show()


def plot_mean_var_against_num_people(show_strands=True, chromosome="all", save=True):
    """Plots the mean and variance against the number of people that survived the downsampling in each tile for that tile."""
    for variant in VARIANTS:
        print(variant)
        fig, axs = plt.subplots(1, 2, figsize=(7, 4))

        plot_title = "num samples"
        if chromosome == "all":
            chrs_to_enumerate = chromosomes
        else:
            chrs_to_enumerate = chromosome
            plot_title += " {}".format(chromosome)

        for strand in ["+", "-"]:
            for chromosome in chrs_to_enumerate:
                chr_tile_df = tile_df.loc[
                    (tile_df["chromosome"] == chromosome)
                    & (tile_df["strand"] == strand)
                ]
                xs = chr_tile_df.index
                nums = np.zeros(len(xs))
                for x, i in zip(range(len(xs)), chr_tile_df.index):
                    data_df = tile_data_df(i)
                    data_df = data_df.loc[data_df["variant"] == variant]
                    nums[x] = len(data_df.index)

                means = np.zeros(len(chr_tile_df.index))
                variances = np.zeros(len(chr_tile_df.index))
                for i, index in enumerate(chr_tile_df.index):
                    print(i)
                    means[i], variances[i] = mean_var(index, variant)

                if show_strands:
                    marker = "${}$".format(strand)
                else:
                    marker = "${}$".format(chromosome)

                axs[0].plot(
                    nums, means, label="means", marker=marker, linestyle="None",
                )
                axs[1].plot(
                    nums, variances, label="variances", marker=marker, linestyle="None",
                )
        axs[0].set(
            title="{} ".format(variant) + plot_title,
            xlabel="number of samples",
            ylabel="mean",
            yscale="log",
        )
        axs[1].set(
            title="{} ".format(variant) + plot_title,
            xlabel="number of samples",
            ylabel="variance",
            yscale="log",
        )
        if save:
            file_name = "plots\\means_and_variances\\num_people\\{}".format(variant)
            if chromosome != "all":
                file_name += "_num_people_{}".format(chromosome)
            else:
                file_name += "_num_people"
            fig.savefig(file_name + ".png", dpi=600)
            fig.savefig(file_name + ".svg", dpi=1200)

        else:
            plt.show()


def plot_mean_var_against_consensus(show_strands=True, chromosome="all", save=True):
    """Plots the mean and variance against the mean total number of consensus molecules per position for that tile."""
    for variant in VARIANTS:
        print(variant)
        fig, axs = plt.subplots(1, 2, figsize=(7, 4))

        plot_title = "num samples"
        if chromosome == "all":
            chrs_to_enumerate = chromosomes
        else:
            chrs_to_enumerate = chromosome
            plot_title += " {}".format(chromosome)

        for strand in ["+", "-"]:
            for chromosome in chrs_to_enumerate:
                chr_tile_df = tile_df.loc[
                    (tile_df["chromosome"] == chromosome)
                    & (tile_df["strand"] == strand)
                ]
                xs = chr_tile_df.index
                cons = np.zeros(len(xs))
                for x, i in zip(range(len(xs)), chr_tile_df.index):
                    data_df = tile_data_df(i, group_by="position")
                    if not data_df.empty:
                        data_df = data_df.loc[data_df["variant"] == variant]
                        cons[x] = np.mean(data_df["num consensus molecules"])
                    else:
                        cons[x] = 0

                means = np.zeros(len(chr_tile_df.index))
                variances = np.zeros(len(chr_tile_df.index))
                for i, index in enumerate(chr_tile_df.index):
                    print(i)
                    means[i], variances[i] = mean_var(index, variant)

                if show_strands:
                    marker = "${}$".format(strand)
                else:
                    marker = "${}$".format(chromosome)

                axs[0].plot(
                    cons, means, label="means", marker=marker, linestyle="None",
                )
                axs[1].plot(
                    cons, variances, label="variances", marker=marker, linestyle="None",
                )
        axs[0].set(
            title="{} ".format(variant) + plot_title,
            xlabel="mean number of consensus molecules per position",
            ylabel="mean",
            yscale="log",
        )
        axs[1].set(
            title="{} ".format(variant) + plot_title,
            xlabel="mean number of consensus molecules per position",
            ylabel="variance",
            yscale="log",
        )
        if save:
            file_name = "plots\\means_and_variances\\num_consensus\\{}".format(variant)
            if chromosome != "all":
                file_name += "_num_consensus_{}".format(chromosome)
            else:
                file_name += "_num_consensus"
            fig.savefig(file_name + ".png", dpi=600)
            fig.savefig(file_name + ".svg", dpi=1200)

        else:
            plt.show()


def plot_mean_var_against_num_tiles_in_exon(
    show_strands=True, chromosome="all", save=True
):
    """Plots the mean and variance against the number of tiles in the exon that tile is from."""
    for variant in VARIANTS:
        print(variant)
        fig, axs = plt.subplots(1, 2, figsize=(7, 4))

        plot_title = "num samples"
        if chromosome == "all":
            chrs_to_enumerate = chromosomes
        else:
            chrs_to_enumerate = chromosome
            plot_title += " {}".format(chromosome)

        for strand in ["+", "-"]:
            for chromosome in chrs_to_enumerate:
                chr_tile_df = tile_df.loc[
                    (tile_df["chromosome"] == chromosome)
                    & (tile_df["strand"] == strand)
                ]
                xs = chr_tile_df.index
                ys = np.zeros(len(xs))
                for x, i in zip(range(len(xs)), chr_tile_df.index):
                    ys[x] = len(exon_tiles_map[chr_tile_df.at[i, "exon"]].index)

                means = np.zeros(len(chr_tile_df.index))
                variances = np.zeros(len(chr_tile_df.index))
                for i, index in enumerate(chr_tile_df.index):
                    print(i)
                    means[i], variances[i] = mean_var(index, variant)

                if show_strands:
                    marker = "${}$".format(strand)
                else:
                    marker = "${}$".format(chromosome)

                axs[0].plot(
                    ys, means, label="means", marker=marker, linestyle="None",
                )
                axs[1].plot(
                    ys, variances, label="variances", marker=marker, linestyle="None",
                )
        axs[0].set(
            title="{} ".format(variant) + plot_title,
            xlabel="number of tiles in that exon",
            ylabel="mean",
            yscale="log",
        )
        axs[1].set(
            title="{} ".format(variant) + plot_title,
            xlabel="number of tiles in that exon",
            ylabel="variance",
            yscale="log",
        )
        if save:
            file_name = "plots\\means_and_variances\\num_tiles_in_exon\\{}".format(
                variant
            )
            if chromosome != "all":
                file_name += "_num_tiles_in_exon_{}".format(chromosome)
            else:
                file_name += "_num_tiles_in_exon"
            fig.savefig(file_name + ".png", dpi=600)
            fig.savefig(file_name + ".svg", dpi=1200)

        else:
            plt.show()


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


def plot_juicy_hist(juicy_tile_number=None, fit=None, bins_to_fit=-1):
    """Plots histograms of the juiciest positions and variants."""

    for i, juicy_row in juicy_df.iterrows():
        if juicy_tile_number is not None and i != juicy_tile_number:
            continue
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

        ax.set(title=plot_title, yscale="log")
        bottom, top = ax.get_ylim()  # to avoid axes being f*cked by the Poisson fit

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

        ax.set_ylim(bottom, top)  # to avoid axes beung fucked by the Poisson fit

        ax.plot([xs[0], xs[-1]], [1.0, 1.0], label="1/N")

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


def plot_found_hist(
    chromosome, tile_number, position, variant, fit=None, bins_to_fit=-1
):
    """Plots histograms of the positions I've found that I want to plot."""

    # for i, juicy_row in juicy_df.iterrows():
    #     juicy_data_df = tile_data_df(juicy_row["tile"])
    #     juicy_data_df = juicy_data_df.loc[
    #         (juicy_data_df["position"] == juicy_row["position"])
    #         & (juicy_data_df["variant"] == juicy_row["variant"])
    #     ]

    df = tile_data_df(tile_number)
    df = df.loc[(df["position"] == position) & (df["variant"] == variant)]
    print(df)
    plot_title = "chr {}, position {}, {}".format(chromosome, position, variant)

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
        ax.plot(xs, ys, marker="+", color="xkcd:piss yellow", label="Poisson fit mean")

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

    if fit == "beta-binomial fix mean" or fit == "both beta-binomial" or fit == "all":

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

    # def func(x):
    #     return betabinom.pmf(x, n, a, b) * N - 1

    # for x in xs:
    #     if func(x) <= 0:
    #         x
    #         break

    # print(xs)
    # print(func(xs))
    # print("first unlikely data: {}".format(x))
    plt.show()


def plot_hist_mean_by_tile(
    save=True,
    trim_and_flip=True,
    age="all",
    chromosome="all",
    xscale="log",
    yscale="log",
    context=False,
    cutoff_percentile=100,
    show_zeros=True,
):
    """Plots histogram of the means averaged over each tile."""
    if yscale == "log":
        log = True
    else:
        log = False
    if context:
        contexts_to_iterate = CONTEXTS
    else:
        contexts_to_iterate = ["no context"]

    for context_string in contexts_to_iterate:

        for variant in VARIANTS:
            print(variant)

            fig, axs = plt.subplots(3, figsize=(8, 8))

            if chromosome == "all":
                chrs_to_use = chromosomes
            else:
                chrs_to_use = [chromosome]

            for strand in ["+", "-"]:

                strand_tile_df = tile_df.loc[
                    (tile_df["strand"] == strand)
                    & (tile_df["chromosome"].isin(chrs_to_use))
                ]

                means = np.zeros(len(strand_tile_df.index))
                variances = np.zeros(len(strand_tile_df.index))
                for i, index in enumerate(strand_tile_df.index):
                    means[i], variances[i] = mean_var(
                        index, variant, trim_and_flip, age
                    )

                if strand == "+":
                    pos_means = means
                    pos_vars = variances
                    pos_max_mean = np.max(pos_means)
                    pos_max_var = np.max(pos_vars)
                else:
                    neg_means = means
                    neg_vars = variances
                    neg_max_mean = np.max(neg_means)
                    neg_max_var = np.max(neg_vars)

                    max_mean = np.percentile(
                        [pos_max_mean, neg_max_mean], cutoff_percentile
                    )
                    max_var = np.percentile(
                        [pos_max_var, neg_max_var], cutoff_percentile
                    )

            # throw away all elements for mean=0 so I can calculate D = σ²/μ
            pos_non_zero_means = pos_means[pos_means.nonzero()]
            pos_non_zero_vars = pos_vars[pos_means.nonzero()]
            pos_Ds = pos_non_zero_vars / pos_non_zero_means

            neg_non_zero_means = neg_means[neg_means.nonzero()]
            neg_non_zero_vars = neg_vars[neg_means.nonzero()]
            neg_Ds = neg_non_zero_vars / neg_non_zero_means

            n_bins = 100

            if xscale == "log":
                mean_bins = np.logspace(-4, 4, n_bins + 1)
                var_bins = np.logspace(-4, 4, n_bins + 1)
            else:
                mean_bins = np.linspace(0, max_mean * (n_bins + 1) / n_bins, n_bins + 1)
                var_bins = np.linspace(0, max_var * (n_bins + 1) / n_bins, n_bins + 1)

            if show_zeros:

                def f(x):
                    """Turns 0 to 1e-4."""
                    if x == 0:
                        return 1e-4
                    return x

                f = np.vectorize(f)
                pos_means = f(pos_means)
                neg_means = f(neg_means)
                pos_vars = f(pos_vars)
                neg_vars = f(neg_vars)

            axs[0].hist(
                [pos_means, neg_means],
                bins=mean_bins,
                stacked=True,
                log=log,
                color=["blue", "orange"],
                label=["+ means", "- means"],
            )

            axs[1].hist(
                [pos_vars, neg_vars],
                bins=var_bins,
                stacked=True,
                log=log,
                color=["blue", "orange"],
                label=["+ variances", "- variances"],
            )
            axs[2].hist(
                [pos_Ds, neg_Ds],
                bins=int(np.amax(np.append(pos_Ds, neg_Ds))),
                stacked=True,
                log=log,
                color=["blue", "orange"],
                label=["+ Ds", "- Ds"],
            )

            axs0_title = "{} means".format(variant)
            if context:
                axs0_title += " {}".format(context_string)
            if trim_and_flip:
                axs0_title += " t&f"
            if age != "all":
                axs0_title += " age {}".format(age)
            if chromosome != "all":
                axs0_title += " {}".format(chromosome)

            axs[0].set(
                title=axs0_title, xlabel=r"$\mu$", ylabel="frequency", xscale=xscale
            )
            axs[1].set(
                title="variances",
                xlabel=r"$\sigma^2$",
                ylabel="frequency",
                xscale=xscale,
            )
            axs[2].set(
                title=r"$D = \frac{\sigma^2}{\mu}$",
                xlabel=r"$D$",
                ylabel="frequency",
                xscale=xscale,
            )
            for ax in axs:
                ax.set_ylim(bottom=10 ** -0.5)
                ax.legend()

            fig.tight_layout()
            if save:
                location = "plots\\means_and_variances\\histograms_tile\\"

                file_name = "{}_mean_var_tile_histogram".format(variant)
                if context:
                    file_name += "_{}".format(context_string)
                if trim_and_flip:
                    file_name += "_t&f"
                if age != "all":
                    file_name += "_age_{}".format(age)
                if chromosome != "all":
                    file_name += "_{}".format(chromosome)
                fig.savefig(location + file_name + ".png", dpi=600)
                fig.savefig(location + file_name + ".svg", dpi=1200)
            else:
                plt.show()
            plt.close("all")


def plot_hist_mean_by_position(
    save=True,
    trim_and_flip=True,
    age="all",
    chromosome="all",
    xscale="log",
    yscale="log",
    context=False,
    cutoff_percentile=100,
    show_zeros=True,
):
    """Plots histogram of the means for each position."""
    if yscale == "log":
        log = True
    else:
        log = False
    if context:
        contexts_to_iterate = CONTEXTS
    else:
        contexts_to_iterate = ["no context"]

    for context_string in contexts_to_iterate:

        for variant in VARIANTS:
            print(variant)

            fig, axs = plt.subplots(3, figsize=(8, 8))

            if chromosome == "all":
                chrs_to_use = chromosomes
            else:
                chrs_to_use = [chromosome]

            for strand in ["+", "-"]:

                strand_tile_df = tile_df.loc[
                    (tile_df["strand"] == strand)
                    & (tile_df["chromosome"].isin(chrs_to_use))
                ]

                means = []
                variances = []
                for tile_number in strand_tile_df.index:
                    data_df = tile_data_df(tile_number)
                    data_df = data_df.loc[data_df["variant"] == variant]
                    positions = data_df["position"].unique()
                    tile_means = np.zeros(len(positions))
                    tile_vars = np.zeros(len(positions))
                    for i, position in enumerate(positions):
                        df = data_df.loc[data_df["position"] == position]
                        tile_means[i] = np.mean(df["downsample"])
                        tile_vars[i] = np.var(df["downsample"])
                    means.append(tile_means)
                    variances.append(tile_vars)
                means = np.hstack(means)
                variances = np.hstack(variances)

                if strand == "+":
                    pos_means = means
                    pos_vars = variances
                    pos_max_mean = np.max(pos_means)
                    pos_max_var = np.max(pos_vars)
                else:
                    neg_means = means
                    neg_vars = variances
                    neg_max_mean = np.max(neg_means)
                    neg_max_var = np.max(neg_vars)

                    max_mean = np.percentile(
                        [pos_max_mean, neg_max_mean], cutoff_percentile
                    )
                    max_var = np.percentile(
                        [pos_max_var, neg_max_var], cutoff_percentile
                    )
                    print(max_mean)
                    print(max_var)

            # throw away all elements for mean=0 so I can calculate D = σ²/μ
            pos_non_zero_means = pos_means[pos_means.nonzero()]
            pos_non_zero_vars = pos_vars[pos_means.nonzero()]
            pos_Ds = pos_non_zero_vars / pos_non_zero_means

            neg_non_zero_means = neg_means[neg_means.nonzero()]
            neg_non_zero_vars = neg_vars[neg_means.nonzero()]
            neg_Ds = neg_non_zero_vars / neg_non_zero_means

            n_bins = 100

            if xscale == "log":
                mean_bins = np.logspace(-4, 4, n_bins + 1)
                var_bins = np.logspace(-4, 4, n_bins + 1)
            else:
                mean_bins = np.linspace(0, max_mean * (n_bins + 1) / n_bins, n_bins + 1)
                var_bins = np.linspace(0, max_var * (n_bins + 1) / n_bins, n_bins + 1)

            if show_zeros:

                def f(x):
                    """Turns 0 to 1e-4."""
                    if x == 0:
                        return 10e-4
                    return x

                f = np.vectorize(f)
                pos_means = f(pos_means)
                neg_means = f(neg_means)
                pos_vars = f(pos_vars)
                neg_vars = f(neg_vars)

            axs[0].hist(
                [pos_means, neg_means],
                bins=mean_bins,
                stacked=True,
                log=log,
                color=["blue", "orange"],
                label=["+ means", "- means"],
            )

            axs[1].hist(
                [pos_vars, neg_vars],
                bins=var_bins,
                stacked=True,
                log=log,
                color=["blue", "orange"],
                label=["+ variances", "- variances"],
            )
            axs[2].hist(
                [pos_Ds, neg_Ds],
                bins=int(np.amax(np.append(pos_Ds, neg_Ds))),
                stacked=True,
                log=log,
                color=["blue", "orange"],
                label=["+ Ds", "- Ds"],
            )

            axs0_title = "{} means".format(variant)
            if context:
                axs0_title += " {}".format(context_string)
            if trim_and_flip:
                axs0_title += " t&f"
            if age != "all":
                axs0_title += " age {}".format(age)
            if chromosome != "all":
                axs0_title += " {}".format(chromosome)

            axs[0].set(
                title=axs0_title, xlabel=r"$\mu$", ylabel="frequency", xscale=xscale
            )
            axs[1].set(
                title="variances",
                xlabel=r"$\sigma^2$",
                ylabel="frequency",
                xscale=xscale,
            )
            axs[2].set(
                title=r"$D = \frac{\sigma^2}{\mu}$",
                xlabel=r"$D$",
                ylabel="frequency",
                xscale=xscale,
            )
            for ax in axs:
                ax.set_ylim(bottom=10 ** -0.5)
                ax.legend()

            fig.tight_layout()
            if save:
                location = "plots\\means_and_variances\\histograms_position\\"

                file_name = "{}_mean_var_position_histogram".format(variant)
                if context:
                    file_name += "_{}".format(context_string)
                if trim_and_flip:
                    file_name += "_t&f"
                if age != "all":
                    file_name += "_age_{}".format(age)
                if chromosome != "all":
                    file_name += "_{}".format(chromosome)
                fig.savefig(location + file_name + ".png", dpi=600)
                fig.savefig(location + file_name + ".svg", dpi=1200)
            else:
                plt.show()
            plt.close("all")


def plot_hist_D_by_position(
    save=True,
    trim_and_flip=True,
    age="all",
    chromosome="all",
    xscale="log",
    yscale="log",
    cutoff_percentile=100,
    show_zeros=True,
):
    """Plots histogram of the means for each position."""
    if yscale == "log":
        log = True
    else:
        log = False

    fig, axs = plt.subplots(4, 3, figsize=(16, 20))

    for variant, ax in zip(VARIANTS, axs.flatten()):
        print(variant)

        if chromosome == "all":
            chrs_to_use = chromosomes
        else:
            chrs_to_use = [chromosome]

        for strand in ["+", "-"]:

            strand_tile_df = tile_df.loc[
                (tile_df["strand"] == strand)
                & (tile_df["chromosome"].isin(chrs_to_use))
            ]

            means = []
            variances = []
            for tile_number in strand_tile_df.index:
                data_df = tile_data_df(tile_number)
                data_df = data_df.loc[data_df["variant"] == variant]
                positions = data_df["position"].unique()
                tile_means = np.zeros(len(positions))
                tile_vars = np.zeros(len(positions))
                for i, position in enumerate(positions):
                    df = data_df.loc[data_df["position"] == position]
                    tile_means[i] = np.mean(df["downsample"])
                    tile_vars[i] = np.var(df["downsample"])
                means.append(tile_means)
                variances.append(tile_vars)
            means = np.hstack(means)
            variances = np.hstack(variances)

            if strand == "+":
                pos_means = means
                pos_vars = variances
                pos_max_mean = np.max(pos_means)
                pos_max_var = np.max(pos_vars)
            else:
                neg_means = means
                neg_vars = variances
                neg_max_mean = np.max(neg_means)
                neg_max_var = np.max(neg_vars)

                max_mean = np.percentile(
                    [pos_max_mean, neg_max_mean], cutoff_percentile
                )
                max_var = np.percentile([pos_max_var, neg_max_var], cutoff_percentile)
                print(max_mean)
                print(max_var)

        # throw away all elements for mean=0 so I can calculate D = σ²/μ
        pos_non_zero_means = pos_means[pos_means.nonzero()]
        pos_non_zero_vars = pos_vars[pos_means.nonzero()]
        pos_Ds = pos_non_zero_vars / pos_non_zero_means

        neg_non_zero_means = neg_means[neg_means.nonzero()]
        neg_non_zero_vars = neg_vars[neg_means.nonzero()]
        neg_Ds = neg_non_zero_vars / neg_non_zero_means

        n_bins = 100

        Ds = np.append(pos_Ds, neg_Ds)
        bins = np.logspace(np.log10(np.amin(Ds[Ds > 0])), np.log10(np.amax(Ds)), n_bins)
        print("min D: {}".format(np.amin(Ds)))
        n, bins, patches = ax.hist(
            [pos_Ds, neg_Ds],
            bins=bins,
            stacked=True,
            log=log,
            color=["blue", "orange"],
            label=["+ Ds", "- Ds"],
        )
        print(bins)

        ax_title = "{} means".format(variant)
        if trim_and_flip:
            ax_title += " t&f"
        if age != "all":
            ax_title += " age {}".format(age)
        if chromosome != "all":
            ax_title += " {}".format(chromosome)

        ax.set(
            title=variant,
            xlabel=r"$D = \frac{\sigma ^2}{\mu}$",
            ylabel="frequency",
            xscale=xscale,
        )
        ax.set_ylim(bottom=10 ** -0.5)
        ax.legend()

        fig.tight_layout()
    if save:
        location = "plots\\means_and_variances\\index_of_dispersion_hist\\"

        file_name = "D_histogram"

        if trim_and_flip:
            file_name += "_t&f"
        if age != "all":
            file_name += "_age_{}".format(age)
        if chromosome != "all":
            file_name += "_{}".format(chromosome)
        fig.savefig(location + file_name + ".png", dpi=600)
        fig.savefig(location + file_name + ".svg", dpi=1200)
    else:
        plt.show()
    plt.close("all")
