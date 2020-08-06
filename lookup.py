"""
# 2020 Blundell lab internship

Contains functions for accessing the ALSPAC data using Pandas dataframes.
"""
import numpy as np
import pandas as pd
import os
import gc
from constants import (
    SUBS,
    CHUNKSIZE,
    file_names,
    LANES,
    PEOPLE,
    sub_complement_map,
    vec_sorter,
)


sample_column_names = [
    "chromosome",
    "position",
    "sub",
    "num subs",
    "num consensus molecules",
    "sample ID",
]

### Preparing files ###


def percentile(q):
    """
    Finds the qth percentile of the number of consensus molecules.
    """
    Nss = []
    for i, chunk in enumerate(
        pd.read_csv(
            file_names["data"],
            chunksize=CHUNKSIZE,
            header=None,
            names=sample_column_names,
            sep="\t",
            usecols=["num consensus molecules"],
        )
    ):
        print(i)
        Ns = chunk["num consensus molecules"].to_numpy()
        Ns = Ns.astype(int)
        Nss.append(Ns)
        gc.collect()
    Nss = np.hstack(Nss)
    return np.percentile(Nss, q)


def downsample(q=50):
    """
    Downsamples full_data.txt to qth percentile of number of consensus molecules.

    Defaults to 50th percentile
    """
    if os.path.isfile(file_names["downsampled data"]):
        os.remove(file_names["downsampled data"])
    rng = np.random.default_rng()
    header = True
    N_0 = 5000  # percentile(q)
    print("Downsampled to {} consensus molecules".format(N_0))
    for i, chunk in enumerate(
        pd.read_csv(
            file_names["data"],
            chunksize=CHUNKSIZE,
            header=None,
            names=sample_column_names,
            sep="\t",
        )
    ):
        print(i)
        chunk["downsample"] = rng.binomial(
            n=N_0, p=chunk["num subs"] / chunk["num consensus molecules"]
        )
        chunk.to_csv(file_names["downsampled data"], mode="a", header=header)
        header = False
        gc.collect()


def sort_caroline_exons():
    """
    Orders data from Caroline by chromosome.
    """
    df = pd.read_csv(
        file_names["Caroline exons"],
        header=None,
        names=["chromosome", "start", "end", "who tf knows", "length", "strand"],
        sep="\t",
    )
    df["chromosome"] = df["chromosome"].str[3:]
    df = df.sort_values(
        by="chromosome", kind="mergesort", key=vec_sorter, ignore_index=True,
    )  # use mergesort for stability

    df.to_csv(file_names["Caroline exons sorted"], sep="\t")


def separating_exons(exon_numbers):
    """
    Separates the exons in full_data.txt by the exons in Caroline's file.

    exon_numbers = list of exons to produce files for.

    Note: cannot do all exons at once.
    """
    reduced_exon_df = exon_df.loc[exon_numbers]
    exon_dfs = []
    for i in range(len(exon_numbers)):
        exon_dfs.append(pd.DataFrame(columns=sample_column_names))

    for j, chunk in enumerate(
        pd.read_csv(file_names["downsampled data"], chunksize=CHUNKSIZE, index_col=0)
    ):
        print("chunk {}".format(j))
        for i, (k, exon) in zip(range(exon_numbers.size), reduced_exon_df.iterrows()):
            # print("exon {}".format(i))

            exon_dfs[i] = exon_dfs[i].append(
                chunk.loc[
                    (chunk["position"] >= exon["start"])
                    & (chunk["position"] <= exon["end"])
                ],
                ignore_index=True,
            )
        gc.collect()

    for i, df in zip(exon_numbers, exon_dfs):
        df.to_csv(file_names["exon"].format(i))


aggregation_functions = {
    "num subs": "sum",
    "num consensus molecules": "sum",
    "downsample": "sum",
}


def exon_data_df(exon_number, group_by=None, trim_and_flip=True):
    """
    Returns DataFrame from exon_(exon_number).csv.

    group_by = "position" combines samples at the same position.
    group_by = "ID" combines samples with the same ID (person and age)
    Both grouping options use the trimmed and flipped data by default.
    """
    if trim_and_flip:
        if group_by == "ID":
            return pd.read_csv(file_names["exon group IDs t&f"].format(exon_number))
        elif group_by == "position":
            return pd.read_csv(
                file_names["exon group positions t&f"].format(exon_number)
            )
        else:
            return pd.read_csv(file_names["exon t&f"].format(exon_number), index_col=0,)
    else:
        if group_by == "ID":
            return pd.read_csv(file_names["exon group IDs"].format(exon_number))
        elif group_by == "position":
            return pd.read_csv(file_names["exon group positions"].format(exon_number))
        else:
            return pd.read_csv(file_names["exon"].format(exon_number), index_col=0)


def group_by_position(exon_number, trim_and_flip=True):
    """
    Groups the specified exon by position.
    """
    df = exon_data_df(exon_number, trim_and_flip=trim_and_flip)
    df = df.drop(columns=["sample ID"])
    df = df.groupby(["position", "chromosome", "sub"]).agg(aggregation_functions)
    if trim_and_flip:
        df.to_csv(file_names["exon group positions t&f"].format(exon_number))
    else:
        df.to_csv(file_names["exon group positions"].format(exon_number))


def group_by_ID(exon_number, trim_and_flip=True):
    """
    Groups the specified exon by ID (person and age).
    """
    df = exon_data_df(exon_number, trim_and_flip=trim_and_flip)
    df = df.drop(columns=["position"])
    df = df.groupby(["sample ID", "chromosome", "sub"]).agg(aggregation_functions)
    if trim_and_flip:
        df.to_csv(file_names["exon group IDs t&f"].format(exon_number))
    else:
        df.to_csv(file_names["exon group IDs"].format(exon_number))


def trim_and_flip(gene_number):
    """
    Flips neg data to be what the actual change was.

    The called changes are all relative to the top strand. This function creates files which identify the actual sub seen.
    """
    exons = gene_exons_map[gene_number]
    next_df = exon_data_df(exons.index[0], trim_and_flip=False)
    if len(exons.index) >= 2:
        for i in exons.index[:-1]:
            df = next_df
            next_df = exon_data_df(i + 1, trim_and_flip=False)
            cond = ~df["position"].isin(next_df["position"])
            next_cond = ~next_df["position"].isin(df["position"])
            df = df.loc[cond, :]
            next_df = next_df.loc[next_cond, :]
            if exon_df.at[i, "strand"] == "-":
                df = df.replace({"sub": sub_complement_map})
            df.to_csv(file_names["exon t&f"].format(i))
    next_df.to_csv(file_names["exon t&f"].format(exons.index[-1]))


### Wrappers ###


def group_by_position_wrapper(trim_and_flip=True):
    """
    Repeats group by position for all exons.
    """
    for i in np.arange(0, 1063):
        group_by_position(i, trim_and_flip)
        print(i)


def group_by_ID_wrapper(trim_and_flip=True):
    """
    Repeats group by ID for all exons.
    """
    for i in np.arange(0, 1063):
        group_by_ID(i, trim_and_flip)
        print(i)


def separating_exons_wrapper():
    """
    Separates exons in chunks of 100 to avoid memory issues.
    """
    print("Separating exons 0 to 99")
    separating_exons(np.arange(0, 100))
    print("Separating exons 100 to 199")
    separating_exons(np.arange(100, 200))
    print("Separating exons 200 to 299")
    separating_exons(np.arange(200, 300))
    print("Separating exons 300 to 399")
    separating_exons(np.arange(300, 400))
    print("Separating exons 400 to 499")
    separating_exons(np.arange(400, 500))
    print("Separating exons 500 to 599")
    separating_exons(np.arange(500, 600))
    print("Separating exons 600 to 699")
    separating_exons(np.arange(600, 700))
    print("Separating exons 700 to 799")
    separating_exons(np.arange(700, 800))
    print("Separating exons 800 to 899")
    separating_exons(np.arange(800, 900))
    print("Separating exons 900 to 999")
    separating_exons(np.arange(900, 1000))
    print("Separating exons 1000 to 1062")
    separating_exons(np.arange(1000, 1063))


def group_strands_wrapper():
    """
    Repeats group strands for all genes
    """
    for i in np.arange(len(gene_df.index)):
        print(i)
        group_strands(i)


def trim_and_flip_wrapper():
    for i in gene_df.index:
        trim_and_flip(i)
        print("Gene {}".format(i))


def refresh_data(redownsample=False):
    """
    Runs all the functions in turn to freshen up the data
    """
    if redownsample:
        downsample()
        separating_exons_wrapper()
    group_by_ID_wrapper(trim_and_flip=False)
    group_by_position_wrapper(trim_and_flip=False)
    trim_and_flip_wrapper()
    group_by_ID_wrapper()
    group_by_position_wrapper()


### Dataframes ###


# Dataframe of genes
gene_df = pd.read_csv(
    file_names["Wing genes"],
    header=None,
    names=["chromosome", "start", "end"],
    sep="\t",
)

# DataFrame of exon information
exon_df = pd.read_csv(file_names["Caroline exons sorted"], sep="\t", index_col=0)

gene_exons_map = {}
for i, gene in gene_df.iterrows():
    gene_exons_map[i] = exon_df.loc[
        (exon_df["start"] >= gene["start"]) & (exon_df["end"] <= gene["end"])
    ]


### Miscellaneous ###


def empty_exons():
    """
    Identifies any empty exons.
    """
    for i, gene in gene_df.iterrows():
        for j, exon in gene_exons_map[i].iterrows():
            if 0 == len(exon_data_df(j).index):
                print("Empty exon found, gene {} exon {}".format(i, j))


def exons_not_in_genes():
    """
    Identifies any exons that don't appear in any genes
    """
    exons_in_genes = []
    for i in gene_df.index:
        exons_in_genes.extend(gene_exons_map[i].index)

    for i in exon_df.index:
        if ~np.isin(i, exons_in_genes):
            print("Sequence {} is not in any gene".format(i))

    print(exons_in_genes)
