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
    N_0 = percentile(q)
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


def sort_caroline_seqs():
    """
    Orders data from Caroline by chromosome.
    """
    df = pd.read_csv(
        file_names["Caroline seqs"],
        header=None,
        names=["chromosome", "start", "end", "who tf knows", "length", "strand"],
        sep="\t",
    )
    df["chromosome"] = df["chromosome"].str[3:]
    df = df.sort_values(
        by="chromosome", kind="mergesort", key=vec_sorter, ignore_index=True,
    )  # use mergesort for stability

    df.to_csv(file_names["Caroline seqs sorted"], sep="\t")


def separating_seqs(seq_numbers):
    """
    Separates the seqs in full_data.txt by the seqs in Caroline's file.

    seq_numbers = list of seqs to produce files for.

    Note: cannot do all seqs at once.
    """
    reduced_seq_df = seq_df.loc[seq_numbers]
    seq_dfs = []
    for i in range(len(seq_numbers)):
        seq_dfs.append(pd.DataFrame(columns=sample_column_names))

    for j, chunk in enumerate(
        pd.read_csv(file_names["downsampled data"], chunksize=CHUNKSIZE, index_col=0)
    ):
        print("chunk {}".format(j))
        for i, (k, seq) in zip(range(seq_numbers.size), reduced_seq_df.iterrows()):
            # print("seq {}".format(i))

            seq_dfs[i] = seq_dfs[i].append(
                chunk.loc[
                    (chunk["position"] >= seq["start"])
                    & (chunk["position"] <= seq["end"])
                ],
                ignore_index=True,
            )
        gc.collect()

    for i, df in zip(seq_numbers, seq_dfs):
        df.to_csv(file_names["seq"].format(i))


aggregation_functions = {
    "num subs": "sum",
    "num consensus molecules": "sum",
    "downsample": "sum",
}


def seq_data_df(seq_number, group_by=None, trim_and_flip=True):
    """
    Returns DataFrame from seq_(seq_number).csv.

    group_by = "position" combines samples at the same position.
    group_by = "ID" combines samples with the same ID (person and age)
    Both grouping options use the trimmed and flipped data by default.
    """
    if trim_and_flip:
        if group_by == "ID":
            return pd.read_csv(file_names["seq group IDs t&f"].format(seq_number))
        elif group_by == "position":
            return pd.read_csv(file_names["seq group positions t&f"].format(seq_number))
        else:
            return pd.read_csv(file_names["seq t&f"].format(seq_number), index_col=0,)
    else:
        if group_by == "ID":
            return pd.read_csv(file_names["seq group IDs"].format(seq_number))
        elif group_by == "position":
            return pd.read_csv(file_names["seq group positions"].format(seq_number))
        else:
            return pd.read_csv(file_names["seq"].format(seq_number), index_col=0)


def group_by_position(seq_number, trim_and_flip=True):
    """
    Groups the specified seq by position.
    """
    df = seq_data_df(seq_number, trim_and_flip=trim_and_flip)
    df = df.drop(columns=["sample ID"])
    df = df.groupby(["position", "chromosome", "sub"]).agg(aggregation_functions)
    if trim_and_flip:
        df.to_csv(file_names["seq group positions t&f"].format(seq_number))
    else:
        df.to_csv(file_names["seq group positions"].format(seq_number))


def group_by_ID(seq_number, trim_and_flip=True):
    """
    Groups the specified seq by ID (person and age).
    """
    df = seq_data_df(seq_number, trim_and_flip=trim_and_flip)
    df = df.drop(columns=["position"])
    df = df.groupby(["sample ID", "chromosome", "sub"]).agg(aggregation_functions)
    if trim_and_flip:
        df.to_csv(file_names["seq group IDs t&f"].format(seq_number))
    else:
        df.to_csv(file_names["seq group IDs"].format(seq_number))


def trim_and_flip(exon_number):
    """
    Flips neg data to be what the actual change was.

    The called changes are all relative to the top strand. This function creates files which identify the actual sub seen.
    """
    seqs = exon_seqs_map[exon_number]
    next_df = seq_data_df(seqs.index[0], trim_and_flip=False)
    if len(seqs.index) >= 2:
        for i in seqs.index[:-1]:
            df = next_df
            next_df = seq_data_df(i + 1, trim_and_flip=False)
            cond = ~df["position"].isin(next_df["position"])
            next_cond = ~next_df["position"].isin(df["position"])
            df = df.loc[cond, :]
            next_df = next_df.loc[next_cond, :]
            if seq_df.at[i, "strand"] == "-":
                df = df.replace({"sub": sub_complement_map})
            df.to_csv(file_names["seq t&f"].format(i))
    next_df.to_csv(file_names["seq t&f"].format(seqs.index[-1]))


### Wrappers ###


def group_by_position_wrapper(trim_and_flip=True):
    """
    Repeats group by position for all seqs.
    """
    for i in np.arange(0, 1063):
        group_by_position(i, trim_and_flip)
        print(i)


def group_by_ID_wrapper(trim_and_flip=True):
    """
    Repeats group by ID for all seqs.
    """
    for i in np.arange(0, 1063):
        group_by_ID(i, trim_and_flip)
        print(i)


def separating_seqs_wrapper():
    """
    Separates seqs in chunks of 100 to avoid memory issues.
    """
    print("Separating seqs 0 to 99")
    separating_seqs(np.arange(0, 100))
    print("Separating seqs 100 to 199")
    separating_seqs(np.arange(100, 200))
    print("Separating seqs 200 to 299")
    separating_seqs(np.arange(200, 300))
    print("Separating seqs 300 to 399")
    separating_seqs(np.arange(300, 400))
    print("Separating seqs 400 to 499")
    separating_seqs(np.arange(400, 500))
    print("Separating seqs 500 to 599")
    separating_seqs(np.arange(500, 600))
    print("Separating seqs 600 to 699")
    separating_seqs(np.arange(600, 700))
    print("Separating seqs 700 to 799")
    separating_seqs(np.arange(700, 800))
    print("Separating seqs 800 to 899")
    separating_seqs(np.arange(800, 900))
    print("Separating seqs 900 to 999")
    separating_seqs(np.arange(900, 1000))
    print("Separating seqs 1000 to 1062")
    separating_seqs(np.arange(1000, 1063))


def group_strands_wrapper():
    """
    Repeats group strands for all exons
    """
    for i in np.arange(len(exon_df.index)):
        print(i)
        group_strands(i)


def trim_and_flip_wrapper():
    for i in exon_df.index:
        trim_and_flip(i)
        print("Gene {}".format(i))


def refresh_data(redownsample=False):
    """
    Runs all the functions in turn to freshen up the data
    """
    if redownsample:
        downsample()
        separating_seqs_wrapper()
    group_by_ID_wrapper(trim_and_flip=False)
    group_by_position_wrapper(trim_and_flip=False)
    trim_and_flip_wrapper()
    group_by_ID_wrapper()
    group_by_position_wrapper()


### Dataframes ###


# Dataframe of exons
exon_df = pd.read_csv(
    file_names["Wing exons"],
    header=None,
    names=["chromosome", "start", "end"],
    sep="\t",
)

# DataFrame of seq information
seq_df = pd.read_csv(file_names["Caroline seqs sorted"], sep="\t", index_col=0)

exon_seqs_map = {}
for i, exon in exon_df.iterrows():
    exon_seqs_map[i] = seq_df.loc[
        (seq_df["start"] >= exon["start"]) & (seq_df["end"] <= exon["end"])
    ]


### Miscellaneous ###


def empty_seqs():
    """
    Identifies any empty seqs.
    """
    for i, exon in exon_df.iterrows():
        for j, seq in exon_seqs_map[i].iterrows():
            if 0 == len(seq_data_df(j).index):
                print("Empty seq found, exon {} seq {}".format(i, j))


def seqs_not_in_exons():
    """
    Identifies any seqs that don't appear in any exons
    """
    seqs_in_exons = []
    for i in exon_df.index:
        seqs_in_exons.extend(exon_seqs_map[i].index)

    for i in seq_df.index:
        if ~np.isin(i, seqs_in_exons):
            print("Sequence {} is not in any exon".format(i))

    print(seqs_in_exons)
