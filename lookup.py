"""
# 2020 Blundell lab internship

Contains functions for accessing the ALSPAC data using Pandas dataframes.
"""
import numpy as np
import pandas as pd
import os
from constants import CHANGES, CHUNKSIZE, LANES, PEOPLE


sample_column_names = [
    "chromosome",
    "position",
    "sub",
    "num subs",
    "num consensus molecules",
    "sample ID",
]


id_df = pd.read_csv(
    "data_files\\id.txt",
    header=None,
    names=["lane1", "lane2", "lane3", "lane4", "ID"],
    index_col=-1,
)

gene_df = pd.read_csv(
    "data_files\\Wing_genes.bed",
    header=None,
    names=["chromosome", "start", "end"],
    sep="\t",
)

seq_df = pd.read_csv(
    "data_files\\Caroline_sequences.bed",
    header=None,
    names=["chromosome", "start", "end", "who tf knows", "length", "strand"],
    sep="\t",
)

seq_df.sort_values(
    "chromosome", inplace=True, kind="mergesort", ignore_index=True
)  # use mergesort for stability

gene_seqs_map = {}
for i, gene in gene_df.iterrows():
    gene_seqs_map[i] = seq_df.loc[
        (seq_df["start"] >= gene["start"]) & (seq_df["end"] <= gene["end"])
    ]


def seq_data_df(sequence_number, group_by=None):
    """
    Returns DataFrame from seq_(sequence_number).csv.

    group_by = "position" combines samples at the same position.
    group_by = "ID" combines samples with the same ID (person and age)
    """
    if group_by == "position":
        return pd.read_csv(
            "data_files\\sequences_by_position\\seq_{}_group_positions.csv".format(
                sequence_number
            )
        )
    elif group_by == "ID":
        return pd.read_csv(
            "data_files_sequences_by_ID\\seq_{}_group_ID.csv".format(sequence_number)
        )
    else:
        return pd.read_csv(
            "data_files\\sequences\\seq_{}.csv".format(sequence_number), index_col=0
        )


def separating_sequences(sequence_numbers):
    """
    Separates the sequences in full_data.txt by the sequences in Caroline's file.

    sequence_numbers = list of sequences to produce files for.

    Note: cannot do all sequences at once.
    """
    reduced_seq_df = seq_df.loc[sequence_numbers]
    seq_dfs = []
    for i in range(len(sequence_numbers)):
        seq_dfs.append(pd.DataFrame(columns=sample_column_names))

    for j, chunk in enumerate(
        pd.read_csv(
            "data_files\\downsampled_data.txt", chunksize=CHUNKSIZE, index_col=0
        )
    ):
        print("chunk {}".format(j))
        for i, (k, sequence) in zip(
            range(sequence_numbers.size), reduced_seq_df.iterrows()
        ):
            # print("sequence {}".format(i))

            seq_dfs[i] = seq_dfs[i].append(
                chunk.loc[
                    (chunk["position"] >= sequence["start"])
                    & (chunk["position"] <= sequence["end"])
                ],
                ignore_index=True,
            )

    for i, df in zip(sequence_numbers, seq_dfs):
        df.to_csv("data_files\\sequences\\seq_{}.csv".format(i))


def percentile(q):
    """
    Finds the qth percentile of the number of consensus molecules.
    """
    Nss = []
    for i, chunk in enumerate(
        pd.read_csv(
            "data_files\\full_data.txt",
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
    Nss = np.hstack(Nss)
    return np.percentile(Nss, q)


def downsample(q):
    """
    Downsamples full_data.txt to qth percentile of number of consensus molecules.
    """
    if os.path.isfile("data_files\\downsampled_data.txt"):
        os.remove("data_files\\downsampled_data.txt")
    rng = np.random.default_rng()
    header = True
    N_0 = percentile(q)
    for i, chunk in enumerate(
        pd.read_csv(
            "data_files\\full_data.txt",
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
        chunk.to_csv("data_files\\downsampled_data.txt", mode="a", header=header)
        header = False


aggregation_functions = {
    "num subs": "sum",
    "num consensus molecules": "sum",
    "downsample": "sum",
}


def group_by_position(sequence_number):
    """
    Groups the specified sequence by position.
    """
    df = seq_data_df(sequence_number)
    df = df.drop(columns=["sample ID"])
    df = df.groupby(["position", "chromosome", "sub"]).agg(aggregation_functions)
    df.to_csv(
        "data_files\\sequences_by_position\\seq_{}_group_positions.csv".format(
            sequence_number
        )
    )


def group_by_position_wrapper():
    """
    Repeats group by position for all sequences.
    """
    for i in np.arange(0, 1063):
        group_by_position(i)
        print(i)


def group_by_ID(sequence_number):
    """
    Groups the specified sequence by ID (person and age).
    """
    df = seq_data_df(sequence_number)
    df = df.drop(columns=["position"])
    df = df.groupby(["sample ID", "chromosome", "sub"]).agg(aggregation_functions)
    df.to_csv(
        "data_files\\sequences_by_ID\\seq_{}_group_ID.csv".format(sequence_number)
    )


def group_by_ID_wrapper():
    """
    Repeats group by ID for all sequences.
    """
    for i in np.arange(0, 1063):
        group_by_ID(i)
        print(i)


def separating_sequences_wrapper():
    """
    Separates sequences in chunks of 100 to avoid memory issues.
    """
    separating_sequences(np.arange(0, 100))
    separating_sequences(np.arange(100, 200))
    separating_sequences(np.arange(200, 300))
    separating_sequences(np.arange(300, 400))
    separating_sequences(np.arange(400, 500))
    separating_sequences(np.arange(500, 600))
    separating_sequences(np.arange(600, 700))
    separating_sequences(np.arange(700, 800))
    separating_sequences(np.arange(800, 900))
    separating_sequences(np.arange(900, 1000))
    separating_sequences(np.arange(1000, 1063))


def empty_sequences():
    """
    Identifies any empty sequences.
    """
    for i, gene in gene_df.iterrows():
        for j, seq in gene_seqs_map[i].iterrows():
            if 0 == len(seq_data_df(j).index):
                print("Empty sequence found, gene {} sequence {}".format(i, j))
