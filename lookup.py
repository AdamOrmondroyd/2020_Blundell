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
    # print("{}: {}".format(i, len(gene_seqs_map[i].index)))


def seq_data_df(sequence_number, group_by=None):
    """
    Returns DataFrame from seq_(sequence_number).csv

    group_by = "position" combines samples at the same position
    """
    if group_by == "position":
        return pd.read_csv(
            "data_files\\sequences_by_position\\seq_{}_group_positions.csv".format(
                sequence_number
            )
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


def figuring_out_the_data():
    seq_df = pd.DataFrame(columns=["chromosome", "start", "end", "length"])
    start_position = 0
    current_position = 0
    current_chr = ""
    counter = 0
    df = pd.read_csv(
        "data_files\\full_data.txt", header=None, names=sample_column_names, sep="\t",
    )
    start_position = df.at[0, "position"]
    current_position = start_position
    current_chr = df.at[0, "chromosome"]
    for i in np.arange(len(df.index)):
        chromosome = df.at[i, "chromosome"]
        position = df.at[i, "position"]
        if chromosome != current_chr or (
            position != current_position and position != (current_position + 1)
        ):
            print("{}: {} to {}".format(current_chr, start_position, current_position))
            seq_df.append(
                {
                    "chromosome": current_chr,
                    "start": start_position,
                    "end": current_position,
                    "length": current_position - start_position + 1,
                },
                ignore_index=True,
            )

            current_chr = chromosome
            start_position = position
            current_position = start_position
        elif position == current_position + 1:
            current_position += 1
    print("{};{};{}".format(current_chr, start_position, current_position))
    seq_df.append(
        {
            "chromosome": current_chr,
            "start": start_position,
            "end": current_position,
            "length": current_position - start_position + 1,
        },
        ignore_index=True,
    )
    seq_df.to_csv("data_files\\sequences2.txt")
    return


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
        chunk["sub rate"] = chunk["num subs"] / chunk["num consensus molecules"]
        chunk["downsample"] = rng.binomial(n=N_0, p=chunk["sub rate"])
        chunk.to_csv("data_files\\downsampled_data.txt", mode="a", header=header)
        header = False


def group_by_position(sequence_number):
    df = seq_data_df(sequence_number)
    df = df.drop(columns=["sample ID", "sub rate"])
    aggregation_functions = {
        "chromosome": "first",
        "sub": "first",
        "num subs": "sum",
        "num consensus molecules": sum,
        "downsample": "sum",
    }
    df = df.groupby(df["position"]).aggregate(aggregation_functions)
    df.to_csv(
        "data_files\\sequences_by_position\\seq_{}_group_positions.csv".format(
            sequence_number
        )
    )


def group_by_position_wrapper():
    for i in np.arange(0, 1063):
        group_by_position(i)
        print(i)


# downsample(0.1)


def separating_sequences_wrap():
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


# separating_sequences_wrap()

