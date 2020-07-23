"""
# 2020 Blundell lab internship

Contains functions for accessing the ALSPAC data using Pandas dataframes.
"""
import numpy as np
import pandas as pd
import os
from constants import CHANGES, CHUNKSIZE, LANES, PEOPLE, sub_complement_map


sample_column_names = [
    "chromosome",
    "position",
    "sub",
    "num subs",
    "num consensus molecules",
    "sample ID",
]


# Dataframe of genes
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


def seq_data_df(sequence_number, group_by=None, trimmed_and_flipped=True):
    """
    Returns DataFrame from seq_(sequence_number).csv.

    group_by = "position" combines samples at the same position.
    group_by = "ID" combines samples with the same ID (person and age)
    Both grouping options use the trimmed and flipped data by default.
    """
    if trimmed_and_flipped:
        if group_by == "position":
            return pd.read_csv(
                "data_files\\sequences_by_position_t&f\\seq_{}_group_positions_t&f.csv".format(
                    sequence_number
                )
            )
        elif group_by == "ID":
            return pd.read_csv(
                "data_files\\sequences_by_ID_t&f\\seq_{}_group_ID_t&f.csv".format(
                    sequence_number
                )
            )
        else:
            return pd.read_csv(
                "data_files\\sequences_t&f\\seq_{}_t&f.csv".format(sequence_number),
                index_col=0,
            )
    else:
        if group_by == "position":
            return pd.read_csv(
                "data_files\\sequences_by_position\\seq_{}_group_positions.csv".format(
                    sequence_number
                )
            )
        elif group_by == "ID":
            return pd.read_csv(
                "data_files\\sequences_by_ID\\seq_{}_group_ID.csv".format(
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
        Ns = chunk["num consensus molecules"].to_numpy()
        Ns = Ns.astype(int)
        Nss.append(Ns)
    Nss = np.hstack(Nss)
    return np.percentile(Nss, q)


def downsample(q=50):
    """
    Downsamples full_data.txt to qth percentile of number of consensus molecules.

    Defaults to 50th percentile
    """
    if os.path.isfile("data_files\\downsampled_data.txt"):
        os.remove("data_files\\downsampled_data.txt")
    rng = np.random.default_rng()
    header = True
    N_0 = percentile(q)
    print("Downsampled to {} consensus molecules".format(N_0))
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


def group_by_position(sequence_number, trimmed_and_flipped=True):
    """
    Groups the specified sequence by position.
    """
    df = seq_data_df(sequence_number, trimmed_and_flipped=trimmed_and_flipped)
    df = df.drop(columns=["sample ID"])
    df = df.groupby(["position", "chromosome", "sub"]).agg(aggregation_functions)
    if trimmed_and_flipped:
        df.to_csv(
            "data_files\\sequences_by_position_t&f\\seq_{}_group_positions_t&f.csv".format(
                sequence_number
            )
        )
    else:
        df.to_csv(
            "data_files\\sequences_by_position\\seq_{}_group_positions.csv".format(
                sequence_number
            )
        )


def group_by_ID(sequence_number, trimmed_and_flipped=True):
    """
    Groups the specified sequence by ID (person and age).
    """
    df = seq_data_df(sequence_number, trimmed_and_flipped=trimmed_and_flipped)
    df = df.drop(columns=["position"])
    df = df.groupby(["sample ID", "chromosome", "sub"]).agg(aggregation_functions)
    if trimmed_and_flipped:
        df.to_csv(
            "data_files\\sequences_by_ID_t&f\\seq_{}_group_ID_t&f.csv".format(
                sequence_number
            )
        )
    else:
        df.to_csv(
            "data_files\\sequences_by_ID\\seq_{}_group_ID.csv".format(sequence_number)
        )


def trim_and_flip(gene_number):
    """
    Flips neg data to be what the actual change was.

    The called changes are all relative to the top strand. This function creates files which identify the actual sub seen.
    """
    seqs = gene_seqs_map[gene_number]
    next_df = seq_data_df(seqs.index[0], trimmed_and_flipped=False)
    if len(seqs.index) >= 2:
        for i in seqs.index[:-1]:
            df = next_df
            next_df = seq_data_df(i + 1, trimmed_and_flipped=False)
            cond = ~df["position"].isin(next_df["position"])
            next_cond = ~next_df["position"].isin(df["position"])
            df = df.loc[cond, :]
            next_df = next_df.loc[next_cond, :]
            if seq_df.at[i, "strand"] == "-":
                df = df.replace({"sub": sub_complement_map})
            df.to_csv("data_files\\sequences_t&f\\seq_{}_t&f.csv".format(i))
    next_df.to_csv("data_files\\sequences_t&f\\seq_{}_t&f.csv".format(seqs.index[-1]))


### Wrappers ###


def group_by_position_wrapper(trimmed_and_flipped=True):
    """
    Repeats group by position for all sequences.
    """
    for i in np.arange(0, 1063):
        group_by_position(i, trimmed_and_flipped)
        print(i)


def group_by_ID_wrapper(trimmed_and_flipped=True):
    """
    Repeats group by ID for all sequences.
    """
    for i in np.arange(0, 1063):
        group_by_ID(i, trimmed_and_flipped)
        print(i)


def separating_sequences_wrapper():
    """
    Separates sequences in chunks of 100 to avoid memory issues.
    """
    print("Separating sequences 0 to 99")
    separating_sequences(np.arange(0, 100))
    print("Separating sequences 100 to 199")
    separating_sequences(np.arange(100, 200))
    print("Separating sequences 200 to 299")
    separating_sequences(np.arange(200, 300))
    print("Separating sequences 300 to 399")
    separating_sequences(np.arange(300, 400))
    print("Separating sequences 400 to 499")
    separating_sequences(np.arange(400, 500))
    print("Separating sequences 500 to 599")
    separating_sequences(np.arange(500, 600))
    print("Separating sequences 600 to 699")
    separating_sequences(np.arange(600, 700))
    print("Separating sequences 700 to 799")
    separating_sequences(np.arange(700, 800))
    print("Separating sequences 800 to 899")
    separating_sequences(np.arange(800, 900))
    print("Separating sequences 900 to 999")
    separating_sequences(np.arange(900, 1000))
    print("Separating sequences 1000 to 1062")
    separating_sequences(np.arange(1000, 1063))


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
        separating_sequences_wrapper()
    group_by_ID_wrapper(trimmed_and_flipped=False)
    group_by_position_wrapper(trimmed_and_flipped=False)
    trim_and_flip_wrapper()
    group_by_ID_wrapper()
    group_by_position_wrapper()


### Miscellaneous ###


def empty_sequences():
    """
    Identifies any empty sequences.
    """
    for i, gene in gene_df.iterrows():
        for j, seq in gene_seqs_map[i].iterrows():
            if 0 == len(seq_data_df(j).index):
                print("Empty sequence found, gene {} sequence {}".format(i, j))


def sequences_not_in_genes():
    """
    Identifies any sequences that don't appear in any genes
    """
    seqs_in_genes = []
    for i in gene_df.index:
        seqs_in_genes.extend(gene_seqs_map[i].index)

    for i in seq_df.index:
        if ~np.isin(i, seqs_in_genes):
            print("Sequence {} is not in any gene".format(i))

    print(seqs_in_genes)
