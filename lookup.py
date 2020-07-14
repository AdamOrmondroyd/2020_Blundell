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

seq_df = pd.read_csv(
    "data_files\\illumina_80Genes_panel.bed",
    header=None,
    names=["chromosome", "start", "end", "who tf knows", "length", "strand"],
    sep="\t",
)
seq_df.sort_values(
    "chromosome", inplace=True, kind="mergesort", ignore_index=True
)  # use mergesort for stability
seq_df.to_csv("spam.csv")


def seq_data_df(sequence_number):
    """
    Returns DataFrame from seq_(sequence_number).csv
    """
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
        pd.read_csv("data_files\\downsampled_data.txt", chunksize=CHUNKSIZE,)
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
