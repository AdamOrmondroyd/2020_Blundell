import numpy as np
import pandas as pd
from constants import CHANGES, LANES, PEOPLE

column_names = [
    "chromosome",
    "position",
    "change",
    "num changes",
    "num consensus molecules",
    "sample ID",
]


id_df = pd.read_csv(
    "data_files\\id.txt",
    header=None,
    names=["lane1", "lane2", "lane3", "lane4", "ID"],
    sep=";",
    index_col=-1,
)

seq_df = pd.read_csv(
    "data_files\\sequences.txt",
    header=None,
    names=["chromosome", "start", "end"],
    sep=";",
)


def lookup(chromosome, positions, change=CHANGES, people=PEOPLE, lanes=LANES):
    """
    Looks up a given people, ages (given by lanes) and transitions (e.g. AC).
    "lane1" = age0, "lane2" = age7, "lane3" = age17, "lane4" = age24
    """
    sample_ids = np.zeros(
        people.size * lanes.size, dtype="U14"
    )  # string length 14, keep an eye on this

    for j in range(people.size):
        for i in range(lanes.size):
            sample_ids[j * lanes.size + i] = id_df.at[people[j], lanes[i]][:14]

    # Empty dataframe
    df = pd.DataFrame(columns=column_names)

    df = pd.read_csv(
        "data_files\\full_data.txt", header=None, names=column_names, sep="\t",
    )
    df = df.loc[
        np.isin(df["sample ID"], sample_ids)
        & np.isin(df["position"], positions)
        & (df["chromosome"] == chromosome)
        & np.isin(df["change"], change)
    ]
    return df


def figuring_out_the_data():
    start_position = 0
    current_position = 0
    current_chr = ""
    counter = 0
    df = pd.read_csv(
        "data_files\\full_data.txt", header=None, names=column_names, sep="\t",
    )
    for i in np.arange(len(df.index)):
        chromosome = df.at[i, "chromosome"]
        position = df.at[i, "position"]
        if chromosome != current_chr or (
            position != current_position and position != (current_position + 1)
        ):
            print("{}: {} to {}".format(current_chr, start_position, current_position))

            current_chr = chromosome
            start_position = position
            current_position = start_position
        elif position == current_position + 1:
            current_position += 1
    print("{};{};{}".format(current_chr, start_position, current_position))
    return
