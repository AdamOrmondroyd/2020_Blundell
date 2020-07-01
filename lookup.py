import numpy as np
import pandas as pd

CHUNKSIZE = 10 ** 6  # number of rows per chunk

AGES = np.array([0, 7, 17, 24])
LANES = np.array(["lane1", "lane2", "lane3", "lane4"])
age_lane_map = {0: "lane1", 7: "lane2", 17: "lane3", 24: "lane4"}
lane_age_map = {"lane1": 0, "lane2": 7, "lane3": 17, "lane4": 24}
change_map = {
    "A": np.array(["AC", "AG", "AT"]),
    "C": np.array(["CA", "CG", "CT"]),
    "G": np.array(["GA", "GC", "GT"]),
    "T": np.array(["TA", "TC", "TG"]),
}

CHANGES = np.array(
    ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
)

a = np.full(30, "als", dtype="U3")
b = np.arange(1, 31).astype(dtype="U2")
PEOPLE = np.char.add(a, b)

column_names = [
    "chromosome",
    "position",
    "change",
    "num changes",
    "num consensus molecules",
    "sample ID",
]


def id_age_map(sample_id):
    return lane_age_map[sample_id[:5]]


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


def lookup(chromosome, position, people=PEOPLE, lanes=LANES):
    """
    Looks up a given people, ages (given by lanes) and transitions (e.g. AC).
    "lane1" = age0, "lane2" = age7, "lane3" = age17, "lane4" = age24
    """
    sample_ids = np.zeros(
        (people.size, lanes.size), dtype="U14"
    )  # string length 14, keep an eye on this

    for j in range(people.size):
        for i in range(lanes.size):
            sample_ids[j, i] = id_df.at[people[j], lanes[i]][:14]

    # Empty dataframe
    df = pd.DataFrame(columns=column_names)

    for chunk in pd.read_csv(
        "data_files\\full_data.txt",
        chunksize=CHUNKSIZE,
        header=None,
        names=column_names,
        sep="\t",
    ):

        chunk = chunk.loc[
            (np.any(np.in1d(sample_ids, chunk["sample ID"])))
            & (chunk["position"] == position)
            & (chunk["chromosome"] == chromosome)
        ]
        print(chunk)
        df = df.append(chunk, ignore_index=True)

    return df


def figuring_out_the_data():
    start_position = 0
    current_position = 0
    current_chr = ""
    counter = 0
    for chunk in pd.read_csv(
        "data_files\\full_data.txt",
        chunksize=CHUNKSIZE,
        header=None,
        names=column_names,
        sep="\t",
    ):
        for i in np.arange(len(chunk)) + counter * CHUNKSIZE:
            chromosome = chunk.at[i, "chromosome"]
            position = chunk.at[i, "position"]
            if chromosome != current_chr or (
                position != current_position and position != (current_position + 1)
            ):
                print(
                    "{}: {} to {}".format(current_chr, start_position, current_position)
                )

                current_chr = chromosome
                start_position = position
                current_position = start_position
            elif position == current_position + 1:
                current_position += 1
        print("new chonker: {}, position = {}".format(current_chr, current_position))
        counter += 1
    print("{}: {} to {}".format(current_chr, start_position, current_position))
    return
