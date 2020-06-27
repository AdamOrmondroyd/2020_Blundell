import numpy as np
import pandas as pd

AGES = np.array([0, 7, 17, 24])
LANES = np.array(["lane1", "lane2", "lane3", "lane4"])
age_lane_map = {0: "lane1", 7: "lane2", 17: "lane3", 24: "lane4"}
lane_age_map = {"lane1": 0, "lane2": 7, "lane3": 17, "lane4": 24}

CHANGES = np.array(
    ["TA", "TC", "TG", "AC", "AG", "AT", "CA", "CG", "CT", "AC", "AG", "AT",]
)

id_df = pd.read_csv(
    "data_files\\id.txt",
    header=None,
    names=["lane1", "lane2", "lane3", "lane4", "ID"],
    sep=";",
    index_col=-1,
)


def lookup(person, lanes, chromosome, position, changes):
    """
    Looks up a given person, ages (given by lanes) and transitions (e.g. AC).
    "lane1" = age0, "lane2" = age7, "lane3" = age17, "lane4" = age24
    """
    sample_ids = np.zeros(
        lanes.size, dtype="U14"
    )  # string length 14, keep an eye on this

    for i in range(lanes.size):
        print(id_df.at[person, lanes[i]][:14])
        sample_ids[i] = id_df.at[person, lanes[i]][:14]

    print(sample_ids)
    chunksize = 10 ** 6  # number of rows per chunk
    column_names = [
        "chromosome",
        "position",
        "change",
        "frequency",
        "num consensus molecules",
        "sample ID",
    ]
    # Empty dataframe
    df = pd.DataFrame(columns=column_names)

    for chunk in pd.read_csv(
        "data_files\\full_data.txt",
        chunksize=chunksize,
        header=None,
        names=column_names,
        sep="\t",
    ):

        df_chunk = chunk.loc[
            (chunk["sample ID"].isin(sample_ids))
            & (chunk["change"].isin(changes))
            & (chunk["position"] == position)
            & (chunk["chromosome"] == chromosome)
        ]
        print(df_chunk)
        df = df.append(df_chunk, ignore_index=True)

    return df
