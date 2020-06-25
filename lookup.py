import numpy as np
import pandas as pd

AGES = np.array([0, 7, 17, 24])
age_lane_map = {0: "lane1", 7: "lane2", 17: "lane3", 24: "lane4"}

CHANGES = np.array(
    [
        "TA",
        "TC",
        "TG",
        "AC",
        "AG",
        "AT",
        "CA",
        "CG",
        "CT",
        "AC",
        "AG",
        "AT",
        "AC",
        "AG",
        "AT",
    ]
)

id_df = pd.read_csv(
    "data_files\\id.txt",
    header=None,
    names=["lane1", "lane2", "lane3", "lane4", "ID"],
    sep=";",
    index_col=-1,
)


def lookup_frequency(id, lane, chromosome, position, change):
    """
    Looks up a given person, age (given by lane) and transition (e.g. AC)
    "lane1" = age0, "lane2" = age7, "lane3" = age17, "lane4" = age24
    """
    sample_id = id_df.at[id, lane][:14]

    chunksize = 10 ** 6  # number of rows per chunk

    for chunk in pd.read_csv(
        "data_files\\full_data.txt",
        chunksize=chunksize,
        header=None,
        names=[
            "chromosome",
            "position",
            "change",
            "frequency",
            "num consensus molecules",
            "sample ID",
        ],
        sep="\t",
    ):

        a = chunk.loc[
            (chunk["sample ID"] == sample_id)
            & (chunk["change"] == change)
            & (chunk["position"] == position)
            & (chunk["chromosome"] == chromosome)
        ]["frequency"].values.tolist()

        if a:
            return a[0]
    return
