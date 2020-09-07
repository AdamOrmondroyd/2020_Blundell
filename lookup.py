"""
# 2020 Blundell lab internship

Contains functions for accessing the ALSPAC data using Pandas dataframes.
"""
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import os
import gc
from constants import (
    CHUNKSIZE,
    file_names,
    LANES,
    PEOPLE,
    sorter,
)


sample_column_names = [
    "chromosome",
    "position",
    "variant",
    "num variants",
    "num consensus molecules",
    "sample ID",
]

### Preparing files ###


def percentile(q):
    """Finds the qth percentile of the number of consensus molecules."""
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

    Defaults to 50th percentile.
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
        chunk = chunk.loc[chunk["num consensus molecules"] >= N_0]
        chunk["downsample"] = rng.binomial(
            n=N_0, p=chunk["num variants"] / chunk["num consensus molecules"]
        )
        chunk.to_csv(file_names["downsampled data"], mode="a", header=header)
        header = False
        gc.collect()


def sort_caroline_tiles():
    """Orders data from Caroline by chromosome."""
    df = pd.read_csv(
        file_names["Caroline tiles"],
        header=None,
        names=["chromosome", "start", "end", "who tf knows", "length", "strand"],
        sep="\t",
    )
    df["chromosome"] = df["chromosome"].str[3:]
    df = df.sort_values(
        by="chromosome", kind="mergesort", key=sorter, ignore_index=True,
    )  # use mergesort for stability

    def assign_exon_number_to_tile(df_row):
        """Assigns the exon number to the tiles in tile_df."""
        for j, exon in exon_df.iterrows():
            if (
                df_row["chromosome"] == exon["chromosome"]
                and df_row["start"] >= exon["start"]
                and df_row["end"] <= exon["end"]
            ):
                return j
        return -1  # Error value if no matching exon is found

    # assign_exon_number_to_tile = np.vectorize(assign_exon_number_to_tile)

    df["exon"] = df.apply(assign_exon_number_to_tile, axis=1)
    df.to_csv(file_names["Caroline tiles sorted"], sep="\t")


def sort_wing_exons():
    """Orders exon data from Wing."""
    df = pd.read_csv(
        file_names["Wing exons"],
        header=None,
        names=["chromosome", "start", "end"],
        sep="\t",
    )
    df["chromosome"] = df["chromosome"].str[3:]
    df = df.sort_values(
        by="chromosome", kind="mergesort", key=sorter, ignore_index=True,
    )  # use mergesort for stability
    df.to_csv(file_names["Wing exons sorted"], sep="\t")


def separating_tiles(tile_numbers):
    """
    Separates the tiles in full_data.txt by the tiles in Caroline's file.

    tile_numbers = list of tiles to produce files for.

    Note: cannot do all tiles at once.
    """
    reduced_tile_df = tile_df.loc[tile_numbers]
    tile_dfs = []
    for i in range(len(tile_numbers)):
        tile_dfs.append(pd.DataFrame(columns=sample_column_names))

    for j, chunk in enumerate(
        pd.read_csv(file_names["downsampled data"], chunksize=CHUNKSIZE, index_col=0)
    ):
        print("chunk {}".format(j))
        for i, (k, tile) in zip(range(tile_numbers.size), reduced_tile_df.iterrows()):
            # print("tile {}".format(i))

            tile_dfs[i] = tile_dfs[i].append(
                chunk.loc[
                    (chunk["position"] >= tile["start"])
                    & (chunk["position"] <= tile["end"])
                ],
                ignore_index=True,
            )
        gc.collect()

    for i, df in zip(tile_numbers, tile_dfs):
        df.to_csv(file_names["tile"].format(i))


aggregation_functions = {
    "num variants": "sum",
    "num consensus molecules": "sum",
    "downsample": "sum",
}


def tile_data_df(tile_number, group_by=None, trim_and_flip=True):
    """
    Returns DataFrame from tile_(tile_number).csv.

    group_by = "position" combines samples at the same position.
    group_by = "ID" combines samples with the same ID (person and age)
    Both grouping options use the trimmed and flipped data by default.
    """
    if trim_and_flip:
        if group_by == "ID":
            return pd.read_csv(file_names["tile group IDs t&f"].format(tile_number))
        elif group_by == "position":
            return pd.read_csv(
                file_names["tile group positions t&f"].format(tile_number)
            )
        else:
            return pd.read_csv(file_names["tile t&f"].format(tile_number), index_col=0,)
    else:
        if group_by == "ID":
            return pd.read_csv(file_names["tile group IDs"].format(tile_number))
        elif group_by == "position":
            return pd.read_csv(file_names["tile group positions"].format(tile_number))
        else:
            return pd.read_csv(file_names["tile"].format(tile_number), index_col=0)


def group_by_position(tile_number, trim_and_flip=True):
    """Groups the specified tile by position."""
    df = tile_data_df(tile_number, trim_and_flip=trim_and_flip)
    df = df.drop(columns=["sample ID"])
    df = df.groupby(["position", "chromosome", "variant"]).agg(aggregation_functions)
    if trim_and_flip:
        df.to_csv(file_names["tile group positions t&f"].format(tile_number))
    else:
        df.to_csv(file_names["tile group positions"].format(tile_number))


def group_by_ID(tile_number, trim_and_flip=True):
    """Groups the specified tile by ID (person and age)."""
    df = tile_data_df(tile_number, trim_and_flip=trim_and_flip)
    df = df.drop(columns=["position"])
    df = df.groupby(["sample ID", "chromosome", "variant"]).agg(aggregation_functions)
    if trim_and_flip:
        df.to_csv(file_names["tile group IDs t&f"].format(tile_number))
    else:
        df.to_csv(file_names["tile group IDs"].format(tile_number))


def trim_and_flip(exon_number):
    """
    Flips neg data to be what the actual change was.

    The called changes are all relative to the top strand. This function creates files which identify the actual variant seen.
    """
    downsample_limit = 1500

    def flip(sequence_string):
        """Returns complement of the variant as a Biopython Seq"""
        return str(Seq(sequence_string).complement())

    flip = np.vectorize(flip)

    def read_genome(tile_number, tile_data_df):
        """A version of read_genome which just returns the string when given the tile directly."""
        genome_string = ""
        df = tile_data_df.drop(columns=["sample ID"])
        df = df.groupby(["position", "chromosome", "variant"]).agg(
            aggregation_functions
        )
        df = df.reset_index()
        for i, row in df.iterrows():
            if i % 3 == 0:
                genome_string += row["variant"][0]
        tile_df.at[tile_number, "genome"] = genome_string
        tile_df.to_csv(file_names["Caroline tiles sorted"], sep="\t")
        return genome_string

    def context(tile_number, tile_data_df):
        """Will think of an accurate description later."""

        start = np.min(tile_data_df["position"])
        end = np.max(tile_data_df["position"])

        genome_string = read_genome(tile_number, tile_data_df)

        for position in pd.unique(tile_data_df["position"]):
            if position <= start or position >= end:
                context_string = "?"
            else:
                index = position - start
                context_string = genome_string[index - 1 : index + 2]
            tile_data_df.loc[
                tile_data_df.position == position, "context"
            ] = context_string
        return tile_data_df

    tiles = exon_tiles_map[exon_number]
    next_df = tile_data_df(tiles.index[0], trim_and_flip=False)
    if len(tiles.index) >= 2:
        for i in tiles.index[:-1]:
            df = next_df
            next_df = tile_data_df(i + 1, trim_and_flip=False)

            cond = ~df["position"].isin(next_df["position"])
            next_cond = ~next_df["position"].isin(df["position"])
            df = df.loc[cond, :]
            next_df = next_df.loc[next_cond, :]

            if tile_df.at[i, "strand"] == "-" and len(df.index) != 0:
                df["variant"] = flip(df["variant"])

            df = df.loc[df["downsample"] <= downsample_limit]

            if not df.empty:
                df = context(i, df)

            df.to_csv(file_names["tile t&f"].format(i))

    if len(tiles.index) > 0:
        i = tiles.index[-1]

        if not next_df.empty:
            if tile_df.at[i, "strand"] == "-" and len(next_df.index) != 0:
                next_df["variant"] = flip(next_df["variant"])

            next_df = next_df.loc[next_df["downsample"] <= downsample_limit]

            next_df = context(i, next_df)

        next_df.to_csv(file_names["tile t&f"].format(tiles.index[-1]))


### Wrappers ###


def group_by_position_wrapper(trim_and_flip=True):
    """Repeats group by position for all tiles."""
    for i in np.arange(0, 1063):
        group_by_position(i, trim_and_flip)
        print(i)


def group_by_ID_wrapper(trim_and_flip=True):
    """Repeats group by ID for all tiles."""
    for i in np.arange(0, 1063):
        group_by_ID(i, trim_and_flip)
        print(i)


def separating_tiles_wrapper():
    """Separates tiles in chunks of 100 to avoid memory issues."""
    print("Separating tiles 0 to 99")
    separating_tiles(np.arange(0, 100))
    print("Separating tiles 100 to 199")
    separating_tiles(np.arange(100, 200))
    print("Separating tiles 200 to 299")
    separating_tiles(np.arange(200, 300))
    print("Separating tiles 300 to 399")
    separating_tiles(np.arange(300, 400))
    print("Separating tiles 400 to 499")
    separating_tiles(np.arange(400, 500))
    print("Separating tiles 500 to 599")
    separating_tiles(np.arange(500, 600))
    print("Separating tiles 600 to 699")
    separating_tiles(np.arange(600, 700))
    print("Separating tiles 700 to 799")
    separating_tiles(np.arange(700, 800))
    print("Separating tiles 800 to 899")
    separating_tiles(np.arange(800, 900))
    print("Separating tiles 900 to 999")
    separating_tiles(np.arange(900, 1000))
    print("Separating tiles 1000 to 1062")
    separating_tiles(np.arange(1000, 1063))


def group_strands_wrapper():
    """Repeats group strands for all exons"""
    for i in np.arange(len(exon_df.index)):
        print(i)
        group_strands(i)


def trim_and_flip_wrapper():
    for i in exon_df.index:
        trim_and_flip(i)
        print("Exon {}".format(i))


def refresh_data(redownsample=False, just_trim_and_flip=False):
    """Runs all the functions in turn to freshen up the data"""
    if not just_trim_and_flip:
        if redownsample:
            downsample()
            separating_tiles_wrapper()
        print("Grouping by ID")
        group_by_ID_wrapper(trim_and_flip=False)
        print("Grouping by position")
        group_by_position_wrapper(trim_and_flip=False)
    print("Trimming and flipping")
    trim_and_flip_wrapper()
    print("Grouping by ID (t&f)")
    group_by_ID_wrapper()
    print("Grouping by position (t&f)")
    group_by_position_wrapper()


### Dataframes ###


# Dataframe of exons
exon_df = pd.read_csv(file_names["Wing exons sorted"], sep="\t", index_col=0,)

# DataFrame of tile information
tile_df = pd.read_csv(file_names["Caroline tiles sorted"], sep="\t", index_col=0)

exon_tiles_map = {}
for i, exon in exon_df.iterrows():
    exon_tiles_map[i] = tile_df.loc[
        (tile_df["start"] >= exon["start"]) & (tile_df["end"] <= exon["end"])
    ]


chromosome_tiles_map = {}
chromosome_exon_map = {}
chromosomes = tile_df["chromosome"].unique()
for chromosome in chromosomes:
    chromosome_tiles_map[chromosome] = tile_df.loc[tile_df["chromosome"] == chromosome]
    chromosome_exon_map[chromosome] = exon_df.loc[exon_df["chromosome"] == chromosome]

juicy_df = pd.read_csv(file_names["juicy tiles"], index_col=0)


### Miscellaneous ###


def empty_tiles():
    """Identifies any empty tiles."""
    for i, exon in exon_df.iterrows():
        for j, tile in exon_tiles_map[i].iterrows():
            if 0 == len(tile_data_df(j).index):
                print("Empty tile found, exon {} tile {}".format(i, j))


def tiles_not_in_exons():
    """Identifies any tiles that don't appear in any exons"""
    tiles_in_exons = []
    for i in exon_df.index:
        tiles_in_exons.extend(exon_tiles_map[i].index)

    for i in tile_df.index:
        if ~np.isin(i, tiles_in_exons):
            print("Sequence {} is not in any exon".format(i))

    print(tiles_in_exons)


def read_genome():
    """
    Uses the first letter of the first variant at each position to work out what the "correct" genome is.

    Uses the trimmed and flipped data, so the genome is on the negative strand for those.
    """
    tile_genomes = np.empty(len(tile_df.index), dtype=object)
    for j in tile_df.index:
        tile_genomes[j] = ""
        df = tile_data_df(j, group_by="position", trim_and_flip=True)
        for i, row in df.iterrows():
            if i % 3 == 0:
                tile_genomes[j] += row["variant"][0]
        print(tile_genomes[j])
    tile_df["genome"] = tile_genomes
    tile_df.to_csv(file_names["Caroline tiles sorted"], sep="\t")


def variants_per_position(threshold=1000):
    """Spits out locations with juicy number of variants above threshold."""
    juicy_df = pd.DataFrame(
        columns=["tile", "chromosome", "position", "variant", "downsample"]
    )
    index = 0
    for j, tile in tile_df.iterrows():
        print(j)
        df = tile_data_df(j, group_by="position")
        for i, row in df.iterrows():
            if row["downsample"] >= threshold:
                juicy_df.loc[index] = [
                    j,
                    tile["chromosome"],
                    row["position"],
                    row["variant"],
                    row["downsample"],
                ]
                index += 1
    juicy_df.to_csv(file_names["juicy tiles"])


def actual_differences():
    """Looks in the juicy locations to find actual genetic differences between people."""
    for j, juicy_row in juicy_df.iterrows():
        df = tile_data_df(juicy_row["tile"])
        df = df.loc[
            (df["position"] == juicy_row["position"])
            & (df["variant"] == juicy_row["variant"])
            & (df["downsample"] > 0)
        ]
        print(df)

