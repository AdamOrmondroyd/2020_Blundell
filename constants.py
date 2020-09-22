"""
# 2020 Blundell lab intership

Contains constants associated with working with the ALSPAC data.
"""
import numpy as np

CHUNKSIZE = 10 ** 6

AGES = np.array([0, 7, 17, 24])
BASES = np.array(["A", "C", "G", "T"])
LANES = np.array(["lane1", "lane2", "lane3", "lane4"])
age_lane_map = {0: "lane1", 7: "lane2", 17: "lane3", 24: "lane4"}
lane_age_map = {"lane1": 0, "lane2": 7, "lane3": 17, "lane4": 24}
VARIANTS = np.array(
    ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
)
base_variants_map = {
    "A": VARIANTS[:3],
    "C": VARIANTS[3:6],
    "G": VARIANTS[6:9],
    "T": VARIANTS[9:12],
}
base_color_map = {
    "A": "red",
    "C": "lime",
    "G": "magenta",
    "T": "black",
}
variant_color_map = {
    "AC": "lime",
    "AG": "blueviolet",
    "AT": "dodgerblue",
    "CA": "red",
    "CG": "blueviolet",
    "CT": "dodgerblue",
    "GA": "red",
    "GC": "lime",
    "GT": "dodgerblue",
    "TA": "red",
    "TC": "lime",
    "TG": "blueviolet",
}

variant_contexts_map = {
    "AC": np.array(
        [
            "AAA",
            "AAC",
            "AAG",
            "AAT",
            "CAA",
            "CAC",
            "CAG",
            "CAT",
            "GAA",
            "GAC",
            "GAG",
            "GAT",
            "TAA",
            "TAC",
            "TAG",
            "TAT",
        ]
    ),
    "AG": np.array(
        [
            "AAA",
            "AAC",
            "AAG",
            "AAT",
            "CAA",
            "CAC",
            "CAG",
            "CAT",
            "GAA",
            "GAC",
            "GAG",
            "GAT",
            "TAA",
            "TAC",
            "TAG",
            "TAT",
        ]
    ),
    "AT": np.array(
        [
            "AAA",
            "AAC",
            "AAG",
            "AAT",
            "CAA",
            "CAC",
            "CAG",
            "CAT",
            "GAA",
            "GAC",
            "GAG",
            "GAT",
            "TAA",
            "TAC",
            "TAG",
            "TAT",
        ]
    ),
    "CA": np.array(
        [
            "ACA",
            "ACC",
            "ACG",
            "ACT",
            "CCA",
            "CCC",
            "CCG",
            "CCT",
            "GCA",
            "GCC",
            "GCG",
            "GCT",
            "TCA",
            "TCC",
            "TCG",
            "TCT",
        ]
    ),
    "CG": np.array(
        [
            "ACA",
            "ACC",
            "ACG",
            "ACT",
            "CCA",
            "CCC",
            "CCG",
            "CCT",
            "GCA",
            "GCC",
            "GCG",
            "GCT",
            "TCA",
            "TCC",
            "TCG",
            "TCT",
        ]
    ),
    "CT": np.array(
        [
            "ACA",
            "ACC",
            "ACG",
            "ACT",
            "CCA",
            "CCC",
            "CCG",
            "CCT",
            "GCA",
            "GCC",
            "GCG",
            "GCT",
            "TCA",
            "TCC",
            "TCG",
            "TCT",
        ]
    ),
    "GA": np.array(
        [
            "AGA",
            "AGC",
            "AGG",
            "AGT",
            "CGA",
            "CGC",
            "CGG",
            "CGT",
            "GGA",
            "GGC",
            "GGG",
            "GGT",
            "TGA",
            "TGC",
            "TGG",
            "TGT",
        ]
    ),
    "GC": np.array(
        [
            "AGA",
            "AGC",
            "AGG",
            "AGT",
            "CGA",
            "CGC",
            "CGG",
            "CGT",
            "GGA",
            "GGC",
            "GGG",
            "GGT",
            "TGA",
            "TGC",
            "TGG",
            "TGT",
        ]
    ),
    "GT": np.array(
        [
            "AGA",
            "AGC",
            "AGG",
            "AGT",
            "CGA",
            "CGC",
            "CGG",
            "CGT",
            "GGA",
            "GGC",
            "GGG",
            "GGT",
            "TGA",
            "TGC",
            "TGG",
            "TGT",
        ]
    ),
    "TA": np.array(
        [
            "ATA",
            "ATC",
            "ATG",
            "ATT",
            "CTA",
            "CTC",
            "CTG",
            "CTT",
            "GTA",
            "GTC",
            "GTG",
            "GTT",
            "TTA",
            "TTC",
            "TTG",
            "TTT",
        ]
    ),
    "TC": np.array(
        [
            "ATA",
            "ATC",
            "ATG",
            "ATT",
            "CTA",
            "CTC",
            "CTG",
            "CTT",
            "GTA",
            "GTC",
            "GTG",
            "GTT",
            "TTA",
            "TTC",
            "TTG",
            "TTT",
        ]
    ),
    "TG": np.array(
        [
            "ATA",
            "ATC",
            "ATG",
            "ATT",
            "CTA",
            "CTC",
            "CTG",
            "CTT",
            "GTA",
            "GTC",
            "GTG",
            "GTT",
            "TTA",
            "TTC",
            "TTG",
            "TTT",
        ]
    ),
}

a = np.full(30, "als", dtype="U3")
b = np.arange(1, 31).astype(dtype="U2")
PEOPLE = np.char.add(a, b)


def id_age_map(sample_id):
    return lane_age_map[sample_id[:5]]


data_location = "C:\\Users\\Adam\\Programming_files\\2020_Blundell_data_files"
main_data_file = "\\full_data.txt"
# main_data_file = "\\alspac.all.cleaned.sorted.out"
file_names = {
    "data": data_location + main_data_file,
    "downsampled data": data_location + "\\downsampled_data.txt",
    "Wing exons": data_location + "\\Wing_exons.bed",
    "Wing exons sorted": data_location + "\\Wing_exons_sorted.txt",
    "Caroline tiles": data_location + "\\Caroline_tiles.bed",
    "Caroline tiles sorted": data_location + "\\Caroline_tiles_sorted.txt",
    "tile": data_location + "\\tiles\\tile_{}.csv",
    "tile group IDs": data_location + "\\tiles_by_ID\\tile_{}_group_ID.csv",
    "tile group positions": data_location
    + "\\tiles_by_position\\tile_{}_group_positions.csv",
    "tile t&f": data_location + "\\tiles_t&f\\tile_{}_t&f.csv",
    "tile group IDs t&f": data_location + "\\tiles_by_ID_t&f\\tile_{}_group_ID_t&f.csv",
    "tile group positions t&f": data_location
    + "\\tiles_by_position_t&f\\tile_{}_group_positions_t&f.csv",
    "juicy tiles": data_location + "\\juicy_tiles.csv",
    "found means": data_location + "\\positions_to_plot",
}


def sorter(chromosome):
    """Converts X and Y to 23 and 24 for sorting chromosomes."""
    if chromosome == "X":
        return 23
    elif chromosome == "Y":
        return 24
    else:
        return int(chromosome)


sorter = np.vectorize(sorter)


def row_age_map(sample_id):
    """Gets the age from a data row."""
    return lane_age_map[sample_id[:5]]


row_age_map = np.vectorize(row_age_map)
