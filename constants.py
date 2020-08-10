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
SUBS = np.array(
    ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
)
base_variants_map = {
    "A": SUBS[:3],
    "C": SUBS[3:6],
    "G": SUBS[6:9],
    "T": SUBS[9:12],
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
}


def sorter(chromosome):
    """
    Converts X and Y to 23 and 24 for sorting chromosomes.
    """
    if chromosome == "X":
        return 23
    elif chromosome == "Y":
        return 24
    else:
        return int(chromosome)


vec_sorter = np.vectorize(sorter)
