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
CHANGES = np.array(
    ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
)
base_subs_map = {
    "A": CHANGES[:3],
    "C": CHANGES[3:6],
    "G": CHANGES[6:9],
    "T": CHANGES[9:12],
}
base_color_map = {
    "A": "red",
    "C": "lime",
    "G": "magenta",
    "T": "black",
}
sub_color_map = {
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


sub_complement_map = {
    "AC": "TG",
    "AG": "TC",
    "AT": "TA",
    "CA": "GT",
    "CG": "GC",
    "CT": "GA",
    "GA": "CT",
    "GC": "CG",
    "GT": "CA",
    "TA": "AT",
    "TC": "AG",
    "TG": "AC",
}

data_location = "C:\\Users\\Adam\\Programming_files\\2020_Blundell_data_files"
file_names = {
    "data": data_location + "\\full_data.txt",
    "downsampled data": data_location + "\\downsampled_data.txt",
    "Wing genes": data_location + "\\Wing_genes.bed",
    "Caroline seqs": data_location + "\\Caroline_sequences.bed",
    "seq": data_location + "\\sequences\\seq_{}.csv",
    "seq group IDs": data_location + "\\sequences_by_ID\\seq_{}_group_ID.csv",
    "seq group positions": data_location
    + "\\sequences_by_position\\seq_{}_group_positions.csv",
    "seq t&f": data_location + "\\sequences_t&f\\seq_{}_t&f.csv",
    "seq group IDs t&f": data_location
    + "\\sequences_by_ID_t&f\\seq_{}_group_ID_t&f.csv",
    "seq group positions t&f": data_location
    + "\\sequences_by_position_t&f\\seq_{}_group_positions_t&f.csv",
}
