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
base_changes_map = {
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
change_color_map = {
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
