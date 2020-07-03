import numpy as np

AGES = np.array([0, 7, 17, 24])
LANES = np.array(["lane1", "lane2", "lane3", "lane4"])
age_lane_map = {0: "lane1", 7: "lane2", 17: "lane3", 24: "lane4"}
lane_age_map = {"lane1": 0, "lane2": 7, "lane3": 17, "lane4": 24}
CHANGES = np.array(
    ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]
)
base_change_map = {
    "A": CHANGES[:3],
    "C": CHANGES[3:6],
    "G": CHANGES[6:9],
    "T": CHANGES[9:12],
}
change_color_map = {
    "AC": "lightcoral",
    "AG": "red",
    "AT": "darkred",
    "CA": "lightgreen",
    "CG": "lime",
    "CT": "darkgreen",
    "GA": "deeppink",
    "GC": "magenta",
    "GT": "darkmagenta",
    "TA": "lightgrey",
    "TC": "grey",
    "TG": "black",
}

a = np.full(30, "als", dtype="U3")
b = np.arange(1, 31).astype(dtype="U2")
PEOPLE = np.char.add(a, b)


def id_age_map(sample_id):
    return lane_age_map[sample_id[:5]]
