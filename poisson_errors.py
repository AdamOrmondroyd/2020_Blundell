import numpy as np
from lookup import lookup, LANES, PEOPLE, age_lane_map, change_map, seq_df

# df = lookup(PEOPLE, LANES, seq_df.at[0, "chromosome"], seq_df.at[0, "start"])
chromosome = "chr1"
position = 36931702
change = "AG"
df = lookup(chromosome, position)
print(df)
df = df.loc[df["change"] == change]
df.to_csv("data_files\\spam.csv")
print(df)
num_changes = np.sum(df["num changes"])
total_consensus = np.sum(df["num consensus molecules"])
print(
    "{} out of {} consensus mulecules had a {} change in chromosome {}, position {}".format(
        num_changes, total_consensus, change, chromosome, position
    )
)
print("error rate: {}".format(num_changes / total_consensus))
