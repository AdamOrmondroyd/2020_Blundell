import numpy as np
import matplotlib.pyplot as plt
from constants import CHANGES
from lookup import lookup, seq_df

plot_title = "histogram"

sequence = 1
positions = np.arange(seq_df.at[sequence, "start"], seq_df.at[sequence, "end"])
chromosome = seq_df.at[sequence, "chromosome"]

df = lookup(chromosome, positions, downsample=True)

for change in CHANGES:
    print(change)
    fig, ax = plt.subplots()

    change_df = df.loc[change == df["change"]]

    mean = np.mean(change_df["num changes"])
    variance = np.var(change_df["num changes"])
    print(
        "mean:{}, variance: {}, ratio (variance/mean): {}".format(
            mean, variance, variance / mean
        )
    )
    # print(np.max(change_df["num changes"]))
    # bins = np.arange(-0.5, np.max(change_df["num changes"]) + 0.5)
    # print(bins)
    # ax.hist(change_df["num changes"], bins=bins)
    # ax.set(title="{}_{}".format(plot_title, change))
    # plt.show()
