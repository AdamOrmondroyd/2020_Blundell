import numpy as np
import matplotlib.pyplot as plt
from constants import SUBS
from lookup import seq_df

plot_title = "histogram"

seq = 1
positions = np.arange(seq_df.at[seq, "start"], seq_df.at[seq, "end"])
chromosome = seq_df.at[seq, "chromosome"]

df = lookup(chromosome, positions, downsample=True)

for sub in SUBS:
    print(sub)
    fig, ax = plt.subplots()

    change_df = df.loc[sub == df["sub"]]

    mean = np.mean(change_df["num subs"])
    variance = np.var(change_df["num subs"])
    print(
        "mean:{}, variance: {}, ratio (variance/mean): {}".format(
            mean, variance, variance / mean
        )
    )
    # print(np.max(change_df["num subs"]))
    # bins = np.arange(-0.5, np.max(change_df["num subs"]) + 0.5)
    # print(bins)
    # ax.hist(change_df["num subs"], bins=bins)
    # ax.set(title="{}_{}".format(plot_title, sub))
    # plt.show()
