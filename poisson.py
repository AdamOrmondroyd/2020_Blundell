import numpy as np
import matplotlib.pyplot as plt
from constants import VARIANTS
from lookup import tile_df

plot_title = "histogram"

tile = 1
positions = np.arange(tile_df.at[tile, "start"], tile_df.at[tile, "end"])
chromosome = tile_df.at[tile, "chromosome"]

df = lookup(chromosome, positions, downsample=True)

for sub in VARIANTS:
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
