import numpy as np
import matplotlib.pyplot as plt
from constants import (
    CHANGES,
    LANES,
    PEOPLE,
    age_lane_map,
    base_change_map,
    change_color_map,
)
from lookup import lookup, seq_df

# df = lookup(PEOPLE, LANES, seq_df.at[0, "chromosome"], seq_df.at[0, "start"])
chromosome = "chr1"
positions = np.arange(seq_df.at[0, "start"], seq_df.at[0, "end"])

df = lookup(chromosome, positions)

fig, ax = plt.subplots()
plot_title = "{} {} to {}".format(chromosome, positions[0], positions[-1])

for change in CHANGES:
    # for change in base_change_map["A"]:
    print(change)
    df_change = df.loc[df["change"] == change]
    # df.to_csv("data_files\\spam.csv")
    error_rates = np.zeros(positions.size)
    for i, position in enumerate(positions):
        df_position = df_change.loc[df_change["position"] == position]
        if df_position.shape[0] != 0:
            num_changes = np.sum(df_position["num changes"])
            total_consensus = np.sum(df_position["num consensus molecules"])
            error_rates[i] = num_changes / total_consensus
            # if 0 == error_rates[i]:
            #     print("zero at{}".format(position))
    ax.plot(
        positions,
        error_rates,
        label=change,
        linestyle="None",
        marker="+",
        color=change_color_map[change],
    )

# don't use standard form for axes
ax.ticklabel_format(useOffset=False, style="plain")

ax.set(title=plot_title)
ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
plt.savefig("plots\\" + plot_title + ".png")
ax.set(yscale="log")
fig.tight_layout()
plt.savefig("plots\\" + plot_title + "log.png")
plt.show()
