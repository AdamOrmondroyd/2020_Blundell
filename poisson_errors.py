import numpy as np
import matplotlib.pyplot as plt
from constants import (
    CHANGES,
    LANES,
    PEOPLE,
    age_lane_map,
    base_change_map,
    base_color_map,
    change_color_map,
)
from lookup import lookup, seq_df

# df = lookup(PEOPLE, LANES, seq_df.at[0, "chromosome"], seq_df.at[0, "start"])
chromosome = "chr1"
positions = np.arange(seq_df.at[0, "start"], seq_df.at[0, "end"])

df = lookup(chromosome, positions)

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
plot_title = "{}_{}-{}".format(chromosome, positions[0], positions[-1])

for change in CHANGES:
    # for change in base_change_map["A"]:
    print(change)
    df_change = df.loc[df["change"] == change]
    # df.to_csv("data_files\\spam.csv")
    error_rates = np.zeros(positions.size)
    for i, position in enumerate(positions):
        df_position = df_change.loc[df_change["position"] == position]
        if len(df_position.index) != 0:
            num_changes = np.sum(df_position["num changes"])
            total_consensus = np.sum(df_position["num consensus molecules"])
            error_rates[i] = num_changes / total_consensus
            # if 0 == error_rates[i]:
            #     print("zero at{}".format(position))
    ax1.plot(
        positions,
        error_rates,
        label=change,
        linestyle="None",
        marker="+",
        color=change_color_map[change],
    )


ax1.set(title=plot_title, xlabel="position", ylabel="error rate")
ax1.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
fig1.tight_layout()
# don't use standard form for axes, have to set log scale afterwards
ax1.ticklabel_format(useOffset=False, style="plain")
ax1.set(yscale="log")
fig1.savefig("plots\\{}_error_rate.png".format(plot_title))

consensuses = np.zeros(positions.size)
for i, position in enumerate(positions):
    df_position = df.loc[df["position"] == position]
    consensuses[i] = np.sum(df_position["num consensus molecules"])
ax2.plot(positions, consensuses, label="consensus", color="black")
ax2.set(title=plot_title, xlabel="position", ylabel="number of consensus molecules")
ax2.ticklabel_format(useOffset=False, style="plain")
fig2.tight_layout()
fig2.savefig("plots\\{}_total_consensuses.png".format(plot_title))

plt.show()
