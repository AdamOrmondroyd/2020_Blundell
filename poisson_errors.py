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
base = "A"

df = lookup(chromosome, positions)

fig1, axs = plt.subplots(2, 2, figsize=(15, 7))
fig2, ax2 = plt.subplots(figsize=(10, 7))
fig3, ax3 = plt.subplots(figsize=(10, 7))
plot_title = "{}_{}-{}".format(chromosome, positions[0], positions[-1])


for base, ax in zip(np.array(["A", "C", "G", "T"]), axs.flatten()):
    for change in base_change_map[base]:
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
        ax.plot(
            positions,
            error_rates,
            label=change,
            linestyle="None",
            marker="+",
            color=change_color_map[change],
        )
        ax2.plot(
            positions,
            error_rates,
            label=change,
            linestyle="None",
            marker="+",
            color=change_color_map[change],
        )
    ax.set(title=base, xlabel="position", ylabel="error rate")
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")

    # don't use standard form for axes, have to set log scale afterwards
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set(yscale="log")

fig1.tight_layout()
ax2.set(title=plot_title, xlabel="position", ylabel="error rate")
ax2.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
fig2.tight_layout()
# don't use standard form for axes, have to set log scale afterwards
ax2.ticklabel_format(useOffset=False, style="plain")
ax2.set(yscale="log")
fig1.savefig("plots\\{}_error_rate.png".format(plot_title))
fig2.savefig("plots\\{}_error_rate_together.png".format(plot_title))

consensuses = np.zeros(positions.size)
for i, position in enumerate(positions):
    df_position = df.loc[df["position"] == position]
    consensuses[i] = np.sum(df_position["num consensus molecules"])
ax3.plot(positions, consensuses, label="consensus", color="black")
ax3.set(title=plot_title, xlabel="position", ylabel="number of consensus molecules")
ax3.ticklabel_format(useOffset=False, style="plain")
fig3.tight_layout()
fig3.savefig("plots\\{}_total_consensuses.png".format(plot_title))
ax3.set(yscale="log")
fig3.savefig("plots\\{}_total_consensuses_log.png".format(plot_title))

plt.show()
