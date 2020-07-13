import numpy as np
import matplotlib.pyplot as plt
from constants import (
    BASES,
    CHANGES,
    LANES,
    PEOPLE,
    age_lane_map,
    base_changes_map,
    base_color_map,
    change_color_map,
)
from lookup import lookup, seq_df, seq_data_df

# df = lookup(PEOPLE, LANES, seq_df.at[0, "chromosome"], seq_df.at[0, "start"])
sequence = 1
positions = np.arange(seq_df.at[sequence, "start"], seq_df.at[sequence, "end"])
chromosome = seq_df.at[sequence, "chromosome"]
strand = seq_df.at[sequence, "strand"]

change_error_rates_map = {}

df = seq_data_df(sequence)

base_fig_map = {}
base_axs_map = {}
fig2, ax2 = plt.subplots(figsize=(10, 7))
plot_title = "{}_{}-{}_({})".format(chromosome, positions[0], positions[-1], strand)


for base in BASES:
    base_fig_map[base], base_axs_map[base] = plt.subplots(2, 2, figsize=(15, 8))
    for change in base_changes_map[base]:
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

        change_error_rates_map[change] = error_rates

        ax2.plot(
            positions,
            error_rates,
            label=change,
            linestyle="None",
            marker="+",
            color=change_color_map[change],
        )

consensuses = np.zeros(positions.size)
for i, position in enumerate(positions):
    df_position = df.loc[df["position"] == position]
    consensuses[i] = np.sum(df_position["num consensus molecules"])


for base in BASES:
    axs = base_axs_map[base]
    axs = axs.flatten()
    for j, ax in enumerate(axs[:3]):
        for i, change in enumerate(base_changes_map[base]):
            if i == j:
                color = change_color_map[change]
                ax.set(title=change)
            else:
                color = "lightgrey"
            ax.plot(
                positions,
                change_error_rates_map[change],
                label=change,
                linestyle="None",
                marker="o",
                color=color,
                alpha=0.5,
            )
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.set(xlabel="position", ylabel="error rate", yscale="log")
    for change in base_changes_map[base]:
        axs[-1].plot(
            positions,
            change_error_rates_map[change],
            label=change,
            linestyle="None",
            marker="o",
            color=change_color_map[change],
            alpha=0.5,
        )
    axs[-1].legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
    axs[-1].set(title=base, xlabel="position", ylabel="error rate", yscale="log")

    for ax in axs:
        axtwin = ax.twinx()
        axtwin.plot(
            positions, consensuses, label="consensus", color="black", alpha=0.25
        )
        axtwin.set(ylabel="number of consensus molecules")

    fig = base_fig_map[base]
    fig.suptitle("{} {}".format(plot_title, base), size=16, y=0.52)
    fig.subplots_adjust(top=0.8)
    fig.tight_layout()
    fig.savefig("plots\\{}_{}_error_rate.png".format(plot_title, base))


ax2.set(title=plot_title, xlabel="position", ylabel="error rate")
ax2.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
fig2.tight_layout()
# don't use standard form for axes, have to set log scale afterwards
ax2.ticklabel_format(useOffset=False, style="plain")
ax2.set(yscale="log")
fig2.savefig("plots\\{}_error_rate_together.png".format(plot_title))
