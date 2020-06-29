# chr1	36931702	AT	10	10097	lane4:ACCCAGCA
# als5

import numpy as np
import matplotlib.pyplot as plt
from lookup import lookup, AGES, CHANGES, LANES, age_lane_map, change_map, id_age_map

fig, ax = plt.subplots()

person = "als5"
chromosome = "chr1"
position = 36931702
plot_title = person + chromosome + str(position)

mutation_frequencies = np.zeros(AGES.size)

# for change in CHANGES:
#     print(change)
#     for i in range(AGES.size):
#         mutation_frequencies[i] = lookup_frequency(
#             id, age_lane_map[AGES[i]], chromosome, position, change
#         )

#     ax.plot(AGES, mutation_frequencies, label=change, marker="+", linestyle="None")
# ax.set(title=plot_title, xlabel="age", ylabel="frequency of mutation")
# ax.set_xticks(AGES)
# ax.legend()

df = lookup(person, LANES, chromosome, position, CHANGES)
print(df)
df["relative frequency"] = df["frequency"] / df["num consensus molecules"]
print(df)

for change in change_map["A"]:
    df_chunk = df.loc[df["change"] == change]
    ages = df_chunk["sample ID"].map(id_age_map)
    ax.bar(ages, df_chunk["relative frequency"], label=change)

ax.set(title=plot_title, yscale="log")
ax.legend()
plt.savefig("plots\\" + plot_title + ".png")
plt.show()
