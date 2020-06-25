# chr1	36931702	AT	10	10097	lane4:ACCCAGCA
# als5

import numpy as np
import matplotlib.pyplot as plt
from lookup import lookup_frequency, AGES, CHANGES, age_lane_map

fig, ax = plt.subplots()

id = "als5"
chromosome = "chr1"
position = 36931702
plot_title = id + chromosome + str(position)

mutation_frequencies = np.zeros(AGES.size)

for change in CHANGES:
    print(change)
    for i in range(AGES.size):
        mutation_frequencies[i] = lookup_frequency(
            id, age_lane_map[AGES[i]], chromosome, position, change
        )

    ax.plot(AGES, mutation_frequencies, label=change, marker="+", linestyle="None")
ax.set(title=plot_title, xlabel="age", ylabel="frequency of mutation")
ax.set_xticks(AGES)
ax.legend()

plt.savefig("plots\\" + plot_title + ".png")
plt.show()
