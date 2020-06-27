# chr1	36931702	AT	10	10097	lane4:ACCCAGCA
# als5

import numpy as np
import matplotlib.pyplot as plt
from lookup import lookup, AGES, CHANGES, LANES, age_lane_map

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
print("something")
print(df)

# plt.savefig("plots\\" + plot_title + ".png")
# plt.show()
