import numpy as np
import matplotlib.pyplot as plt
from lookup import gene_seqs_map

tiles = np.zeros(len(gene_seqs_map))

for i in np.arange(len(tiles)):
    tiles[i] = len(gene_seqs_map[i].index)

bins = np.arange(-0.5, np.max(tiles) + 1.5)

fig, ax = plt.subplots()

ax.hist(tiles, bins=bins, color="cyan", edgecolor="black")
ax.set(title="Number of tiles per gene", xlabel="number of tiles", ylabel="frequency")
fig.savefig("plots\\gene_sequences_histogram.png")
plt.show()
