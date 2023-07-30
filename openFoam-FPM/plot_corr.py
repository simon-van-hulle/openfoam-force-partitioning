import numpy as np
import matplotlib.pyplot as plt


plt.rcParams["mathtext.fontset"] = 'stix'


a = np.array([
    100,
    500,
    1000,
    2000,
    4000,
    8000,
    ])

b = np.array([
    0.7901,
    0.6668,
    0.6885,
    0.5971,
    0.5678,
    0.5660,
    ])

c = np.array([
    0.9974,
    0.9988,
    0.9993,
    0.9911,
    0.9952,
    0.9980,
    ])

plt.figure(figsize=(5.0, 2.0))
plt.scatter(a, b, marker='+', color="#ffa000", label="Moving tower")
plt.scatter(a, c, marker='o', color="#505050", label="Fixed tower")

plt.xlabel(r"$Re$")
plt.ylabel(r"$NCC$", rotation=0, labelpad=15)

plt.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.01,1))
plt.tight_layout()

plt.ylim(ymin=0)

plt.savefig("tower_corr.pdf")
plt.show()
