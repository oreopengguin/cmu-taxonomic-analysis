import matplotlib.pyplot as plt
import numpy as np


# -----------------------------
# Data from your HTML
# -----------------------------

seasons = ["Fall", "Winter", "Spring", "Summer"]

totals = {
    "Fall": 31149,
    "Winter": 60057,
    "Spring": 66335,
    "Summer": 52311
}


data = {
    "Bacteriodota": [
        210, 51796, 33888, 7199
    ],

    "Pseudomonadota": [
        21602, 1445, 4186, 9703
    ],

    "Pseudomonadota (low % ID)": [
        109, 2953, 5439, 2012
    ],

    "Actinomycetota": [
        92, 259, 3344, 863
    ],

    "Actinomycetota (low % ID)": [
        2860, 1231, 10169, 20699
    ],

    "Cyanobacteriota (low % ID)": [
        5892, 0, 0, 7495
    ],

    "Unresolved (low % ID)": [
        384, 2373, 9309, 4340
    ]
}


# -----------------------------
# Convert reads -> percentages
# -----------------------------

percent_data = {}

for phylum, reads in data.items():
    percent_data[phylum] = [
        (reads[i] / totals[seasons[i]]) * 100
        for i in range(len(seasons))
    ]


# -----------------------------
# Plot settings
# -----------------------------

fig, ax = plt.subplots(figsize=(9,6))


x = np.arange(len(seasons))

bottom = np.zeros(len(seasons))


# colors
colors = {
    "Bacteriodota": "#1baf7a",
    "Pseudomonadota": "#2a78d6",
    "Pseudomonadota (low % ID)": "#6da7ec",
    "Actinomycetota": "#eda100",
    "Actinomycetota (low % ID)": "#f3c259",
    "Cyanobacteriota (low % ID)": "#008300",
    "Unresolved (low % ID)": "#4a3aa7"
}


# hatch low confidence groups
low_groups = [
    "Pseudomonadota (low % ID)",
    "Actinomycetota (low % ID)",
    "Cyanobacteriota (low % ID)",
    "Unresolved (low % ID)"
]


# -----------------------------
# Create stacked bars
# -----------------------------

for phylum, values in percent_data.items():

    hatch = "//" if phylum in low_groups else None

    bars = ax.bar(
        x,
        values,
        bottom=bottom,
        label=phylum,
        color=colors[phylum],
        hatch=hatch,
        edgecolor="white",
        linewidth=1
    )


    # add labels for large sections
    for i, value in enumerate(values):
        if value >= 14:
            ax.text(
                i,
                bottom[i] + value/2,
                f"{value:.1f}%",
                ha="center",
                va="center",
                fontsize=10,
                color="white"
            )

    bottom += values



# -----------------------------
# Formatting
# -----------------------------

ax.set_ylim(0,100)

ax.set_ylabel("Relative abundance (%)")

ax.set_title(
    "Phylum-level community composition shifts by season\nHouse4 Rivers 2019 — Top 20 OTUs",
    fontsize=14,
    weight="bold"
)

ax.set_xticks(x)
ax.set_xticklabels(
    [
        f"{s}\n({totals[s]:,} reads)"
        for s in seasons
    ]
)


ax.grid(
    axis="y",
    linestyle="--",
    alpha=0.3
)

ax.set_axisbelow(True)


# legend
ax.legend(
    bbox_to_anchor=(1.05,1),
    loc="upper left",
    fontsize=9
)


plt.tight_layout()

plt.show()