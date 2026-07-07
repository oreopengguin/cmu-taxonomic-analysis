# ============================================================
# Genus-level (or highest confident taxonomy) community profile
# House4 Rivers 2019 Top OTUs
# ============================================================


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



# ============================================================
# 1. READ OTU TABLE
# ============================================================

otu_table = pd.read_csv(
    "top20_otu_table.tsv",
    sep="\t"
)


# ============================================================
# 2. ADD TAXONOMIC ASSIGNMENTS
#    Highest confident level from BLAST
# ============================================================


taxonomy_map = {

    # Bacteriodota
    "4e08bbe015e4e6fe2ec5f7d633ab287b":
        "Flectobacillaceae",

    "6f68e7660a5f7b5f44965bf1685607b2":
        "Flectobacillaceae",


    # Pseudomonadota
    "aca380d3d06d653c6296f75364e4e6fb":
        "Candidatus Pelagibacterium",


    # Pseudomonadota Erythrobacter group
    "d98dc146ecb153beb762e1cc3deff3e5":
        "Erythrobacter",


    # Actinomycetota unresolved
    "7c1dfb1e6b92aee1901909b4986b05e7":
        "Actinomycetota (unresolved genus)",


    # Flavobacterium
    "928b5af9732384984cfe2aca1d05c02f":
        "Flavobacterium",

    "cdcd9d3865dfc531fb006b7f1173637f":
        "Flavobacterium",


    # Cyanobacteriota unresolved
    "81ab86d2a43760ccffa6ef9df05fcccd":
        "Cyanobacteriota (unresolved genus)",


    # Actinomycetota Ilumatobacter group
    "66e6c667b517985b356e6fbd9a8fb2e6":
        "Actinomycetota (Ilumatobacter group)",


    # Flavobacteriaceae
    "e4bcb71092b80a020ab4a05a8873fe80":
        "Flavobacteriaceae",


    # Unresolved
    "d7dc5fbbdb32a59504730b7f2e8ac140":
        "Unresolved bacterium",


    # Limnohabitans
    "b24603aa8d4f64ab88572121e12582cb":
        "Limnohabitans",


    # Actinomycetota
    "701133fe4fa48059a020d5a0b0beee08":
        "Actinomycetota (unresolved genus)",


    # Pseudomonadota low confidence
    "d480331f53a14248512e658ad495481a":
        "Pseudomonadota (unresolved genus)",


    # Flavobacterium
    "4c0e147b36e20d2407b00fb454e66bce":
        "Flavobacterium",


    # Actinomycetota Corynebacterium
    "d6624a09afb4e1fe68a18beb5b205b2f":
        "Actinomycetota (Corynebacterium)",


    # Limnohabitans
    "2afad548b03e6a422929e8491c18c891":
        "Limnohabitans",


    # Rhodoluna
    "a7b90f4be204034a8203f2c38b49b8f8":
        "Rhodoluna",


    # Cyanobacteria unresolved
    "87343061d740c5d06ed07102233bb937":
        "Cyanobacteriota (unresolved genus)",


    # Limnohabitans
    "9a26fd7617f2cd0be3008af7e4f6e337":
        "Limnohabitans"

}



# add taxonomy column

otu_table["Taxon"] = (
    otu_table["OTU ID"]
    .map(taxonomy_map)
    .fillna("Unresolved")
)



# ============================================================
# 3. SUM OTUs INTO TAXONOMIC GROUPS
# ============================================================


season_cols = [
    "Season Reads: Fall",
    "Season Reads: Winter",
    "Season Reads: Spring",
    "Season Reads: Summer"
]


taxon_table = (
    otu_table
    .groupby("Taxon")[season_cols]
    .sum()
)



# ============================================================
# 4. KEEP TOP TAXA
# ============================================================


taxon_table["Total"] = (
    taxon_table.sum(axis=1)
)


top_taxa = (
    taxon_table
    .sort_values(
        "Total",
        ascending=False
    )
    .drop(columns="Total")
)



# optional:
# keep top 15 taxa
top_taxa = top_taxa.head(15)



# ============================================================
# 5. CONVERT TO RELATIVE ABUNDANCE
# ============================================================


relative_abundance = (
    top_taxa
    .div(
        top_taxa.sum(axis=0),
        axis=1
    )
    * 100
)



# ============================================================
# 6. STACKED BAR GRAPH
# ============================================================


fig, ax = plt.subplots(
    figsize=(11,7)
)


bottom = np.zeros(
    len(relative_abundance.columns)
)


colors = plt.cm.tab20.colors



for i, taxon in enumerate(relative_abundance.index):

    values = relative_abundance.loc[taxon]


    ax.bar(
        relative_abundance.columns,
        values,
        bottom=bottom,
        label=taxon,
        color=colors[i % len(colors)]
    )


    # label large sections
    for j, value in enumerate(values):

        if value >= 10:

            ax.text(
                j,
                bottom[j] + value/2,
                f"{value:.1f}%",
                ha="center",
                va="center",
                fontsize=9,
                color="white"
            )


    bottom += values



# ============================================================
# 7. FORMAT FIGURE
# ============================================================


ax.set_ylabel(
    "Relative abundance (%)",
    fontsize=12
)


ax.set_xlabel(
    "Season",
    fontsize=12
)


ax.set_title(
    "Genus-level (highest confident taxonomy) composition\nHouse4 Rivers 2019 Top OTUs",
    fontsize=14,
    weight="bold"
)


ax.set_ylim(
    0,
    100
)


ax.legend(
    bbox_to_anchor=(1.05,1),
    loc="upper left",
    fontsize=9
)


ax.grid(
    axis="y",
    linestyle="--",
    alpha=0.3
)


plt.xticks(
    rotation=45
)


plt.tight_layout()


plt.show()



# ============================================================
# 8. SAVE OUTPUT
# ============================================================


relative_abundance.to_csv(
    "genus_level_relative_abundance.csv"
)


plt.savefig(
    "genus_level_composition.png",
    dpi=300,
    bbox_inches="tight"
)