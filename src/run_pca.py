# run_pca.py --- Run PCA
#
# Filename: run_pca.py
# Description: Run PCA on a provided counts file
# Author: Student Zachary Maas <zama8258@colorado.edu>
# Maintainer: Student Zachary Maas <zama8258@colorado.edu>
# Created: Thu Jan 17 11:20:45 2019 (-0700)
#

# Commentary:
#
# This file takes a counts file (output from featureCounts), and
# performs PCA using scipy.stats. This data is then formatted as a
# figure and saved to the provided directory.
#

# Code:

import argparse
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

## Set up argument parsing
PARSER = argparse.ArgumentParser()
PARSER.add_argument('-c', '--counts')
PARSER.add_argument('-o', '--output')
ARGS = PARSER.parse_args()

COUNTS = ARGS.counts
OUTFILE = ARGS.output

if COUNTS is None or OUTFILE is None:
    print("Missing arguments. Try --help")
    quit()


def main():
    """
    @brief Runs PCA on the provided file 'COUNTS'

    @details Uses sklearn to perform a dimension 2 PCA on a counts
    table generated using featureCounts.

    @param None

    @return None
    """
    try:
        counts_table = pd.read_csv(COUNTS, sep='\t')
    except FileNotFoundError:
        print("File Not Found. Check the filename to '-c'")
        quit()
    counts_only = counts_table.drop(
        ["Geneid", "Chr", "Start", "End", "Strand", "Length"], axis=1).T
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(counts_only)
    principal_data = pd.DataFrame(
        data=principal_components, columns=["PC 1", "PC 2"])
    principal_data.index = list(counts_only.index)
    principal_data["Replicate"] = ["PO 1", "PO 2", "Treatment 1", "Treatment 2"]
    sns.set_style("darkgrid")
    graph = sns.scatterplot(
        x="PC 1", y="PC 2", hue="Replicate", data=principal_data)
    plt.title("Principal Components")
    # plt.show(graph)
    fig = graph.get_figure()
    fig.savefig(OUTFILE, dpi=1000)


if __name__ == "__main__":
    main()
#
# run_pca.py ends here
