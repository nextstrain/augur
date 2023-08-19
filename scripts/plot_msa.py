import argparse
from Bio import AlignIO
from Bio.Align import AlignInfo
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def prepare_msa_heatmap(msa_path, consensus_threshold):
    """Plots a heatmap for the given heatmap.
    """
    msa = AlignIO.read(msa_path, "fasta")
    summary_align = AlignInfo.SummaryInfo(msa)
    consensus = summary_align.dumb_consensus(threshold=consensus_threshold)
    consensus_array = np.asarray(consensus)
    matches = np.apply_along_axis(lambda row: row == consensus_array, 1, np.asarray(msa)).astype(int)
    sorted_matches = np.array(sorted(matches, key=lambda row: row.sum(), reverse=True))

    return sorted_matches


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Finds the consensus for a given MSA and plots a heatmap with a mismatches to the consensus shown in black against white."
    )
    parser.add_argument("msa", help="FASTA file containing a multiple sequence alignment to plot")
    parser.add_argument("output", help="Binary heatmap of MSA (e.g., msa.pdf)")
    parser.add_argument("--consensus_threshold", type=float, default=0.3, help="Proportion of sequences to require for consensus")
    args = parser.parse_args()

    sorted_matches = prepare_msa_heatmap(args.msa, args.consensus_threshold)
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    axes = plt.imshow(sorted_matches, vmin=0, vmax=1, aspect="auto", cmap=plt.get_cmap("gray"))
    ax.set_xlabel("Alignment position")
    ax.set_ylabel("Strains")
    ax.set_yticklabels([])

    plt.tight_layout()
    plt.savefig(args.output)
