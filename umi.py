#!/usr/bin/env python3

import pysam
import pandas as pd
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()

parser.add_argument("--bam", required=True)
parser.add_argument("--vcf", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument(
    "--min-umi",
    type=int,
    default=5,
    help="Minimum total (Hap0+Hap1) UMI count for a cell to be included "
         "in the per-cell fraction plot (default: 5)."
)
parser.add_argument(
    "--top-n-genes",
    type=int,
    default=25,
    help="Number of top genes (by total UMI) to show in the per-gene plot "
         "(default: 25)."
)

args = parser.parse_args()

import os
os.makedirs(args.out, exist_ok=True)


print("Loading phased SNPs...")


# -------------------------
# Read phased VCF
# -------------------------
# We store REF/ALT identity per haplotype (not just the raw GT), so that
# downstream we can compare the actual base observed in each read against
# the allele expected on each haplotype. Using GT alone (as in the
# original script) ignores what base the read actually carries and just
# labels every overlapping read as both Hap0 and Hap1 evidence, which is
# not allele-specific counting.

vcf = pysam.VariantFile(args.vcf)

variants = {}

for rec in vcf.fetch():

    if not rec.samples[0].phased:
        continue

    gt = rec.samples[0]["GT"]

    # Only handle biallelic SNPs cleanly; skip anything else (multiallelic,
    # indels slipped through, missing alleles) to avoid mis-assigning reads.
    if len(rec.ref) != 1:
        continue
    if rec.alts is None or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
        continue
    if gt[0] is None or gt[1] is None:
        continue
    if gt[0] == gt[1]:
        # Homozygous phased sites carry no haplotype information for ASE.
        continue

    ref = rec.ref
    alt = rec.alts[0]

    variants[(rec.chrom, rec.pos)] = {
        "hap0_allele": ref if gt[0] == 0 else alt,
        "hap1_allele": ref if gt[1] == 0 else alt,
    }


print(
    "Informative phased SNPs:",
    len(variants)
)


# -------------------------
# Open BAM
# -------------------------

bam = pysam.AlignmentFile(
    args.bam,
    "rb"
)


cell_counts = defaultdict(
    lambda:
    {
        "Hap0": set(),
        "Hap1": set()
    }
)


gene_counts = defaultdict(
    lambda:
    {
        "Hap0": set(),
        "Hap1": set()
    }
)


print("Scanning reads...")


for read in bam.fetch(
    until_eof=True
):

    if read.is_unmapped:
        continue

    # Skip secondary/supplementary alignments to avoid double-counting the
    # same UMI/molecule from multi-mapped or split alignments.
    if read.is_secondary or read.is_supplementary:
        continue

    # Cell barcode / UMI
    try:
        cell = read.get_tag("CB")
        umi = read.get_tag("UB")
    except KeyError:
        continue

    gene = "Unknown"

    if read.has_tag("GN"):
        gene = read.get_tag("GN")

    chrom = bam.get_reference_name(
        read.reference_id
    )

    for qpos, rpos in read.get_aligned_pairs():

        # rpos is None -> insertion relative to reference (no ref position)
        # qpos is None -> deletion relative to reference (no read base)
        # Both must be present for us to read an actual base at a real
        # reference position.
        if rpos is None or qpos is None:
            continue

        pos = rpos + 1

        key = (chrom, pos)

        if key not in variants:
            continue

        base = read.query_sequence[qpos]

        v = variants[key]

        hap0_allele = v["hap0_allele"]
        hap1_allele = v["hap1_allele"]

        # Assign based on which haplotype's allele the observed base
        # actually matches. A base matching neither (sequencing error /
        # third allele) or matching both (shouldn't happen at a real het
        # site) is dropped rather than guessed at.
        if base == hap0_allele and base != hap1_allele:
            cell_counts[cell]["Hap0"].add(umi)
            gene_counts[gene]["Hap0"].add(umi)

        elif base == hap1_allele and base != hap0_allele:
            cell_counts[cell]["Hap1"].add(umi)
            gene_counts[gene]["Hap1"].add(umi)


bam.close()


# -------------------------
# Save cell counts
# -------------------------

cells = []

for c, v in cell_counts.items():

    cells.append(
        [
            c,
            len(v["Hap0"]),
            len(v["Hap1"])
        ]
    )


cell_df = pd.DataFrame(
    cells,
    columns=[
        "Cell",
        "Hap0_UMI",
        "Hap1_UMI"
    ]
)

cell_df["Total_UMI"] = cell_df["Hap0_UMI"] + cell_df["Hap1_UMI"]

# Avoid div-by-zero for cells with no allele-specific UMIs at all
cell_df["Hap0_fraction"] = cell_df["Hap0_UMI"] / cell_df["Total_UMI"].replace(0, pd.NA)

cell_df.to_csv(
    args.out +
    "/haplotype_counts_per_cell.csv",
    index=False
)


# -------------------------
# Allele-specific UMI plot (totals)
# -------------------------

total_hap0 = cell_df["Hap0_UMI"].sum()
total_hap1 = cell_df["Hap1_UMI"].sum()

labels = [
    "Haplotype 0",
    "Haplotype 1"
]

values = [
    total_hap0,
    total_hap1
]

plt.figure(figsize=(5, 4))

plt.bar(
    labels,
    values
)

plt.ylabel(
    "Total allele-specific UMI counts"
)

plt.title(
    args.sample +
    " X-linked allele-specific expression"
)

for i, v in enumerate(values):
    plt.text(
        i,
        v,
        str(v),
        ha="center",
        va="bottom"
    )

plt.tight_layout()

plt.savefig(
    args.out +
    "/Hap0_vs_Hap1_UMI_barplot.png",
    dpi=300
)

plt.close()

print(
    "Saved:",
    args.out + "/Hap0_vs_Hap1_UMI_barplot.png"
)


# -------------------------
# Gene counts
# -------------------------

genes = []

for g, v in gene_counts.items():

    genes.append(
        [
            g,
            len(v["Hap0"]),
            len(v["Hap1"])
        ]
    )

gene_df = pd.DataFrame(
    genes,
    columns=[
        "Gene",
        "Hap0_UMI",
        "Hap1_UMI"
    ]
)

gene_df["Total_UMI"] = gene_df["Hap0_UMI"] + gene_df["Hap1_UMI"]
gene_df["Hap0_fraction"] = gene_df["Hap0_UMI"] / gene_df["Total_UMI"].replace(0, pd.NA)

gene_df.to_csv(
    args.out +
    "/haplotype_counts_per_gene.csv",
    index=False
)


# -------------------------
# Per-gene Hap0 vs Hap1 plot
# -------------------------
# Stacked horizontal bar of the top-N genes by total allele-specific UMI
# count, showing the Hap0/Hap1 split for each. This is the standard way to
# spot genes that escape X-inactivation (both haplotypes ~equally
# expressed) versus genes that are cleanly monoallelic (skewed to one
# haplotype), assuming clonal skewing of inactivation in the population/
# sample profiled.

gene_plot_df = (
    gene_df[gene_df["Gene"] != "Unknown"]
    .sort_values("Total_UMI", ascending=False)
    .head(args.top_n_genes)
    .sort_values("Total_UMI", ascending=True)  # so largest ends up on top
)

if len(gene_plot_df) > 0:

    fig_height = max(4, 0.35 * len(gene_plot_df))
    plt.figure(figsize=(7, fig_height))

    y_pos = range(len(gene_plot_df))

    plt.barh(
        y_pos,
        gene_plot_df["Hap0_UMI"],
        label="Haplotype 0"
    )

    plt.barh(
        y_pos,
        gene_plot_df["Hap1_UMI"],
        left=gene_plot_df["Hap0_UMI"],
        label="Haplotype 1"
    )

    plt.yticks(y_pos, gene_plot_df["Gene"])
    plt.xlabel("Allele-specific UMI counts")
    plt.title(
        args.sample +
        f" top {len(gene_plot_df)} genes: Hap0 vs Hap1 UMI counts"
    )
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        args.out +
        "/Hap0_vs_Hap1_per_gene.png",
        dpi=300
    )

    plt.close()

    print(
        "Saved:",
        args.out + "/Hap0_vs_Hap1_per_gene.png"
    )

else:
    print(
        "Skipped per-gene plot: no genes with allele-specific UMI counts."
    )


# -------------------------
# Per-cell Hap0 fraction plot
# -------------------------
# Histogram of each cell's Hap0 fraction (Hap0_UMI / Total_UMI), restricted
# to cells with at least --min-umi total allele-specific UMIs so that noisy,
# low-count cells don't dominate the distribution. For clonal X-inactivation
# you expect a bimodal distribution near 0 and 1 (cells committed to one
# haplotype); skewing toward the middle suggests escape from inactivation,
# doublets, or ambient/index contamination.

cell_plot_df = cell_df[cell_df["Total_UMI"] >= args.min_umi].dropna(subset=["Hap0_fraction"])

if len(cell_plot_df) > 0:

    plt.figure(figsize=(6, 4))

    plt.hist(
        cell_plot_df["Hap0_fraction"],
        bins=20,
        range=(0, 1),
        edgecolor="black"
    )

    plt.xlabel("Hap0 fraction (Hap0_UMI / Total_UMI)")
    plt.ylabel("Number of cells")
    plt.title(
        args.sample +
        f" per-cell Hap0 fraction (n={len(cell_plot_df)} cells, "
        f"min {args.min_umi} UMI)"
    )
    plt.tight_layout()

    plt.savefig(
        args.out +
        "/Hap0_vs_Hap1_per_cell.png",
        dpi=300
    )

    plt.close()

    print(
        "Saved:",
        args.out + "/Hap0_vs_Hap1_per_cell.png"
    )

else:
    print(
        f"Skipped per-cell plot: no cells with >= {args.min_umi} "
        "total allele-specific UMIs."
    )


print("Finished")
