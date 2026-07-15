#!/usr/bin/env python3

import argparse
import sys
import pysam
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from scipy.stats import binomtest

#python xci_analysis.py   --bam GSM7148/possorted_genome_bam.bam  --vcf /GSM7148/chrX.phased.vcf.gz   --gtf /GSM7148/gencode.v44.annotation.gtf


# =========================================================
# CLI ARGUMENTS
# =========================================================
parser = argparse.ArgumentParser(description="Haplotype-based XCI Escape Detection Pipeline")

parser.add_argument("--bam", required=True)
parser.add_argument("--vcf", required=True)
parser.add_argument("--gtf", required=True)

parser.add_argument("--min_mapq", type=int, default=10)
parser.add_argument("--min_cells_per_gene", type=int, default=5)
parser.add_argument("--min_hap_reads", type=int, default=2)

# NEW: SNP control
parser.add_argument("--min_snps_per_gene", type=int, default=3)

args = parser.parse_args()

# =========================================================
# HELPERS
# =========================================================
def norm_chrom(c):
    return str(c).replace("chr", "")

def get_cb(read):
    return read.get_tag("CB") if read.has_tag("CB") else None

def get_umi(read):
    return read.get_tag("UB") if read.has_tag("UB") else None

def get_base(read, pos):
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if rpos == pos - 1:
            if qpos is None:
                return None
            return read.query_sequence[qpos].upper()
    return None

# =========================================================
# LOAD SNPs
# =========================================================
def load_phased_snps(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    snps = {}

    for rec in vcf.fetch():
        if len(rec.alts) != 1 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue

        gt = rec.samples[0]["GT"]
        if gt is None or None in gt:
            continue
        if set(gt) != {0, 1}:
            continue

        phased = rec.samples[0].phased if hasattr(rec.samples[0], "phased") else False
        if not phased:
            continue

        chrom = norm_chrom(rec.contig)

        snps[(chrom, rec.pos)] = {
            "ref": rec.ref.upper(),
            "alt": rec.alts[0].upper()
        }

    print(f"[✓] Loaded phased SNPs: {len(snps)}")
    return snps

# =========================================================
# LOAD GENES
# =========================================================
def load_genes(gtf_file):
    gr = pr.read_gtf(gtf_file)
    genes = gr.df
    genes = genes[genes.Feature == "gene"].copy()
    genes["Chromosome"] = genes["Chromosome"].apply(norm_chrom)
    return genes[["Chromosome", "Start", "End", "gene_name"]]

def map_gene(chrom, pos, genes):
    hits = genes[
        (genes["Chromosome"] == chrom) &
        (genes["Start"] <= pos) &
        (genes["End"] >= pos)
    ]
    if hits.empty:
        return None
    return hits.iloc[0]["gene_name"]

# =========================================================
# BUILD MATRIX (NOW WITH SNP ID)
# =========================================================
def build_hap_ase(bam_file, snps, genes):

    bam = pysam.AlignmentFile(bam_file, "rb")
    records = []

    print("[-] Scanning BAM...")

    for read in tqdm(bam.fetch(until_eof=True)):

        if read.is_unmapped or read.mapping_quality < args.min_mapq:
            continue

        cell = get_cb(read)
        umi = get_umi(read)
        if cell is None or umi is None:
            continue

        chrom = norm_chrom(bam.get_reference_name(read.reference_id))

        for pos0 in read.get_reference_positions():
            pos = pos0 + 1
            key = (chrom, pos)

            if key not in snps:
                continue

            base = get_base(read, pos)
            if base is None:
                continue

            snp = snps[key]
            gene = map_gene(chrom, pos, genes)
            if gene is None:
                continue

            if base == snp["ref"]:
                hap = 0
            elif base == snp["alt"]:
                hap = 1
            else:
                continue

            snp_id = f"{chrom}:{pos}"

            records.append((cell, umi, gene, hap, snp_id))

    bam.close()

    df = pd.DataFrame(records,
                      columns=["cell", "umi", "gene", "hap", "snp"])

    print(f"[✓] SNP overlaps: {len(df)}")

    return df

# =========================================================
# WEIGHTED UMI + SNP AGGREGATION
# =========================================================
def aggregate(df):

    if df.empty:
        return pd.DataFrame()

    print("[-] Weighted SNP aggregation...")

    # UMI consensus
    umi = (
        df.groupby(["cell", "umi", "gene", "snp"])["hap"]
        .agg(lambda x: x.mode()[0])
        .reset_index()
    )

    # SNP weights (informativeness)
    snp_weight = umi.groupby(["gene", "snp"]).size().reset_index(name="weight")

    umi = umi.merge(snp_weight, on=["gene", "snp"], how="left")

    out = []

    for (gene, cell), sub in umi.groupby(["gene", "cell"]):

        hap1 = sub[sub["hap"] == 1]["weight"].sum()
        hap0 = sub[sub["hap"] == 0]["weight"].sum()

        total = hap1 + hap0

        if total < args.min_hap_reads:
            continue

        out.append((gene, cell, hap1, hap0, total, hap1 / total))

    return pd.DataFrame(out,
                        columns=["gene", "cell", "hap1", "hap0", "total", "ratio"])

# =========================================================
# ESCAPE STATS (WITH SNP FILTER)
# =========================================================
def compute_escape(gdf):

    if gdf.empty:
        return pd.DataFrame()

    results = []

    for gene, sub in gdf.groupby("gene"):

        if len(sub) < args.min_cells_per_gene:
            continue

        # NEW: SNP-support proxy
        snp_support = sub["total"].sum()

        if snp_support < args.min_snps_per_gene:
            continue

        hap1 = sub["hap1"].sum()
        hap0 = sub["hap0"].sum()
        total = hap1 + hap0

        if total == 0:
            continue

        p_hat = hap1 / total
        pval = binomtest(hap1, total, p=0.5).pvalue
        escape_score = 1.0 - abs(p_hat - 0.5)

        results.append({
            "gene": gene,
            "n_cells": len(sub),
            "hap1_total": hap1,
            "hap0_total": hap0,
            "p_hat": p_hat,
            "escape_score": escape_score,
            "p_value": pval,
            "snp_support": snp_support
        })

    return pd.DataFrame(results)

# =========================================================
# MAIN
# =========================================================
def main():

    print("=========================================================")
    print("XCI HAPLOTYPE PIPELINE (WEIGHTED SNP VERSION)")
    print("=========================================================")

    snps = load_phased_snps(args.vcf)
    genes = load_genes(args.gtf)

    df = build_hap_ase(args.bam, snps, genes)
    df.to_csv("ase1.csv", index=False)

    if df.empty:
        print("No signal")
        sys.exit(1)

    gdf = aggregate(df)
    gdf.to_csv("agg1.csv", index=False)

    stats = compute_escape(gdf)

    if stats.empty:
        print("No genes passed filters")
        sys.exit(1)

    stats = stats.sort_values("escape_score", ascending=False)
    stats.to_csv("xci_haplotype_escape1.tsv", sep="\t", index=False)

    print(stats.head(20).to_string(index=False))

if __name__ == "__main__":
    main()