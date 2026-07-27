#!/usr/bin/env python3

"""
In-silico haplotype phasing for scRNA-seq X-inactivation data.

Uses per-cell allele co-occurrence (clonal XCI means every transcript
in a cell reflects the same haplotype) instead of read-backed or
population-based phasing. See the original docstring version of this
script for the full explanation of the method.

KEY CHANGE from the first version: BAM scanning (slow) is now separate
from thresholding + phasing (fast), via --save-counts / --load-counts.
Scan the BAM once, then experiment with --min-umi-per-call,
--min-call-ratio, --min-joint-cells cheaply and repeatedly using
--load-counts, without re-reading the BAM each time.

Also: defaults are relaxed compared to the first version
(min-umi-per-call 2->1, min-call-ratio 0.8->0.7, min-joint-cells 3->1).
scRNA-seq per-site UMI depth per cell is usually 0-1, so requiring 2+
UMIs agreeing before trusting a single-cell call throws away almost
everything. With thousands of cells, many individually-noisy
single-UMI calls still carry real aggregate signal through the
spectral step -- that's the whole point of using the cell dimension.
Noisier thresholds mean noisier edges, not wrong ones on average, as
long as per-base sequencing error rate is low relative to signal.

Usage (first run -- scans BAM, caches raw counts, phases once):
  python insilico_phase.py \
      --bam possorted_genome_bam.bam \
      --vcf results/GSM7148/chrX.het.vcf.gz \
      --out results/GSM7148/GSM7148.chrX.insilico_phased.vcf.gz \
      --sample GSM7148 \
      --save-counts results/GSM7148/insilico_counts_cache

Usage (later runs -- skip BAM scan, just re-threshold/re-phase):
  python insilico_phase.py \
      --load-counts results/GSM7148/insilico_counts_cache \
      --out results/GSM7148/GSM7148.chrX.insilico_phased.v2.vcf.gz \
      --sample GSM7148 \
      --min-umi-per-call 1 --min-call-ratio 0.6 --min-joint-cells 2
"""

import argparse
import pickle
import subprocess
from collections import defaultdict

import numpy as np
import pysam
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


parser = argparse.ArgumentParser()
parser.add_argument("--bam", help="Required unless --load-counts is given")
parser.add_argument("--vcf", help="Required unless --load-counts is given. Unphased het-SNP VCF.")
parser.add_argument("--out", required=True, help="Output phased VCF path (.vcf.gz)")
parser.add_argument("--sample", required=True)
parser.add_argument("--save-counts", default=None,
                     help="Path prefix to cache raw per-cell/per-site UMI counts, for fast re-thresholding later")
parser.add_argument("--load-counts", default=None,
                     help="Path prefix previously written with --save-counts. Skips BAM scanning entirely.")
parser.add_argument("--min-umi-per-call", type=int, default=1,
                     help="Minimum UMIs supporting the majority allele at a site in a cell for that call to be trusted (default: 1)")
parser.add_argument("--min-call-ratio", type=float, default=0.7,
                     help="Minimum fraction of UMIs at a site in a cell that must agree on one allele (default: 0.7). With 1 total UMI this is trivially satisfied.")
parser.add_argument("--min-joint-cells", type=int, default=1,
                     help="Minimum number of cells jointly covering a SNP pair for that pair's edge to be trusted (default: 1)")
args = parser.parse_args()

if args.load_counts is None and (args.bam is None or args.vcf is None):
    parser.error("--bam and --vcf are required unless --load-counts is given")


# -------------------------
# 1. Get raw per-cell, per-site UMI counts: either scan the BAM, or load
#    a previously cached scan.
# -------------------------

if args.load_counts:

    print("Loading cached counts from", args.load_counts, "...")

    with open(args.load_counts + "_meta.pkl", "rb") as f:
        meta = pickle.load(f)
    sites = meta["sites"]
    cells = meta["cells"]

    ref_counts = sp.load_npz(args.load_counts + "_ref.npz")
    alt_counts = sp.load_npz(args.load_counts + "_alt.npz")

    n_sites = len(sites)
    n_cells = len(cells)
    print("Sites:", n_sites, " Cells:", n_cells)

else:

    print("Loading het SNPs...")

    vcf = pysam.VariantFile(args.vcf)

    sites = []          # list of (chrom, pos, ref, alt) in file order
    site_index = {}      # (chrom,pos) -> index into `sites`

    for rec in vcf.fetch():
        if len(rec.ref) != 1:
            continue
        if rec.alts is None or len(rec.alts) != 1 or len(rec.alts[0]) != 1:
            continue
        idx = len(sites)
        sites.append((rec.chrom, rec.pos, rec.ref, rec.alts[0]))
        site_index[(rec.chrom, rec.pos)] = idx

    n_sites = len(sites)
    print("Candidate het SNPs:", n_sites)

    print("Scanning reads for per-cell allele calls (this is the slow step)...")

    # site_idx -> cell -> {"ref": set(umi), "alt": set(umi)}
    site_cell_umis = defaultdict(lambda: defaultdict(lambda: {"ref": set(), "alt": set()}))

    bam = pysam.AlignmentFile(args.bam, "rb")

    for read in bam.fetch(until_eof=True):

        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        try:
            cell = read.get_tag("CB")
            umi = read.get_tag("UB")
        except KeyError:
            continue

        chrom = bam.get_reference_name(read.reference_id)

        for qpos, rpos in read.get_aligned_pairs():

            if rpos is None or qpos is None:
                continue

            pos = rpos + 1
            key = (chrom, pos)

            if key not in site_index:
                continue

            site_idx = site_index[key]
            ref, alt = sites[site_idx][2], sites[site_idx][3]

            base = read.query_sequence[qpos]

            if base == ref and base != alt:
                site_cell_umis[site_idx][cell]["ref"].add(umi)
            elif base == alt and base != ref:
                site_cell_umis[site_idx][cell]["alt"].add(umi)

    bam.close()

    cells = sorted({c for calls in site_cell_umis.values() for c in calls})
    cell_index = {c: i for i, c in enumerate(cells)}
    n_cells = len(cells)

    print("Cells with >=1 covered site:", n_cells)

    rows, cols, ref_vals, alt_vals = [], [], [], []
    for site_idx, cell_calls in site_cell_umis.items():
        for cell, umis in cell_calls.items():
            n_ref, n_alt = len(umis["ref"]), len(umis["alt"])
            if n_ref == 0 and n_alt == 0:
                continue
            rows.append(cell_index[cell])
            cols.append(site_idx)
            ref_vals.append(n_ref)
            alt_vals.append(n_alt)

    ref_counts = sp.csr_matrix((ref_vals, (rows, cols)), shape=(n_cells, n_sites))
    alt_counts = sp.csr_matrix((alt_vals, (rows, cols)), shape=(n_cells, n_sites))

    if args.save_counts:
        print("Caching raw counts to", args.save_counts, "...")
        sp.save_npz(args.save_counts + "_ref.npz", ref_counts)
        sp.save_npz(args.save_counts + "_alt.npz", alt_counts)
        with open(args.save_counts + "_meta.pkl", "wb") as f:
            pickle.dump({"sites": sites, "cells": cells}, f)
        print("Cached. Next time, rerun with:")
        print(f"  --load-counts {args.save_counts}")
        print("to skip the BAM scan entirely.")


# -------------------------
# 2. Threshold raw counts into confident -1/+1 calls (fast, vectorized)
# -------------------------

print("Applying thresholds: min-umi-per-call=%d min-call-ratio=%.2f"
      % (args.min_umi_per_call, args.min_call_ratio))

ref_d = np.asarray(ref_counts.todense())
alt_d = np.asarray(alt_counts.todense())
total_d = ref_d + alt_d

with np.errstate(divide="ignore", invalid="ignore"):
    ref_ratio = np.where(total_d > 0, ref_d / np.maximum(total_d, 1), 0)

is_ref_call = (ref_d >= args.min_umi_per_call) & (ref_ratio >= args.min_call_ratio)
is_alt_call = (alt_d >= args.min_umi_per_call) & ((1 - ref_ratio) >= args.min_call_ratio) & (total_d > 0)
# guard against both being true only when total_d==0 edge case handled above

call = np.zeros_like(total_d, dtype=np.int8)
call[is_ref_call & ~is_alt_call] = -1
call[is_alt_call & ~is_ref_call] = 1

M = sp.csr_matrix(call)
n_calls = M.nnz
print("Confident (cell, site) calls:", n_calls, "out of", ref_counts.nnz + alt_counts.nnz, "raw observations")


# -------------------------
# 3. Signed SNP x SNP co-occurrence matrix and spectral 2-clustering
# -------------------------

print("Computing SNP-SNP co-occurrence and phasing...")

S = (M.T @ M).astype(float)
S = S.tolil()
S.setdiag(0)
S = S.tocsr()

Mabs = M.copy()
Mabs.data = np.abs(Mabs.data)
joint_coverage = (Mabs.T @ Mabs).tolil()
joint_coverage.setdiag(0)
joint_coverage = joint_coverage.tocsr()

S_arr = S.toarray()
cov_arr = joint_coverage.toarray()
S_arr[cov_arr < args.min_joint_cells] = 0
S = sp.csr_matrix(S_arr)

has_edge = np.asarray((np.abs(S) > 0).sum(axis=1)).flatten() > 0
phaseable_idx = np.where(has_edge)[0]

print("Sites with usable linkage:", len(phaseable_idx), "/", n_sites,
      f"(min-joint-cells={args.min_joint_cells})")

if len(phaseable_idx) < 2:
    raise SystemExit(
        "Not enough linked SNPs to phase (need >=2 sites with joint "
        "cell coverage). Try lowering --min-joint-cells / "
        "--min-umi-per-call further, or check that --bam / --vcf match."
    )

S_sub = S[phaseable_idx][:, phaseable_idx]

if S_sub.shape[0] == 2:
    eigvec = np.array([1.0, -1.0 if S_sub[0, 1] < 0 else 1.0])
else:
    _, eigvecs = eigsh(S_sub, k=1, which="LA")
    eigvec = eigvecs[:, 0]

hap_label = {}
for local_i, site_idx in enumerate(phaseable_idx):
    hap_label[site_idx] = 1 if eigvec[local_i] >= 0 else 0


# -------------------------
# 4. Write phased VCF
# -------------------------

print("Writing phased VCF...")

tmp_vcf = args.out.replace(".vcf.gz", "") + ".tmp.vcf"

with open(tmp_vcf, "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    contigs_seen = sorted({sites[i][0] for i in hap_label.keys()})
    for contig in contigs_seen:
        f.write(f"##contig=<ID={contig}>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
        args.sample + "\n"
    )

    n_written = 0
    for site_idx in sorted(hap_label.keys(), key=lambda i: (sites[i][0], sites[i][1])):
        chrom, pos, ref, alt = sites[site_idx]
        hap1_is_alt = (hap_label[site_idx] == 1)
        gt = "0|1" if hap1_is_alt else "1|0"
        f.write(
            f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"
        )
        n_written += 1

print("Phased sites written:", n_written)

subprocess.run(["bcftools", "sort", tmp_vcf, "-Oz", "-o", args.out], check=True)
subprocess.run(["bcftools", "index", "-t", args.out], check=True)
subprocess.run(["rm", tmp_vcf], check=True)

print("Done:", args.out)
