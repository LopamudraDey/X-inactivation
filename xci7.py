#!/usr/bin/env python3

import argparse
import sys
import pysam
import pandas as pd
import pyranges as pr
from tqdm import tqdm
from scipy.stats import binomtest

# =========================================================
# CLI ARGUMENTS
# =========================================================
parser = argparse.ArgumentParser(description="Haplotype-based XCI Escape Detection Pipeline")

parser.add_argument("--bam", required=True, help="Input single-cell BAM file (with CB and UB/UMI tags)")
parser.add_argument("--vcf", required=True, help="Phased, single-sample VCF file containing heterozygous SNPs")
parser.add_argument("--gtf", required=True, help="Gene annotation GTF file matching your reference genome")

parser.add_argument("--min_mapq", type=int, default=10, help="Minimum mapping quality for reads")
parser.add_argument("--min_cells_per_gene", type=int, default=5, help="Minimum qualifying cells needed to evaluate a gene")
parser.add_argument("--min_hap_reads", type=int, default=1, help="Minimum unique UMIs needed per cell-gene pair (set to 1 for sparse data)")

args = parser.parse_args()

# =========================================================
# HELPER FUNCTIONS
# =========================================================
def norm_chrom(c):
    """Standardize chromosome naming conventions (strips 'chr')."""
    return str(c).replace("chr", "")

def get_cb(read):
    """Extract Cell Barcode tag from 10x Genomics read alignment."""
    return read.get_tag("CB") if read.has_tag("CB") else None

def get_umi(read):
    """Extract Unique Molecular Identifier (UB) tag from 10x Genomics read alignment."""
    return read.get_tag("UB") if read.has_tag("UB") else None

def get_base(read, pos):
    """Finds the precise nucleotide base matching the reference genomic coordinate."""
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if rpos == pos - 1:
            if qpos is None:
                return None
            return read.query_sequence[qpos].upper()
    return None

# =========================================================
# DATA LOADING FUNCTIONS
# =========================================================
def load_phased_snps(vcf_file):
    """Loads high-confidence phased, heterozygous single-nucleotide variants."""
    vcf = pysam.VariantFile(vcf_file)
    snps = {}

    for rec in vcf.fetch():
        # Restrict strictly to simple biallelic single nucleotide variants
        if len(rec.alts) != 1 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue

        gt = rec.samples[0]["GT"]
        if gt is None or None in gt:
            continue
        if set(gt) != {0, 1}:  # Must be strictly heterozygous
            continue

        # Enforce phasing constraint to maintain parental tracking matrix stability
        phased = rec.samples[0].phased if hasattr(rec.samples[0], "phased") else False
        if not phased:
            continue

        hp = rec.samples[0].get("HP", None)
        chrom = norm_chrom(rec.contig)

        snps[(chrom, rec.pos)] = {
            "ref": rec.ref.upper(),
            "alt": rec.alts[0].upper(),
            "hp": hp
        }

    print(f"[✓] Loaded phased SNPs: {len(snps)}")
    return snps

def load_genes(gtf_file):
    """Parses structural genomic boundaries from a reference annotation file."""
    print("[-] Parsing GTF boundaries...")
    gr = pr.read_gtf(gtf_file)
    genes = gr.df
    genes = genes[genes.Feature == "gene"].copy()
    genes["Chromosome"] = genes["Chromosome"].apply(norm_chrom)
    return genes[["Chromosome", "Start", "End", "gene_name"]]

def map_gene(chrom, pos, genes):
    """Vectorized coordinate search to locate intersecting gene boundaries."""
    hits = genes[
        (genes["Chromosome"] == chrom) &
        (genes["Start"] <= pos) &
        (genes["End"] >= pos)
    ]
    if hits.empty:
        return None
    return hits.iloc[0]["gene_name"]

# =========================================================
# CORE PROCESSING PIPELINE
# =========================================================
def build_hap_ase(bam_file, snps, genes):
    """Scans the alignment landscape to match reads against phased SNPs and gene loci."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    records = []

    print("[-] Scanning BAM file alignments...")
    for read in tqdm(bam.fetch(until_eof=True), desc="Processing Reads"):
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

            # Establish relative Haplotype alignment values
            if base == snp["ref"]:
                hap = 0
            elif base == snp["alt"]:
                hap = 1
            else:
                continue

            records.append((cell, umi, gene, hap))

    bam.close()
    df = pd.DataFrame(records, columns=["cell", "umi", "gene", "hap"])
    
    # --- DIAGNOSTIC METRICS ---
    print(f"\n[✓] Raw SNP-overlapping alignments tracked: {df.shape[0]}")
    if not df.empty:
        print(f"    - Unique genes identified: {df['gene'].nunique()}")
        print(f"    - Unique cells identified: {df['cell'].nunique()}")
    
    return df

# =========================================================
# UMI CONSENSUS AGGREGATION
# =========================================================
def aggregate(df):
    """Resolves multiple read fragments per molecule to a single consensus allele."""
    if df.empty:
        return pd.DataFrame()

    print("[-] Resolving unique UMI molecules to consensus haplotypes...")
    # Group by the unique molecular tag and compute the modal consensus haplotype.
    # This prevents unphased multi-SNP transcripts from forcing an artificial 50:50 score.
    umi_consensus = df.groupby(["cell", "umi", "gene"])["hap"].agg(lambda x: x.mode()[0]).reset_index()
    print(f"[✓] Deduplicated molecular landscape resolved to {len(umi_consensus)} unique UMI observations.")

    out = []
    print("[-] Collapsing molecular indices to tissue profiles...")
    for (gene, cell), sub in umi_consensus.groupby(["gene", "cell"]):
        hap1 = int((sub["hap"] == 1).sum())
        hap0 = int((sub["hap"] == 0).sum())
        total = hap0 + hap1

        if total < args.min_hap_reads:
            continue

        ratio = hap1 / total
        out.append((gene, cell, hap1, hap0, total, ratio))

    agg_df = pd.DataFrame(out, columns=["gene", "cell", "hap1", "hap0", "total", "ratio"])
    print(f"[✓] Complete. {len(agg_df)} cell-gene expression profiles survived processing thresholds.")
    return agg_df

# =========================================================
# ESCAPE ANALYSIS METRICS
# =========================================================
def compute_escape(gdf):
    """Calculates global allele distributions and executes statistical significance evaluations."""
    if gdf.empty:
        return pd.DataFrame()
        
    results = []
    for gene, sub in gdf.groupby("gene"):
        # Enforce minimum cellular population representation parameters
        if len(sub) < args.min_cells_per_gene:
            continue

        hap1 = int(sub["hap1"].sum())
        hap0 = int(sub["hap0"].sum())
        total = hap1 + hap0
        
        if total == 0:
            continue

        p_hat = hap1 / total
        pval = binomtest(hap1, total, p=0.5).pvalue
        
        # Calculate deviation from structural allele parity 
        escape_score = 1.0 - abs(p_hat - 0.5)

        results.append({
            "gene": gene,
            "n_cells": len(sub),
            "hap1_total": hap1,
            "hap0_total": hap0,
            "p_hat": p_hat,
            "escape_score": escape_score,
            "p_value": pval
        })

    return pd.DataFrame(results)

# =========================================================
# PIPELINE EXECUTION FLOW
# =========================================================
def main():
    print("=========================================================")
    print("STARTING HAPLOTYPE ESCAPE DETECTION PIPELINE")
    print("=========================================================")
    
    snps = load_phased_snps(args.vcf)
    if not snps:
        print("[X] ERROR: Phased VCF matrix is empty or invalid. Terminating workflow.")
        sys.exit(1)
        
    genes = load_genes(args.gtf)
    if genes.empty:
        print("[X] ERROR: GTF file produced an empty structural database. Terminating workflow.")
        sys.exit(1)

    df = build_hap_ase(args.bam, snps, genes)
    df.to_csv("ase.csv", index=False)
    print("[✓] Raw mapping configurations cached: ase.csv")

    if df.empty:
        print("[X] ERROR: Zero shared intersections found between BAM, VCF, and GTF coordinates.")
        print("    Check chromosome naming alignment (e.g., 'chrX' vs 'X') across files.")
        sys.exit(1)

    gdf = aggregate(df)
    gdf.to_csv("agg.csv", index=False)
    print("[✓] High-confidence cell-level profiles cached: agg.csv")

    if gdf.empty:
        print("\n[X] FILTER DROP: No cell-gene pairs passed your count limits.")
        print(f"    Try lowering your --min_hap_reads threshold (Current value: {args.min_hap_reads}).")
        sys.exit(0)

    print("[-] Quantifying molecular skewing and analytical metrics...")
    stats = compute_escape(gdf)

    if stats.empty:
        print("\n[X] FILTER DROP: Data was dropped during gene-level grouping.")
        print(f"    Your genes lacked enough qualifying cells. Try lowering --min_cells_per_gene (Current value: {args.min_cells_per_gene}).")
        sys.exit(0)

    # Sort results to highlight strong biallelic escape expressions first
    stats = stats.sort_values("escape_score", ascending=False)
    stats.to_csv("xci_haplotype_escape.tsv", sep="\t", index=False)
    print("[✓] Output analysis metrics archived: xci_haplotype_escape.tsv")

    print("\n=========================================================")
    print("TOP PREDICTED GENES:")
    print("=========================================================")
    print(stats.head(20).to_string(index=False))
    print("=========================================================\n[✓] PIPELINE RUN COMPLETED SUCCESSFULLY.")

if __name__ == "__main__":
    main()