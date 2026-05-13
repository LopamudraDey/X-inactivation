import pysam
import pandas as pd
import numpy as np
from collections import defaultdict
from tqdm import tqdm

from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

import umap.umap_ as umap
import matplotlib.pyplot as plt


# =========================================================
# INPUTS
# =========================================================

BAM = "clean.bam"
VCF = "chrX.phased.vcf.gz"

MIN_MAPQ = 10
MIN_CELL_SNPS = 5

N_CLUSTERS = 2


# =========================================================
# 1. LOAD PHASED SNPs
# =========================================================

def load_phased_snps(vcf_file):

    vcf = pysam.VariantFile(vcf_file)

    snps = {}

    for rec in vcf.fetch():

        if len(rec.alts) != 1:
            continue

        gt = rec.samples[0]["GT"]

        if gt is None or None in gt:
            continue

        # phased only
        if gt not in [(0,1), (1,0)]:
            continue

        snps[(rec.contig, rec.pos)] = {
            "gt": gt,
            "ref": rec.ref,
            "alt": rec.alts[0]
        }

    return snps


# =========================================================
# 2. GET READ BASE AT SNP
# =========================================================

def get_base_at_pos(read, target_pos):

    for qpos, rpos in read.get_aligned_pairs(matches_only=True):

        if rpos == target_pos - 1:

            if qpos is None:
                return None

            return read.query_sequence[qpos]

    return None


# =========================================================
# 3. BUILD CELL × SNP MATRIX
# =========================================================

def build_ase_matrix(bam_file, snps):

    bam = pysam.AlignmentFile(bam_file, "rb")

    cell_snp = defaultdict(dict)

    for read in tqdm(bam.fetch(until_eof=True)):

        if read.is_unmapped:
            continue

        if read.mapping_quality < MIN_MAPQ:
            continue

        if not read.has_tag("CB"):
            continue

        cell = read.get_tag("CB")

        chrom = bam.get_reference_name(read.reference_id)

        positions = read.get_reference_positions()

        for pos0 in positions:

            pos = pos0 + 1

            key = (chrom, pos)

            if key not in snps:
                continue

            snp = snps[key]

            base = get_base_at_pos(read, pos)

            if base is None:
                continue

            ref = snp["ref"]
            alt = snp["alt"]
            gt = snp["gt"]

            # -----------------------------------------
            # Assign allele to phased haplotype
            # -----------------------------------------

            value = None

            if base == ref:

                if gt == (0,1):
                    value = +1
                elif gt == (1,0):
                    value = -1

            elif base == alt:

                if gt == (0,1):
                    value = -1
                elif gt == (1,0):
                    value = +1

            if value is None:
                continue

            snp_id = f"{chrom}:{pos}"

            # overwrite allowed
            cell_snp[cell][snp_id] = value

    return cell_snp


# =========================================================
# 4. CONVERT TO DATAFRAME
# =========================================================

def build_dataframe(cell_snp):

    all_snps = set()

    for d in cell_snp.values():
        all_snps.update(d.keys())

    all_snps = sorted(all_snps)

    rows = []
    cells = []

    for cell, d in cell_snp.items():

        if len(d) < MIN_CELL_SNPS:
            continue

        row = []

        for snp in all_snps:
            row.append(d.get(snp, np.nan))

        rows.append(row)
        cells.append(cell)

    df = pd.DataFrame(
        rows,
        index=cells,
        columns=all_snps
    )

    return df


# =========================================================
# 5. PCA + UMAP
# =========================================================

def run_embedding(df):

    X = df.values

    # fill missing SNPs with 0
    imp = SimpleImputer(strategy="constant", fill_value=0)
    X = imp.fit_transform(X)

    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    pca = PCA(n_components=20)
    X_pca = pca.fit_transform(X)

    reducer = umap.UMAP(
        n_neighbors=15,
        min_dist=0.3,
        metric="cosine",
        random_state=42
    )

    embedding = reducer.fit_transform(X_pca)

    return embedding


# =========================================================
# 6. CLUSTERING
# =========================================================

def cluster_cells(embedding):

    kmeans = KMeans(
        n_clusters=N_CLUSTERS,
        random_state=42,
        n_init="auto"
    )

    labels = kmeans.fit_predict(embedding)

    return labels


# =========================================================
# 7. PLOT
# =========================================================

def plot_umap(embedding, labels):

    plt.figure(figsize=(7,6))

    plt.scatter(
        embedding[:,0],
        embedding[:,1],
        c=labels,
        s=10
    )

    plt.xlabel("UMAP1")
    plt.ylabel("UMAP2")

    plt.title("Allele-pattern clustering")

    plt.tight_layout()

    plt.savefig("allele_umap_clusters.png", dpi=300)

    plt.show()


# =========================================================
# 8. CLUSTER SUMMARY
# =========================================================

def summarize_clusters(df, labels):

    df["cluster"] = labels

    print("\nCluster sizes:")
    print(df["cluster"].value_counts())

    print("\nCluster fractions:")
    print(df["cluster"].value_counts(normalize=True))

    return df


# =========================================================
# RUN
# =========================================================

print("Loading phased SNPs...")
snps = load_phased_snps(VCF)

print(f"Loaded {len(snps)} phased SNPs")

print("Building ASE matrix...")
cell_snp = build_ase_matrix(BAM, snps)

print("Converting to dataframe...")
df = build_dataframe(cell_snp)

print("\nMatrix shape:")
print(df.shape)

df.to_csv("ase_matrix.tsv", sep="\t")

print("Running UMAP...")
embedding = run_embedding(df)

print("Clustering...")
labels = cluster_cells(embedding)

print("Summarizing clusters...")
df = summarize_clusters(df, labels)

df.to_csv("ase_clusters.tsv", sep="\t")

print("Plotting...")
plot_umap(embedding, labels)
# -----------------------------------
# Haplotype bias per cluster
# -----------------------------------

cluster_means = df.groupby("cluster").mean(numeric_only=True)

cluster_means.to_csv("cluster_haplotype_means.tsv", sep="\t")

print(cluster_means.iloc[:, :10])

hap_bias = df.drop(columns=["cluster"]).mean(axis=1)

df["hap_bias"] = hap_bias

print(df[["cluster", "hap_bias"]].head())
print("\nDONE")

feature_means = df.groupby("cluster").mean(numeric_only=True)

print(feature_means.iloc[:, :20])

feature_means.to_csv(
    "cluster_feature_means.tsv",
    sep="\t"
)

from sklearn.metrics import silhouette_score

score = silhouette_score(embedding, labels)

print("\nSilhouette score:", score)