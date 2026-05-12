#!/bin/bash
set -e

BAM=$1
SAMPLE=$2

THREADS=8

REF=/mnt/c/Users/lopde33/Project/refdata-gex-GRCh38-2024-A
GENOME=$REF/fasta/genome.fa

OUTDIR=results/${SAMPLE}
mkdir -p $OUTDIR

echo "=== INPUT BAM: $BAM ==="


# -------------------------
# 1. CLEAN BAM
# -------------------------
echo "Step 1: BAM filtering"
samtools view -b -F 4 $BAM \
> $OUTDIR/${SAMPLE}.clean.bam

samtools sort $OUTDIR/${SAMPLE}.clean.bam \
-o $OUTDIR/${SAMPLE}.clean.sorted.bam

samtools index $OUTDIR/${SAMPLE}.clean.sorted.bam


# -------------------------
# 2. VARIANT CALLING
# -------------------------
echo "Step 2: variant calling"

bcftools mpileup \
 -Ou \
 -f $GENOME \
 $OUTDIR/${SAMPLE}.clean.sorted.bam \
| bcftools call -mv -Oz -o $OUTDIR/${SAMPLE}.raw.vcf.gz

bcftools index $OUTDIR/${SAMPLE}.raw.vcf.gz


# -------------------------
# 3. FILTER VARIANTS
# -------------------------
echo "Step 3: filtering"

bcftools filter \
 -i 'QUAL>20 && DP>5' \
 $OUTDIR/${SAMPLE}.raw.vcf.gz \
 -Oz -o $OUTDIR/${SAMPLE}.filtered.vcf.gz

bcftools index $OUTDIR/${SAMPLE}.filtered.vcf.gz


# -------------------------
# 4. HET SNPs ONLY
# -------------------------
bcftools view -g het \
$OUTDIR/${SAMPLE}.filtered.vcf.gz \
-Oz -o $OUTDIR/${SAMPLE}.het.vcf.gz

bcftools index $OUTDIR/${SAMPLE}.het.vcf.gz


# -------------------------
# 5. CHR X EXTRACTION
# -------------------------
bcftools view -r chrX \
$OUTDIR/${SAMPLE}.het.vcf.gz \
-Oz -o $OUTDIR/chrX.het.vcf.gz

bcftools index $OUTDIR/chrX.het.vcf.gz


# -------------------------
# 6. PHASING
# -------------------------
echo "Step 4: phasing"

whatshap phase \
 --reference $GENOME \
 --ignore-read-groups \
 $OUTDIR/chrX.het.vcf.gz \
 $OUTDIR/${SAMPLE}.clean.sorted.bam \
 -o $OUTDIR/chrX.phased.vcf.gz

bcftools index $OUTDIR/chrX.phased.vcf.gz


echo "DONE"
echo "Output: $OUTDIR"