#!/bin/bash
set -e

BAM=$1
SAMPLE=$2

THREADS=8

REF=/mnt/c/Users/lopde33/Project/refdata-gex-GRCh38-2020-A
GENOME=$REF/fasta/genome.fa

OUTDIR=results/${SAMPLE}
mkdir -p $OUTDIR

echo "=== INPUT BAM: $BAM ==="
echo "=== OUTPUT DIR: $OUTDIR ==="


# =========================================================
# 1. CLEAN BAM (improved filtering)
# =========================================================
echo "Step 1: BAM filtering"

samtools view -b \
  -F 2308 \
  -q 10 \
  $BAM \
  > $OUTDIR/${SAMPLE}.clean.bam

samtools sort \
  $OUTDIR/${SAMPLE}.clean.bam \
  -o $OUTDIR/${SAMPLE}.clean.sorted.bam

samtools index $OUTDIR/${SAMPLE}.clean.sorted.bam


# =========================================================
# 2. VARIANT CALLING
# =========================================================
echo "Step 2: variant calling"

bcftools mpileup \
  --threads $THREADS \
  --no-BAQ \
  -d 10000 \
  -Ou \
  -f $GENOME \
  $OUTDIR/${SAMPLE}.clean.sorted.bam \
| bcftools call \
  -mv \
  -Oz \
  -o $OUTDIR/${SAMPLE}.raw.vcf.gz

bcftools index $OUTDIR/${SAMPLE}.raw.vcf.gz


# =========================================================
# 3. FILTER VARIANTS
# =========================================================
echo "Step 3: filtering"

bcftools filter \
  -i 'QUAL>10 && DP>3' \
  $OUTDIR/${SAMPLE}.raw.vcf.gz \
  -Oz \
  -o $OUTDIR/${SAMPLE}.filtered.vcf.gz

bcftools index $OUTDIR/${SAMPLE}.filtered.vcf.gz

echo "Filtered VCF variant count:"
bcftools view -H $OUTDIR/${SAMPLE}.filtered.vcf.gz | wc -l


# =========================================================
# 4. HET SNPs ONLY
# =========================================================
echo "Step 4: heterozygous SNPs"

bcftools view \
  -g het \
  $OUTDIR/${SAMPLE}.filtered.vcf.gz \
  -Oz \
  -o $OUTDIR/${SAMPLE}.het.vcf.gz

bcftools index $OUTDIR/${SAMPLE}.het.vcf.gz


# =========================================================
# 5. CHR X EXTRACTION (IMPORTANT FIX: X not chrX)
# =========================================================
echo "Step 5: chrX extraction"

bcftools view \
  -r chrX \
  $OUTDIR/${SAMPLE}.het.vcf.gz \
  -Oz \
  -o $OUTDIR/chrX.het.vcf.gz

bcftools index $OUTDIR/chrX.het.vcf.gz

echo "chrX SNP count:"
bcftools view -H $OUTDIR/chrX.het.vcf.gz | wc -l


# =========================================================
# 6. PHASING
# =========================================================
echo "Step 6: phasing"

whatshap phase \
  --reference $GENOME \
  --ignore-read-groups \
  $OUTDIR/chrX.het.vcf.gz \
  $OUTDIR/${SAMPLE}.clean.sorted.bam \
  -o $OUTDIR/chrX.phased.vcf.gz

bcftools index $OUTDIR/chrX.phased.vcf.gz


# =========================================================
# 7. STATS
# =========================================================
echo "Step 7: WhatsHap stats"

whatshap stats $OUTDIR/chrX.phased.vcf.gz


# =========================================================
# DONE
# =========================================================
echo "DONE"
echo "Output directory: $OUTDIR"
