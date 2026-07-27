#!/bin/bash
set -euo pipefail


################################################################################
# Single-cell RNA-seq X-inactivation analysis pipeline
#
# Steps:
# 1. Clean and sort BAM
# 2. Call variants
# 3. Extract heterozygous SNPs
# 4. Extract chrX SNPs
# 5. Cell-aware in-silico haplotype phasing
# 6. Count phased haplotypes
# 7. scRNA-seq allele-specific expression analysis
#
################################################################################


BAM=$1
SAMPLE=$2

THREADS=8


################################################################################
# Reference genome
################################################################################

REF=/home/lopde33/refdata-cellranger-GRCh38-3.0.0
GENOME=$REF/fasta/genome.fa


OUTDIR=results/${SAMPLE}

mkdir -p $OUTDIR


echo "================================="
echo "INPUT BAM: $BAM"
echo "SAMPLE: $SAMPLE"
echo "OUTPUT: $OUTDIR"
echo "================================="



################################################################################
# 1. CLEAN BAM
################################################################################

echo "Step 1: BAM filtering"


samtools view \
-F 4 \
-b \
$BAM \
> $OUTDIR/${SAMPLE}.clean.bam


samtools sort \
-@ $THREADS \
$OUTDIR/${SAMPLE}.clean.bam \
-o $OUTDIR/${SAMPLE}.clean.sorted.bam


samtools index \
$OUTDIR/${SAMPLE}.clean.sorted.bam



################################################################################
# Check 10x cell barcode and UMI tags
################################################################################

echo "Checking 10x CB and UB tags"


samtools view \
$OUTDIR/${SAMPLE}.clean.sorted.bam \
| head -1 \
| grep -E "CB:Z|UB:Z" || {

echo "ERROR: CB/UB tags not found in BAM"
echo "10x barcode information is required"
exit 1

}



################################################################################
# 2. VARIANT CALLING
################################################################################

echo "Step 2: variant calling"


bcftools mpileup \
--threads $THREADS \
--no-BAQ \
-d 10000 \
-Ou \
-f $GENOME \
$OUTDIR/${SAMPLE}.clean.sorted.bam \
|
bcftools call \
-mv \
-Oz \
-o $OUTDIR/${SAMPLE}.raw.vcf.gz


bcftools index \
$OUTDIR/${SAMPLE}.raw.vcf.gz



################################################################################
# 3. FILTER VARIANTS
################################################################################

echo "Step 3: filtering variants"


bcftools filter \
-i 'QUAL>10 && DP>3' \
$OUTDIR/${SAMPLE}.raw.vcf.gz \
-Oz \
-o $OUTDIR/${SAMPLE}.filtered.vcf.gz


bcftools index \
$OUTDIR/${SAMPLE}.filtered.vcf.gz



################################################################################
# 4. HETEROZYGOUS SNPs ONLY
################################################################################

echo "Step 4: extracting heterozygous SNPs"


bcftools view \
-g het \
$OUTDIR/${SAMPLE}.filtered.vcf.gz \
-Oz \
-o $OUTDIR/${SAMPLE}.het.vcf.gz


bcftools index \
$OUTDIR/${SAMPLE}.het.vcf.gz



################################################################################
# 5. CHROMOSOME X SNPs
################################################################################

echo "Step 5: extracting chrX heterozygous SNPs"


bcftools view \
-r X \
$OUTDIR/${SAMPLE}.het.vcf.gz \
-Oz \
-o $OUTDIR/chrX.het.vcf.gz


bcftools index \
$OUTDIR/chrX.het.vcf.gz



################################################################################
# 6. CELL-AWARE IN-SILICO PHASING
################################################################################

echo "Step 6: in-silico haplotype phasing"


python insilico_phase_v2.py \
--bam $OUTDIR/${SAMPLE}.clean.sorted.bam \
--vcf $OUTDIR/chrX.het.vcf.gz \
--out $OUTDIR/chrX.insilico.phased.vcf.gz \
--sample $SAMPLE \
--save-counts $OUTDIR/${SAMPLE}_insilico_counts \
--min-umi-per-call 1 \
--min-call-ratio 0.7 \
--min-joint-cells 3


echo "Phasing completed"



################################################################################
# 7. HAPLOTYPE COUNTS
################################################################################

echo "Step 7: counting phased haplotypes"


HAP0=$(bcftools query \
-f '[%GT\n]' \
$OUTDIR/chrX.insilico.phased.vcf.gz \
|
awk -F"|" '$1=="0"{c++} END{print c+0}')


HAP1=$(bcftools query \
-f '[%GT\n]' \
$OUTDIR/chrX.insilico.phased.vcf.gz \
|
awk -F"|" '$1=="1"{c++} END{print c+0}')


echo "Haplotype 0 SNP count: $HAP0"
echo "Haplotype 1 SNP count: $HAP1"



################################################################################
# 8. PLOT HAPLOTYPE DISTRIBUTION
################################################################################

echo "Step 8: creating haplotype plot"


cat > $OUTDIR/haplo_plot.py <<EOF

import matplotlib.pyplot as plt


labels = [
    "Haplo 0",
    "Haplo 1"
]

values = [
    $HAP0,
    $HAP1
]


plt.figure(figsize=(5,4))

plt.bar(
    labels,
    values
)

plt.ylabel(
    "Number of phased SNPs"
)

plt.title(
    "$SAMPLE ChrX haplotype"
)


for i,v in enumerate(values):
    plt.text(
        i,
        v,
        str(v),
        ha="center",
        va="bottom"
    )


plt.tight_layout()

plt.savefig(
"$OUTDIR/${SAMPLE}_Hap0_vs_Hap1.png",
dpi=300
)

EOF


python $OUTDIR/haplo_plot.py



################################################################################
# 9. scRNA-seq ALLELE-SPECIFIC EXPRESSION
################################################################################

echo "Step 9: scRNA-seq allele-specific expression"


mkdir -p $OUTDIR/scRNA_ASE


python umi.py \
--bam $OUTDIR/${SAMPLE}.clean.sorted.bam \
--vcf $OUTDIR/chrX.insilico.phased.vcf.gz \
--out $OUTDIR/scRNA_ASE \
--sample $SAMPLE


echo "ASE completed"



################################################################################
# DONE
################################################################################

echo "================================="
echo "DONE"
echo "Output directory:"
echo $OUTDIR
echo
echo "Phased VCF:"
echo $OUTDIR/chrX.insilico.phased.vcf.gz
echo
echo "ASE output:"
echo $OUTDIR/scRNA_ASE
echo
echo "Plot:"
echo $OUTDIR/${SAMPLE}_Hap0_vs_Hap1.png
echo "================================="
