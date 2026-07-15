#!/bin/bash
set -e

BAM=$1
SAMPLE=$2

THREADS=8

# Reference genome
# Change this path for Linux
REF=/home/lopde33/refdata-cellranger-GRCh38-3.0.0
GENOME=$REF/fasta/genome.fa

OUTDIR=results/${SAMPLE}
mkdir -p $OUTDIR

echo "================================="
echo "INPUT BAM: $BAM"
echo "SAMPLE: $SAMPLE"
echo "OUTPUT: $OUTDIR"
echo "================================="


# -------------------------
# 1. CLEAN BAM
# -------------------------
echo "Step 1: BAM filtering"

samtools view -b \
-F 4 \
$BAM \
> $OUTDIR/${SAMPLE}.clean.bam


samtools sort \
$OUTDIR/${SAMPLE}.clean.bam \
-o $OUTDIR/${SAMPLE}.clean.sorted.bam


samtools index \
$OUTDIR/${SAMPLE}.clean.sorted.bam



# -------------------------
# 2. VARIANT CALLING
# -------------------------
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


bcftools index \
$OUTDIR/${SAMPLE}.raw.vcf.gz



# -------------------------
# 3. FILTER VARIANTS
# -------------------------
echo "Step 3: filtering"


bcftools filter \
 -i 'QUAL>10 && DP>3' \
 $OUTDIR/${SAMPLE}.raw.vcf.gz \
 -Oz \
 -o $OUTDIR/${SAMPLE}.filtered.vcf.gz


bcftools index \
$OUTDIR/${SAMPLE}.filtered.vcf.gz



# -------------------------
# 4. HET SNPs ONLY
# -------------------------
echo "Step 4: heterozygous SNPs"


bcftools view \
-g het \
$OUTDIR/${SAMPLE}.filtered.vcf.gz \
-Oz \
-o $OUTDIR/${SAMPLE}.het.vcf.gz


bcftools index \
$OUTDIR/${SAMPLE}.het.vcf.gz



# -------------------------
# 5. CHR X EXTRACTION
# -------------------------
echo "Step 5: chrX extraction"


bcftools view \
-r X \
$OUTDIR/${SAMPLE}.het.vcf.gz \
-Oz \
-o $OUTDIR/chrX.het.vcf.gz


bcftools index \
$OUTDIR/chrX.het.vcf.gz



# -------------------------
# 6. PHASING
# -------------------------
echo "Step 6: WhatsHap phasing"


whatshap phase \
 --reference $GENOME \
 --ignore-read-groups \
 $OUTDIR/chrX.het.vcf.gz \
 $OUTDIR/${SAMPLE}.clean.sorted.bam \
 -o $OUTDIR/chrX.phased.vcf.gz


bcftools index \
$OUTDIR/chrX.phased.vcf.gz



# -------------------------
# 7. HAPLOTYPE COUNTS
# -------------------------
echo "Step 7: Hap0 vs Hap1 counting"


HAP0=$(bcftools query -f '[%GT\n]' \
$OUTDIR/chrX.phased.vcf.gz \
| awk -F"|" '$1=="0"{c++} END{print c+0}')


HAP1=$(bcftools query -f '[%GT\n]' \
$OUTDIR/chrX.phased.vcf.gz \
| awk -F"|" '$1=="1"{c++} END{print c+0}')


echo "Haplotype 0 SNP count: $HAP0"
echo "Haplotype 1 SNP count: $HAP1"



# -------------------------
# 8. CREATE BAR PLOT
# -------------------------
echo "Step 8: Creating Hap0/Hap1 plot"


cat > $OUTDIR/haplo_plot.py <<EOF

import matplotlib.pyplot as plt

labels = ["Haplo 0","Haplo 1"]
values = [$HAP0,$HAP1]

plt.figure(figsize=(5,4))

plt.bar(labels, values)

plt.ylabel("Number of phased SNPs")
plt.title("$SAMPLE ChrX haplotype")

for i,v in enumerate(values):
    plt.text(i,v,str(v),
             ha="center",
             va="bottom")

plt.tight_layout()

plt.savefig(
"$OUTDIR/${SAMPLE}_Hap0_vs_Hap1.png",
dpi=300
)

EOF


python $OUTDIR/haplo_plot.py

# -------------------------
# 9. scRNA-seq ALLELE-SPECIFIC EXPRESSION
# -------------------------

echo "Step 9: scRNA-seq allele-specific UMI analysis"


mkdir -p $OUTDIR/scRNA_ASE


python allele_specific_umi.py \
--bam $BAM \
--vcf $OUTDIR/chrX.phased.vcf.gz \
--out $OUTDIR/scRNA_ASE \
--sample $SAMPLE


echo "scRNA-seq ASE completed"


echo "================================="
echo "================================="
echo "DONE"
echo "Output:"
echo "$OUTDIR"
echo "Plot:"
echo "$OUTDIR/${SAMPLE}_Hap0_vs_Hap1.png"
echo "================================="