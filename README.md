# X-inactivation
# Step 1: Download file from GEO
prefetch  SRR16922266

fastq-dump SRR16922266 --split-files --gzip \
This will give you 3 files: _1,_2 and _3 \
# Use zcat to check 
Example: SRR16922266_1.fastq.gz | head -20  \
SRR16922266_1.fastq.gz (Length 8, Sample Index or i7 index) \
SRR16922266_2.fastq.gz (Length 26, the Cell Barcode and UMI (Unique Molecular Identifier).\
SRR16922266_3.fastq.gz ( Length 91, the biological read)

# Step 2: Creating the Bam file (Either STAR or Cellranger)
STAR \
--runThreadN 8 \
--genomeDir /local/data1/lopde33/STAR-2.7.11b/source/refdata-gex-GRCh38-2024-A/star \
--readFilesIn \
/media/lopde33/T7/SRR16/SRR16922266_2.fastq.gz \
/media/lopde33/T7/SRR16/SRR16922266_3.fastq.gz \
--readFilesCommand zcat \
--soloType CB_UMI_Simple \
--soloFeatures GeneFull \
--soloCBstart 1 --soloCBlen 16 \
--soloUMIstart 17 --soloUMIlen 10 \
--soloCBwhitelist None \
--outFilterMultimapNmax 1 \
--outFilterMismatchNmax 2 \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI AS nM CB UB GX GN \
--outFileNamePrefix ./SRR16922266_STARsolo/

# Create Conda environment

conda create -n haplo python=3.10 -y 

conda activate haplo  

conda install -c conda-forge -c bioconda -y \
samtools \
bcftools \
htslib \
pysam \
numpy \
pandas \
matplotlib \
tqdm \
whatshap

