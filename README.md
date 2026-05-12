# X-inactivation
# Step 1: Download file from GEO
prefetch  SRR16922266

fastq-dump SRR16922266 --split-files --gzip \
This will give you 3 files: _1,_2 and _3 \
SRR16922266_1.fastq.gz \
SRR16922266_2.fastq.gz\
SRR16922266_3.fastq.gz 



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

