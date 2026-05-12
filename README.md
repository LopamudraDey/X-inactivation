# X-inactivation

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

