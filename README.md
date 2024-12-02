# Identifying cis-eQTLs in 1000 Genomes and Generating Polygenic Risk Scores for Heart Disease

# Overview
## Part 1:
Identify ciseQTLs in 1000 Genomes LCLs, where cis is defined by +/- 500 Mb of the gene body.
a. across all genes genome-wide (n = 20K protein coding genes) Note: might need to subset gene expression data. If so, use gene_annot.txt.gz to figure out which genes are protein coding.

## Part 2: 
Predict genetic risk for different diseases in 1000 Genomes individuals. Using Tiffany’s 23andme data, where does she fall on the distribution? Are the results reasonable?

# Retrieving the Data Locally
(1) Download 1000 Genomes genotype data downloaded from the [LDREF](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2) \
(2) Download the [gene expression](https://drive.google.com/file/d/1EIj8pALKDAsck7JWAk31JZWPQtF0aS9w/view?usp=sharing) and [annotation files](https://drive.google.com/file/d/1db5vdekQsJ0f7aP34hnr6YXVuV-ov-2M/view?usp=sharing) \
(3) Download [Tiffany's 23andme Data](https://drive.google.com/file/d/1m3eVt-FgkfmBW2w2UW_TDtzlXVpiBaCC/view?usp=sharing) \
(4) Download the [GRCh37](https://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz) reference genome
(4) Download [Diabetes GWAS](https://drive.google.com/file/d/1tCAYBW8KPqNlx3AGlbEWcOCZiOsXEQKe/view?usp=sharing) \
(5) Download [Heart Disease GWAS](https://drive.google.com/file/d/1muZDyipGbDmzNhRnNcrSeDL_1faY9l4c/view?usp=sharing) \
(6) Download [Parkinson's GWAS](https://drive.google.com/file/d/1623dNwgbtf4cv0fcQuKKaIOEGErHcRs5/view?usp=sharing) \

# Running the Project
- To install the dependencies, run the following command from the root directory of the project: `pip install -r requirements.txt`
- Download the version of [Plink](https://www.cog-genomics.org/plink/1.9/) that works with your device.
- Download the version of [Plink2](https://www.cog-genomics.org/plink/2.0/) that works with your device.

# Completing Part 1

# Completing Part 2



