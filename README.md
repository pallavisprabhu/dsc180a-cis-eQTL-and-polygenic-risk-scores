# Identification of cis-eQTLs in 1000 Genomes LCLs and Predicting Genetic Risk for Multiple Diseases with Polygenic Risk Scores

# Overview
## Part 1:
Identify ciseQTLs in 1000 Genomes LCLs, where cis is defined by +/- 500 Kb of the gene body across all protein-coding genes genome-wide

## Part 2: 
Predict genetic risk for different diseases in 1000 Genomes individuals. Using Tiffany’s 23andme data, where does she fall on the distribution?

# Retrieving the Data Locally
(1) Download 1000 Genomes genotype data downloaded from the [LDREF](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2) \
(2) Download the [gene expression](https://drive.google.com/file/d/1EIj8pALKDAsck7JWAk31JZWPQtF0aS9w/view?usp=sharing) and [annotation files](https://drive.google.com/file/d/1db5vdekQsJ0f7aP34hnr6YXVuV-ov-2M/view?usp=sharing) \
(3) Download [Tiffany's 23andme Data](https://drive.google.com/file/d/1m3eVt-FgkfmBW2w2UW_TDtzlXVpiBaCC/view?usp=sharing) \
(4) Download the [GRCh37](https://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz) reference genome \
(5) Download [Diabetes GWAS](https://drive.google.com/file/d/1tCAYBW8KPqNlx3AGlbEWcOCZiOsXEQKe/view?usp=sharing) \
(6) Download [Heart Disease GWAS](https://drive.google.com/file/d/1muZDyipGbDmzNhRnNcrSeDL_1faY9l4c/view?usp=sharing) \
(7) Download [Parkinson's GWAS](https://drive.google.com/file/d/1623dNwgbtf4cv0fcQuKKaIOEGErHcRs5/view?usp=sharing) 

# Running the Project
- To install the dependencies, run the following command from the root directory of the project: `pip install -r requirements.txt`
- Download `bcftools` either [here](https://www.htslib.org/download/) or of you have Homebrew installed, by: `brew install bcftools`
- Download the version of [Plink](https://www.cog-genomics.org/plink/1.9/) that works with your device.
- Download the version of [Plink2](https://www.cog-genomics.org/plink/2.0/) that works with your device.

# Completing Part 1
Create a Jupyter Notebook to host your work and import the packages as noted by `requirements.txt` in addition to the `subprocess` and `os` library:
```py
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
```
## Define the function `ciseQTL`
`ciseQTL` takes in three arguments, the chromosome number as an integer, the file path to the gene expression data as a string, and the file path for the annotation data as a string, and outputs a summary statistics file with the Chromosome, Gene, SNP, Beta_0, Beta_1, $R^2$, Standard Error, and P-value for every SNP for every gene on that chromosome.

### Define the function `extract_protein`
The gene expression data contains expression for various types of genes in the genome. Create the helper function `extract_protein` which returns the gene expression data and annotation data for only the protein-coding genes:

Read in the expression data and gene expression data:
```py
expression_df = pd.read_csv('GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', delim_whitespace=True)
annotations = pd.read_csv('gene_annot.txt', delim_whitespace=True)
```

In order to find cis-eQTLs for the protein-coding genes in the genome, subset the gene expression data to only include protein-coding genes using the annotation DataFrame:
```py
protein_coding = annotations[annotations['TYPE'] == 'protein_coding']
intersection = np.intersect1d(np.array(expression_df['TargetID'].apply(lambda x: x.split('.')[0])), np.array(protein_coding['SYM']))
protein_exp = expression_df[expression_df['TargetID'].apply(lambda x: x.split('.')[0]).apply(lambda x: x in intersection)]
```
The function should return the `protein_exp` DataFrame which contains all the expression levels for the protein-coding genes and the `protein_coding` DataFrame which contains the annotations for all the protein-coding genes:
```py
def extract_protein(exp_path, annot_path):
    expression_df = pd.read_csv(exp_path, delim_whitespace=True)
    annotations = pd.read_csv(annot_path, delim_whitespace=True)
    protein_coding = annotations[annotations['TYPE'] == 'protein_coding']
    intersection = np.intersect1d(np.array(expression_df['TargetID'].apply(lambda x: x.split('.')[0])), np.array(protein_coding['SYM']))
    protein_exp = expression_df[expression_df['TargetID'].apply(lambda x: x.split('.')[0]).apply(lambda x: x in intersection)]
    return protein_exp, protein_coding
```

Call the previously defined `extract_protein` function and returns the gene expression and annotation DataFrames with nly the protein-coding genes. The, query the protein_coding DataFrame to only include genes on that input chromosome. 
```py
protein_exp, protein_coding = extract_protein(exp_path, annot_path)
chrom = protein_coding[protein_coding['CHR'] == chromosome]
```

Then, it iterates through all the rows, one for each gene, to create a coordinate row that defines where the gene stops and ends. Since cis is defined as 500Kb (kilobases) of the gene body, we create that range by making the start of the gene body 500Kb upstream the gene's start and the end of the gene body 500Kb downstream the gene's end. Each of this coordinates for each gene is saved all are loaded into a singular DataFrame:
```py
chrom = protein_coding[protein_coding['CHR'] == chromosome]
chr_str = 'chr' + str(chromosome)
gene_snps = []
for i, gene in chrom.iterrows():
  start = gene['START']
  end = gene['STOP']
  range_start = max(0, (start-500000))
  range_end = end + 500000
  gene_snps.append({'chromosome': chr_str, 'start': range_start, 'end': range_end, 'gene': gene['ID'], " ": 0, '': '+' })
  chr = pd.DataFrame(gene_snps)
```

Iterating through each of these coordinate rows of the `chr` DataFrame, write a tab delmited `txt` file for that row. Note that each gene should have its own `txt` file. Using this `txt` file, Plink2 can filter the SNPs within that gene body (the start and stop) using the bfiles from the 1000 Genomes. Plink2 may not be able to complete this for some genes due to being near the start or end of the chromosome, so keep track of the successes and errors to know which genes make it to the final product:
```py
errors = []
successes = []
for _, row in chr.iterrows():
  chrom = row['chromosome']
  start = row['start']
  end = row['end']
  gene = row['gene']

  bfile_prefix = f"./LDREF/1000G.EUR.{chromosome}" 
  folder = f"./coords_{chromosome}"
  file_path = os.makedirs(folder, exist_ok=True)
  coord_filename = f"./coords_{chromosome}/{gene}_coord.txt"
        
  with open(coord_filename, "w") as coord_file:
  coord_file.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t+")
            
  folder = f"./LDREF_{chromosome}"
  file_path = os.makedirs(folder, exist_ok=True)
  output_prefix = f"./LDREF_{chromosome}/1000G.EUR.{chromosome}.{gene}"
  command = [
    "./plink2", 
    "--bfile", f"./LDREF/1000G.EUR.{chromosome}", 
    "--extract", "bed1", coord_filename,
    "--out", output_prefix,
    "--make-bed"]

  try:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    successes.append(gene)
        except subprocess.CalledProcessError as e:
            errors.append(gene)
```

Now each gene has a set of bfiles (`bim/bed/fam`) with the SNPs in the cis-region. Now filter the gene expression DataFrame to include only genes on the chromosome. In order to read the genotype data for each gene now, convert the file into `.raw` format using Plink2:
```py
exp_protein = protein_exp[protein_exp['Chr'] == str(chromosome)]

folder = f"./genotypes_raw_{chromosome}"
file_path = os.makedirs(folder, exist_ok=True)

for gene in chr['gene'].values:
  command = [
  "./plink2",
  "--bfile",
  f"./LDREF_{chromosome}/1000G.EUR.{chromosome}.{gene}",
  "--export", 
  "A",
  "--out",
  f"./genotypes_raw_{chromosome}/genotype_data_{gene}"]
        
  try:
    result = subprocess.run(command, check=True, capture_output=True, text=True)
  except subprocess.CalledProcessError as e:
    pass
```

Finally, compute a linear regression for each SNP on each gene where the features are the genotypes of all individuals in the sample for that SNP and the target vector is gene expression. In order to do this, first read in the genotype DataFrame for that gene as computed above and find the intersection of samples between the expression DataFrame and genotype DataFrame as not all genes may have expression values in this data and vice versa. Note that this is only done on the `successes`, the genes for which bfiles were able to be created:
```py
results = []
for gene in successes:
  gene_id = protein_coding[protein_coding['ID'] == gene]['SYM']
  if gene not in successes:
    continue
  exp = exp_protein[exp_protein['TargetID'].str.split('.').apply(lambda x: x[0]) == gene_id.values[0]]
  if exp.shape[0] == 0:
    continue
  exp_df = exp.iloc[:,4:].transpose()
    
  gen_df = pd.read_csv(f"./genotypes_raw_{chromosome}/genotype_data_{gene}.raw", sep='\t')
  gen_df.set_index('IID', inplace=True)

  common_samples = exp_df.index.intersection(gen_df.index)
  exp_df = exp_df.loc[common_samples]
  gen_df = gen_df.loc[common_samples]

  y = exp_df.iloc[:, 0].values
  X = gen_df.loc[common_samples].iloc[:, 5:]
  for snp in X.columns:
  try:
    slope, intercept, r_value, p_value, std_error = stats.linregress(X[[snp]].squeeze(), y)
    results.append({'Chr': chromosome, 'Gene': gene, 'SNP': snp,'beta_0': intercept, 'beta_1': slope, 'R_sq': (r_value**2), 'Standard Error': std_error,'P-value': p_value}) 
  except:
    continue
```

Save the linear regression results for each SNP for each gene on the chromosome and combine this information into a DataFrame. Save this as a file:
```py
results_df = pd.DataFrame(results)

folder = f"./ciseqtls_sumstats"
file_path = os.makedirs(folder, exist_ok=True)

results_df.to_csv(f'./ciseqtls_sumstats/sumstats_{chromosome}', sep='\t', index=False)
```

The final `cis-eQTL` function should look like this:
```py
def ciseQTL(chromosome): #chromosome as an integer
    chrom = protein_coding[protein_coding['CHR'] == chromosome]
    chr_str = 'chr' + str(chromosome)
    gene_snps = []
    for i, gene in chrom.iterrows():
        start = gene['START']
        end = gene['STOP']
        range_start = max(0, (start-500000))
        range_end = end + 500000
        gene_snps.append({'chromosome': chr_str, 'start': range_start, 'end': range_end, 'gene': gene['ID'], " ": 0, '': '+' })
    chr = pd.DataFrame(gene_snps)
    errors = []
    successes = []
    for _, row in chr.iterrows():
        chrom = row['chromosome']
        start = row['start']
        end = row['end']
        gene = row['gene']
        
        bfile_prefix = f"./LDREF/1000G.EUR.{chromosome}" 
        folder = f"./coords_{chromosome}"
        file_path = os.makedirs(folder, exist_ok=True)
        coord_filename = f"./coords_{chromosome}/{gene}_coord.txt"
        
        with open(coord_filename, "w") as coord_file:
            coord_file.write(f"{chrom}\t{start}\t{end}\t{gene}\t0\t+")
            
        folder = f"./LDREF_{chromosome}"
        file_path = os.makedirs(folder, exist_ok=True)
        output_prefix = f"./LDREF_{chromosome}/1000G.EUR.{chromosome}.{gene}"
        command = [
            "./plink2", 
            "--bfile", f"./LDREF/1000G.EUR.{chromosome}", 
            "--extract", "bed1", coord_filename,
            "--out", output_prefix,
            "--make-bed"]
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            successes.append(gene)
        except subprocess.CalledProcessError as e:
            errors.append(gene)

    exp_protein = protein_exp[protein_exp['Chr'] == str(chromosome)]

    folder = f"./genotypes_raw_{chromosome}"
    file_path = os.makedirs(folder, exist_ok=True)

    for gene in chr['gene'].values:
        command = [
        "./plink2",
        "--bfile",
        f"./LDREF_{chromosome}/1000G.EUR.{chromosome}.{gene}",
        "--export", 
        "A",
        "--out",
        f"./genotypes_raw_{chromosome}/genotype_data_{gene}"]
        
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            pass

    results = []
    for gene in successes:
        gene_id = protein_coding[protein_coding['ID'] == gene]['SYM']
        if gene not in successes:
            continue
        exp = exp_protein[exp_protein['TargetID'].str.split('.').apply(lambda x: x[0]) == gene_id.values[0]]
        if exp.shape[0] == 0:
            continue
        exp_df = exp.iloc[:,4:].transpose()
    
        gen_df = pd.read_csv(f"./genotypes_raw_{chromosome}/genotype_data_{gene}.raw", sep='\t')
        gen_df.set_index('IID', inplace=True)

        common_samples = exp_df.index.intersection(gen_df.index)
        exp_df = exp_df.loc[common_samples]
        gen_df = gen_df.loc[common_samples]

        y = exp_df.iloc[:, 0].values
        X = gen_df.loc[common_samples].iloc[:, 5:]
        for snp in X.columns:
            try:
                slope, intercept, r_value, p_value, std_error = stats.linregress(X[[snp]].squeeze(), y)
                results.append({'Chr': chromosome, 'Gene': gene, 'SNP': snp,'beta_0': intercept, 'beta_1': slope, 'R_sq': (r_value**2), 'Standard Error': std_error,'P-value': p_value}) 
            except: 
                continue
            
    results_df = pd.DataFrame(results)

    folder = f"./ciseqtls_sumstats"
    file_path = os.makedirs(folder, exist_ok=True)

    results_df.to_csv(f'./ciseqtls_sumstats/sumstats_{chromosome}', sep='\t', index=False)
    
    return results_df
```
## Defining the plotting function
In order to plot the cis-eQTL, you need to plot it per gene since plotting the cis-eQTLs for the chromosome in its entirity would not make sense as we are interested in seeing which SNPs are statitsically significant in impacting gene expression for a specific gene.
Create the function `cis_plot` which has two parameters, the file path to the summary statistics file as created from the `cis-eQTL` function above, and a gene on that chromosome. The function filters the summary statistics to include information for only that gene and plots a locus plot, marking the threshold for signficance with a red dotted line for easier visualization. Note that instead of plotting the raw P-values, the $log_{10}(P-value)$ is plotted for better scaling. This means that lower P-values have higher $log_{10}(P-value)$ and data points above the line are significant. Save the figure as a `png`. 

```py
def cis_plot(sum_stats_path, gene):
  sum_stats = pd.read_csv(sum_stats_path, sep='\t')
  chr = sum_stats['Chr'].iloc[0]
  filtered = sum_stats[sum_stats['Gene'] == gene].reset_index()
  sns.set(font_scale=1)
  plt.figure(figsize=(10, 6))
  plt.title(f'cis-eQTLs of {gene} on Chromosome {chr}')
  plt.scatter(filtered.index, -np.log10(filtered['P-value']))
  plt.axhline(y=-np.log10(0.05), color='r', linestyle='-')
  plt.xlabel('SNPs')
  plt.ylabel('-log10(P-value)')

  folder = f"./ciseqtls_graph"
  file_path = os.makedirs(folder, exist_ok=True)
  plt.savefig(f"./ciseqtls_graph/{gene}_chromosome{chr}.png", dpi=300, bbox_inches='tight')
    
  plt.show()
```

# Part 2

## Convert Tiffany's data into Plink file format
First, convert Tiffany;'s data into vcf format. For this you will need the GRCh37 reference genome. Note that the following should be done in terminal.
```bash
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa ref.fa #rename the file

bgzip -c tiffany.txt > tiffany.tab.gz #zip Tiffany's data

bcftools convert --tsv2vcf tiffany.tab.gz -f ref.fa -s Tiffany -Ob -o tiffany.bcf #convert the data into bcf file format

bcftools view tiffany.bcf -Ov -o tiffany.vcf #convert to vcf file format

./plink2 --vcf tiffany.vcf --make-bed --out tiffany_plink #convert to plink files

```

## Merge chromosome data
Since all the chromosome data from 1000 Genomes are in separate files and Tiffany has her own file, for ease of the program, merge the data into one set of Plink files. Create a text file `all_chrs.txt` with the paths to all the 1000 Genomes data and Tiffany's data, exluding the first chromosome as that will be the base file:
```text
./LDREF/1000G.EUR.2
./LDREF/1000G.EUR.3
./LDREF/1000G.EUR.4
...
./LDREF/1000G.EUR.21
./LDREF/1000G.EUR.22
tiffany_plink
```

In terminal, run the following command in Plink to merge the data:
```bash
./plink --bfile ./LDREF/1000G.EUR.1 --merge-list all_chrs.txt --make-bed -out all_chromosomes
```

## Generating PRS scores
To generate PRS scores from GWAS data, define the function `prs` that has two parameters `gwas`, the file path to the `tsv` file, and `disease`, the name of the disease as a string. The first step is the clean the GWAS `tsv` file to make it usuable for to generate PRS scores. Note that the data wrangling done is specific to the GWAS data obtained from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/) and may not transfer to GWAS obtained from other sources and may have different column names.

### Preprocessing the GWAS data
Define the following helper function `extract_allele` to extract the effect allele from the SNP id (rsid). This will be used to clean the GWAS data to extarct the effect allele from the `'STRONGEST SNP-RISK ALLELE'` column of the GWAS.
```py
def extract_allele(rsid):
    al = rsid.split('-')[1]
    if al.isalpha():
        return al
    return '.'
```

The GWAS data from the GWAS catalog needs some preprocessing in order to be used for the PRS score:
- Load in the data in a DataFrame
- Subset the DataFrame to include only these columns: `'CHR_ID'`,`'SNPS'`, `'OR or BETA'`, `'STRONGEST SNP-RISK ALLELE'`,`'P-VALUE'`
- Apply the `extract_allele` function to the `'STRONGEST SNP-RISK ALLELE'` column and save this as the column `'A1'`
- Rename the column `'OR or BETA'` to `'BETA'`
- Drop the columns `'STRONGEST SNP-RISK ALLELE'` and `'OR or BETA'`
- Some of the SNPs have multiple rows. Keep the ones with the lowest P-values
- Create the column `'EFFECT'` that is the natural logarthim applied to the `'BETA'`
- Remove all rows where `'A1'` is a `.` and drop all null values
- Save this as a tab delimited file

```py
def clean_gwas(gwas, dis):
    disease = pd.read_csv(file, sep= '\t')
    disease = disease[['CHR_ID','SNPS', 'OR or BETA', 'STRONGEST SNP-RISK ALLELE','P-VALUE']]
    disease['A1'] = disease['STRONGEST SNP-RISK ALLELE'].apply(extract_allele)
    disease['BETA'] = disease['OR or BETA']
    disease = disease.drop(columns=['STRONGEST SNP-RISK ALLELE', 'OR or BETA' ])
    disease = disease.loc[disease.groupby('SNPS')['P-VALUE'].idxmin()].reset_index(drop=True)
    disease['EFFECT'] = disease['BETA'].apply(np.log)
    disease = disease[disease['A1'] != '.'].dropna()

    diseases.to_csv(f'{dis}.txt', sep='\t', index=False)

    return diseases
```
Call the helper function `clean_gwas` on the file path to the GWAS data which will save the cleaned GWAS as a `txt` file to so it can be used in Plink:
```py
_ =  clean_gwas(gwas_file, disease)
```
To start creating a PRS score, clump and threshold the data. Thresholding is the process of filtering out the SNPs that do not make the threshold of a specified P-value and $\text{R}^2$. Moreover, clumping is the process of removing dependent SNPs such that to limit linkage disequilibrium from closely associated or linked SNPs.
```py
command= [
    './plink2',
    '--bfile', f'all_chromosomes',
    '--clump-p1', '.0001',
    '--clump-r2', '0.2',
    '--clump-kb', '500',
    '--clump', f'{disease}.txt', 
    '--clump-snp-field', 'SNPS', 
    '--clump-field', 'P-VALUE',
    '--out', f'./clumped_{disease}']
    
try:
  result = subprocess.run(command, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
  pass
```

In order to extract the SNP Ids from the clumped output, run the following command:

```py
command = f"awk 'NR!=1{{print $3}}' ./clumped_{disease}.clumps > ./PRS.SNPs.{disease}"
subprocess.run(command, shell=True, check=True)
```

Now we create Plink bfiles (`bed`,`bim`,`fam`) including only those extracted SNPs:
```py
command = f"awk 'NR!=1{{print $3}}' ./clumped_{disease}.clumps > ./PRS.SNPs.{disease}"
subprocess.run(command, shell=True, check=True)

folder = f"./b_{disease}"
file_path = os.makedirs(folder, exist_ok=True)

command2 = [
    './plink2', 
    '--bfile', f'all_chromosomes',
    '--extract', f"./PRS.SNPs.{disease}", 
    '--make-bed',
    '--out', f'./b_{disease}/1000G_eur_PRS']

try:
  result = subprocess.run(command2, check=True, capture_output=True, text=True)
except subprocess.CalledProcessError as e:
  pass 
```

Now we need to create a score file. To do this, we get all the extracted SNPs from output from above. Then, we read in the GWAS data. Find the intersection of SNPs in both files. For the SNPs that are from the clumped output and are in the GWAS, we extract its `'SNP'`, `'A1' `, and `EFFECT` and save these values to the score file

```py
snp_file = f"./PRS.SNPs.{disease}"
y = pd.read_csv(snp_file, header=None, sep="\t")
y = y[y[0].str.len() > 0][0].tolist()
dis = pd.read_csv(f'{disease}.txt', sep='\t')
m = dis['SNPS'].isin(y)
score = pd.DataFrame({
    'SNP': dis['SNPS'][m],
    'A1': dis['A1'][m],
    'BETA': dis['EFFECT'][m]})

score.to_csv(f"./score_file_{disease}.txt", index=False, header=False, sep="\t", quoting=3)
```
Now, we can generate the PRS scores from the score file using Plink and return the PRS scores as a DataFrame:
```py
command3 = [
    './plink2', '--bfile', f'./b_{disease}/1000G_eur_PRS',
    '--out', f'./PRS_{disease}',
    '--score', f"./score_file_{disease}.txt", '1', '2', '3']

try:
    result = subprocess.run(command3, check=True)
except subprocess.CalledProcessError as e:
    pass 

prs_score = pd.read_csv(f'PRS_{disease}.sscore', sep='\t')
```
### Visualizing the PRS Scores
It's more informative and helpful to visulize the distribution of PRS scores on a plot like a histogram, especially since we want to isolate Tiffany's score and see where she falls in the distribution. Define the helper function `prs_plot` that takes in the file path to the prs data from above and the name of the disease as a string and outputs a histogram of the distrubution, labeling where Tiffany falls.
```py
def prs_plot(prs_df, dis):
    prs = pd.read_csv(f'./PRS_{dis}.sscore', sep='\t')
    tiff = prs[prs['IID'] == 'Tiffany']['SCORE1_AVG'].values[0]

    sns.set(font_scale=1)
    plt.figure(figsize=(10, 6))
    plt.hist(prs_df['SCORE1_AVG'], bins=50, alpha=0.7, label='1000 Genomes')
    plt.axvline(tiff_heart, color='r', linestyle='--', label='Tiffany')
    plt.xlabel('PRS Score')
    plt.ylabel('Frequency')
    plt.title(f"PRS Distribution for {dis} Disease")
    plt.legend()

    folder = f"./prs_graphs"
    file_path = os.makedirs(folder, exist_ok=True)
    plt.savefig(f"./prs_graphs/{dis}.png", dpi=300, bbox_inches='tight')

    plt.show()
```
Call this function on the DataFrame `prs_score` that contains the PRS scores for each individual in the sample.

The final `prs` function should like:
```py
def prs(gwas_file, disease):
    def extract_allele(rsid):
        al = rsid.split('-')[1]
        if al.isalpha():
            return al
        return '.'
        
    def clean_gwas(file, dis):
        disease = pd.read_csv(file, sep= '\t')
        disease = disease[['CHR_ID','SNPS', 'OR or BETA', 'STRONGEST SNP-RISK ALLELE','P-VALUE']]
        disease['A1'] = disease['STRONGEST SNP-RISK ALLELE'].apply(extract_allele)
        disease['BETA'] = disease['OR or BETA']
        disease = disease.drop(columns=['STRONGEST SNP-RISK ALLELE', 'OR or BETA' ])
        disease = disease.loc[disease.groupby('SNPS')['P-VALUE'].idxmin()].reset_index(drop=True)
        disease['EFFECT'] = disease['BETA'].apply(np.log)
        disease = disease[disease['A1'] != '.'].dropna()
        disease.to_csv(f'{dis}.txt', sep='\t', index=False)
        return disease
        
    _ = clean_gwas(gwas_file, disease)
    
    command= [
            './plink2',
            '--bfile', f'all_chromosomes',
            '--clump-p1', '.0001',
            '--clump-r2', '0.2',
            '--clump-kb', '500',
            '--clump', f'{disease}.txt', 
            '--clump-snp-field', 'SNPS', 
            '--clump-field', 'P-VALUE',
            '--out', f'./clumped_{disease}']
    
    try:
        result = subprocess.run(command, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        pass

    command = f"awk 'NR!=1{{print $3}}' ./clumped_{disease}.clumps > ./PRS.SNPs.{disease}"
    subprocess.run(command, shell=True, check=True)

    folder = f"./b_{disease}"
    file_path = os.makedirs(folder, exist_ok=True)

    command2 = [
            './plink2', 
            '--bfile', f'all_chromosomes',
            '--extract', f"./PRS.SNPs.{disease}", 
            '--make-bed',
            '--out', f'./b_{disease}/1000G_eur_PRS']

    try:
        result = subprocess.run(command2, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        pass 
    

    snp_file = f"./PRS.SNPs.{disease}"
    y = pd.read_csv(snp_file, header=None, sep="\t")
    y = y[y[0].str.len() > 0][0].tolist()
    dis = pd.read_csv(f'{disease}.txt', sep='\t')
    m = dis['SNPS'].isin(y)
    score = pd.DataFrame({
            'SNP': dis['SNPS'][m],
            'A1': dis['A1'][m],
            'BETA': dis['EFFECT'][m]})

    score.to_csv(f"./score_file_{disease}.txt", index=False, header=False, sep="\t", quoting=3)

    command3 = [
        './plink2', '--bfile', f'./b_{disease}/1000G_eur_PRS',
        '--out', f'./PRS_{disease}',
        '--score', f"./score_file_{disease}.txt", '1', '2', '3']

    try:
        result = subprocess.run(command3, check=True)
    except subprocess.CalledProcessError as e:
        pass #print(e.stdout)

    prs_score = pd.read_csv(f'PRS_{disease}.sscore', sep='\t')

    def prs_plot(dis):
        #prs = pd.read_csv(f'./PRS_{dis}.sscore', sep='\t')
        tiff = prs_score[prs_score['IID'] == 'Tiffany']['SCORE1_AVG'].values[0]
        sns.set(font_scale=1)
        plt.figure(figsize=(10, 6))
        plt.hist(prs_score['SCORE1_AVG'], bins=50, alpha=0.7, label='1000 Genomes')
        plt.axvline(tiff, color='r', linestyle='--', label='Tiffany')
        plt.xlabel('PRS Score')
        plt.ylabel('Frequency')
        t = dis.title()
        plt.title(f"PRS Distribution for {t} Disease")
        plt.legend()
        
        folder = f"./prs_graphs"
        file_path = os.makedirs(folder, exist_ok=True)
        
        plt.savefig(f"./prs_graphs/{dis}.png", dpi=300, bbox_inches='tight')
        plt.show()
        
    prs_plot(disease)
    return prs_score
```

In order to find the PRS and its distribution for each disease, simply  `prs` with GWAS `tsv` file.
```py
diabetes = prs(diabetes_gwas.tsv, 'diabetes')
heart = prs(heart_gwas.tsv, 'heart')
parkinsons = prs(parkinsons_gwas.tsv, 'parkinsons')
```





