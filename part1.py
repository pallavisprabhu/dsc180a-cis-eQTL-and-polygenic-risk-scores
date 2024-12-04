import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os

def ciseQTL(chromosome): #chromosome as an integer
    def extract_protein():
        expression_df = pd.read_csv('GD462.GeneQuantRPKM.50FN.samplename.resk10.txt', delim_whitespace=True)
        annotations = pd.read_csv('gene_annot.txt', delim_whitespace=True)
        protein_coding = annotations[annotations['TYPE'] == 'protein_coding']
        intersection = np.intersect1d(np.array(expression_df['TargetID'].apply(lambda x: x.split('.')[0])), np.array(protein_coding['SYM']))
        protein_exp = expression_df[expression_df['TargetID'].apply(lambda x: x.split('.')[0]).apply(lambda x: x in intersection)]
        return protein_exp, protein_coding
    
    protein_exp, protein_coding = extract_protein()
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
    print(f"Summary statistics filed saved to './ciseqtls_sumstats/sumstats_{chromosome}'")
    
    return results_df

def cis_plot(sum_stats_path, gene):
    sum_stats = pd.read_csv(sum_stats_path, sep='\t')
    chr = sum_stats['Chr'].iloc[0]
    filtered = sum_stats[sum_stats['Gene'] == gene].reset_index()


    sig = " ".join(list(sum_stats[(sum_stats['Gene'] == gene) & (sum_stats['P-value'] < 0.05)]['SNP'].values))
    print(f"Statistically significant SNPs for {gene} on Chromsome {chr} are: {sig}")

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
       