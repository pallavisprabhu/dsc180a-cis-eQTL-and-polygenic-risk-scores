import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os

def prs(gwas_file, disease, color='steelblue'):
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
        plt.hist(prs_score['SCORE1_AVG'], bins=50, alpha=0.7, label='1000 Genomes', color=color)
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