#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-n', type=str, help='number of sampled sequences', required=True)

args = parser.parse_args()
n_sample = int(args.n)

# define paths
project_path = "/hpc/home/jm688/projects/AbMAP_JM"
disease_path = f"{project_path}/diseases_heavy"
model_weights_path = f"{project_path}/ablm/pretrained_models"
embedding_path = f"{disease_path}/Embeddings"
fasta_path = f"{disease_path}/FASTA"
plot_path = f"{project_path}/figures"
database_path = "/hpc/home/jm688/tools/ncbi-igblast-1.22.0/database"

def ReadFASTA(filename, region, as_pd=True):
    if filename.split(".")[-1] not in ["fasta","fna","fa"]:
        raise ValueError('Invalid file format. Expected formats are ["fasta","fna","fa"].')
    file_handle = open(filename,"r")
    seqs = []
    seqid = []
    tmp_seq = ""
    for line in file_handle:
        if (line[0] == ">"):
            if tmp_seq != "":
                seqs.append(tmp_seq.upper())
            seqid.append(line.split("\n")[0][1:])
            tmp_seq = ""
        else:
            tmp_seq+=line.split("\n")[0]
    seqs.append(tmp_seq.upper())
    file_handle.close()
    if as_pd:
        fasta = {}
        for i in range(len(seqs)):
            fasta[seqid[i]] = seqs[i]
        return pd.DataFrame(fasta,index=[f"{region.lower()}_sequence"]).transpose().reset_index()
    else:
        return seqs, seqid

database = {}    
for region in ["V", "D", "J"]:
    tmp = ReadFASTA(f"{database_path}/IGH{region}/Edited_IGH{region}.fasta", region)
    database[region.lower()] = tmp

del tmp


for disease in ['CMV', 'DEN']:
    # load IgBLAST results
    data = pd.read_csv(f"{fasta_path}/{disease}_{n_sample}_imputed_IgBLAST.tsv", sep="\t")
    data = data[['sequence_id', 'sequence',
                'v_call', 'v_sequence_start', 'v_sequence_end', 
                'v_germline_start', 'v_germline_end', 
                'v_germline_alignment', 
                'j_call', 'j_sequence_start', 'j_sequence_end', 
                'j_germline_start', 'j_germline_end', 
                'j_germline_alignment']]
    # merge with reference sequences
    for region in ['v', 'j']:
        data = data.merge(database[region], how='left', 
                          left_on=[f'{region}_call'], right_on=['index'])
        data[f'{region}_coverage'] = (data[f"{region}_germline_end"] - data[f"{region}_germline_start"])/data[f"{region}_sequence"].apply(len)
    # plot to check quality
    fig, axs = plt.subplots(1, 1, figsize=(5, 5))
    sns.kdeplot(data=data, x='v_coverage', y='j_coverage', ax=axs)
    axs.set_title(f"{disease} Alignment")
    plt.savefig(f"{plot_path}/{disease}_{n_sample}_after.png")
    plt.clf()

