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
#n_sample=1000

# define paths
project_path = "/hpc/home/jm688/projects/AbMAP_JM"
disease_path = f"{project_path}/diseases_heavy"
model_weights_path = f"{project_path}/ablm/pretrained_models"
embedding_path = f"{disease_path}/Embeddings"
fasta_path = f"{disease_path}/FASTA"
plot_path = f"{project_path}/figures"
database_path = "/hpc/home/jm688/tools/ncbi-igblast-1.22.0/database"

def LoadDatabase(filename, region, as_pd=True):
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
    fasta = {}
    for i in range(len(seqs)):
        fasta[seqid[i]] = seqs[i]
    fasta = pd.DataFrame(fasta,index=[f"{region.lower()}_germline_sequence"]).transpose().reset_index()
    fasta.columns = [f"{region.lower()}_call", f"{region.lower()}_germline_sequence"]
    return fasta

database = {}    
for region in ["V", "J"]:
    tmp = LoadDatabase(f"{database_path}/IGH{region}/Edited_IGH{region}.fasta", region)
    database[region.lower()] = tmp

del tmp

def LoadQueryFASTA(filename):
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
    fasta = {}
    for i in range(len(seqs)):
        fasta[seqid[i]] = seqs[i]
    fasta = pd.DataFrame(fasta,index=["full_query_sequence"]).transpose().reset_index()
    fasta.columns = ['query_id', 'full_query_sequence']
    return fasta


def LoadAlignment(disease, database):
    vdata = pd.read_csv(f"{fasta_path}/{disease}_{n_sample}_BLASTN_IGHV.tsv", sep="\t")
    jdata = pd.read_csv(f"{fasta_path}/{disease}_{n_sample}_BLASTN_IGHJ.tsv", sep="\t")
    vdata.columns = ['query_id', 'v_query_start', 'v_query_end', 'v_call', 
                    'v_germline_start', 'v_germline_end', 'v_evalue', 'v_pident', 
                    'v_query_frame', 'v_germline_frame']
    jdata.columns = ['query_id', 'j_query_start', 'j_query_end', 'j_call', 
                    'j_germline_start', 'j_germline_end', 'j_evalue', 'j_pident', 
                    'j_query_frame', 'j_germline_frame']
    data = vdata.merge(jdata, on=['query_id'])
    qseqs = LoadQueryFASTA(f"{fasta_path}/{disease}_{n_sample}.fasta")
    data = data.merge(qseqs, on=['query_id'])
    for region in ['v', 'j']:
        data = data.merge(database[region], how='left', on=[f'{region}_call'])
        data[f'{region}_coverage'] = (data[f"{region}_germline_end"] - data[f"{region}_germline_start"])/data[f"{region}_germline_sequence"].apply(len)
    return data


def impute(row):
    head = row['v_germline_sequence'][0:row['v_germline_start']-1]
    body = row['full_query_sequence'][row['v_query_start']-1:row['j_query_end']]
    tail = row['j_germline_sequence'][row['j_germline_end']:]
    return head + body + tail


for disease in ['CMV', 'DEN']:
    # load BLAST results
    data = LoadAlignment(disease, database)

    # plot to check quality
    fig, axs = plt.subplots(1, 1, figsize=(5, 5))
    sns.kdeplot(data=data, x='v_coverage', y='j_coverage', ax=axs)
    axs.set_title(f"{disease} Alignment")
    plt.savefig(f"{plot_path}/{disease}_{n_sample}_before.png")
    plt.clf()
    # impute
    data['imputed_query_sequence'] = data.apply(impute, axis=1)
    data.to_csv(f"{fasta_path}/{disease}_{n_sample}_imputed.csv", index=False)
    fasta = open(f"{fasta_path}/{disease}_{n_sample}_imputed.fasta", "w")
    for line in data[['query_id', 'imputed_query_sequence']].values:
        fasta.write(f">{line[0]}\n")
        fasta.write(f"{line[1]}\n")
    fasta.close()
