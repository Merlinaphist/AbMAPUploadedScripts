#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import pickle
import torch
import abmap
import warnings
from tqdm import tqdm

warnings.filterwarnings("ignore")

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

embeddings = []
for disease in ["DEN", "CMV"]:
    path = f"{disease_path}/{disease}"
    files = os.listdir(path)
    seqs = []
    for file in files:
        seq = pd.read_csv(f"{path}/{file}", compression='gzip', skiprows=1)
        seq = seq[['sequence']]
        seq["header"] = [f"{file.split('.')[0]}|{i}" for i in range(seq.shape[0])]
        seqs.append(seq)
    seqs = pd.concat(seqs)
    seqs = seqs.sample(n=n_sample, random_state=0)
    fasta = open(f"{fasta_path}/{path.split('/')[-1]}_{n_sample}.fasta", "w")
    for line in seqs.values:
        fasta.write(f">{line[1]}\n")
        fasta.write(f"{line[0]}\n")
    fasta.close()