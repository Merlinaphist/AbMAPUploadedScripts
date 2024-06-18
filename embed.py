import os
import pandas as pd
import numpy as np
import pickle
import torch
import abmap
import warnings
from tqdm import tqdm

import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-n', type=str, help='number of sampled sequences', required=True)

args = parser.parse_args()
n_sample = int(args.n)

warnings.filterwarnings("ignore")

# define paths
project_path = "/hpc/home/jm688/projects/AbMAP_JM"
disease_path = f"{project_path}/diseases_heavy"
model_weights_path = f"{project_path}/ablm/pretrained_models"
embedding_path = f"{disease_path}/Embeddings"
fasta_path = f"{disease_path}/FASTA"

chain_type = "H"
embed_type = "beplerberger"
num_mutations = 10
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
abmap_H = abmap.load_abmap(pretrained_path=f"{model_weights_path}/AbMAP_beplerberger_H.pt", 
                           plm_name="beplerberger")
abmap_H.to(device)

def EmbedOne(demo_seq):
    x = abmap.ProteinEmbedding(demo_seq, chain_type="H", embed_device=device)
    x.create_cdr_specific_embedding(embed_type="beplerberger", k=20) 
    with torch.no_grad():
        embed_fl = abmap_H.embed(x.cdr_embedding.unsqueeze(0).to(device), embed_type="fixed")
    return embed_fl.cpu().numpy().flatten().tolist()


embeddings = []
for disease in ["DEN", "CMV"]:
    data = pd.read_csv(f"{fasta_path}/{disease}_{n_sample}_imputed_IgBLAST.tsv", sep="\t")
    seqs = data[['sequence_aa', 'sequence_id']]
    embeds = {}
    for i in tqdm(range(len(seqs.values)), desc=disease):
        line = seqs.values[i]
        try:
            embed = EmbedOne(line[0])
            embeds[line[1]] = embed
        except:
            print(f"{line[1]} ValueError")
            continue
    embeds = pd.DataFrame(embeds).transpose()
    embeds["disease"] = disease
    embeddings.append(embeds)

embeddings = pd.concat(embeddings)
embeddings.to_csv(f"{embedding_path}/embeddings_{n_sample}_imputed.csv")
