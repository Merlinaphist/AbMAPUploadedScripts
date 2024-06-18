import numpy as np
import pandas as pd
import umap.umap_ as umap
import matplotlib.pyplot as plt
import seaborn as sns
project_path = "/hpc/home/jm688/projects/AbMAP_JM"
disease_path = f"{project_path}/diseases_heavy"
model_weights_path = f"{project_path}/ablm/pretrained_models"
embedding_path = f"{disease_path}/Embeddings"
plot_path = f"{project_path}/figures"
nrows = [1000, 10000]

xmin=-2
xmax=13
xticks = np.arange(xmin, xmax+1)
ymin=-2
ymax=11
yticks = np.arange(ymin, ymax+1)
fig, axs = plt.subplots(3, 2, figsize=(15, 15))
for col in range(2):
    i = nrows[col]
    embed = pd.read_csv(f"{embedding_path}/embeddings_{i}_imputed.csv", index_col=0)
    embed = embed.sort_values(by=["disease"])
    print(embed.value_counts("disease").to_string())
    fit = umap.UMAP()
    u = fit.fit_transform(embed.iloc[:,0:512])
    u = pd.DataFrame(u, index=embed.index, columns=['x', 'y'])
    u["disease"] = embed['disease']
    # Scatterplot
    sns.scatterplot(data = u, x='x', y='y', hue='disease', ax=axs[0][col])
    axs[0][col].set_title(f"{i} Antibodies Scatterplot")
    axs[0][col].set_xticks(xticks)
    axs[0][col].set(xlim=(xmin, xmax))
    axs[0][col].set_yticks(yticks)
    axs[0][col].set(ylim=(ymin, ymax))
    # CMV Density
    sns.kdeplot(data = u.query("disease=='CMV'"), x='x', y='y', 
                fill=True, ax=axs[1][col])
    axs[1][col].set_title(f"{embed.value_counts('disease')['CMV']}/{i} Antibodies CMV Density")
    axs[1][col].set_xticks(xticks)
    axs[1][col].set(xlim=(xmin, xmax))
    axs[1][col].set_yticks(yticks)
    axs[1][col].set(ylim=(ymin, ymax))
    # Dengue Density
    sns.kdeplot(data = u.query("disease=='DEN'"), x='x', y='y',
                fill=True, color=sns.color_palette("tab10")[1], ax=axs[2][col])
    axs[2][col].set_title(f"{embed.value_counts('disease')['DEN']}/{i} Antibodies Dengue Density")
    axs[2][col].set_xticks(xticks)
    axs[2][col].set(xlim=(xmin, xmax))
    axs[2][col].set_yticks(yticks)
    axs[2][col].set(ylim=(ymin, ymax))

plt.savefig(f"{plot_path}/embeddings_imputed.png")
plt.clf()