---------------------------------------------------------------scEpiSearch----------------------------------------------------------
import numpy as np
import pandas as pd
from sklearn.cluster import SpectralClustering
from itertools import cycle, islice
import matplotlib.pyplot as plt
from sklearn.neighbors import kneighbors_graph

adj = np.loadtxt("./scepi_adj_embedding_2.txt",delimiter=",")

adj = np.nan_to_num(adj)
adj = np.abs(adj)                                   # no negative weights
adj = adj - np.diag(np.diag(adj))                   # no self loops
adj = np.triu(adj) + np.transpose(np.triu(adj))

import networkx as nx
G = nx.from_numpy_matrix(np.matrix(adj>7))
# color=colors[clustering.labels_]
colors = np.array(list(islice(cycle(['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']),int(max(l)+1))))
# add black color for outliers (if any)
# colors = np.append(colors, ["#000000"])
color=colors[l]
pos=nx.spring_layout(G)

labels = {}    
for i,node in enumerate(G.nodes()):
    labels[node] = ann[i]

fig, ax = plt.subplots(figsize=(5,5))

nodes = nx.draw_networkx_nodes(G, pos=pos, ax=ax,node_color=color,cmap=plt.get_cmap('jet'),node_size=6)
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

p_1 = mpatches.Patch(color='#ff7f00', label='Human-GM')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Bcell')
p_3 = mpatches.Patch(color='#f781bf', label='Human-Tcell')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-Tcell')

plt.legend(handles=[p_1, p_2,p_3,p_4],fontsize=10)
-------------------------------------------------------------------------SCALE-----------------------------------------------------------
#RAN FOLLOWING COMMAND FOR EACH QUERY DATASET TO GET FEATURES IN LATENT SPACE
python SCALE.py --dataset ./DATASET.tsv --n_centroids 1 --batch_size 500 --seed 43 --min_peaks 400 --lr 0.0002 --max_iter 500
######################################################################################

import matplotlib.pyplot
human_bcell = np.loadtxt("./output_scale_embedding3/feature_human_gm.txt",delimiter="\t")
mouse_bcell = np.loadtxt("./output_scale_embedding3/feature_gm_mouse.txt",delimiter="\t")
human_tcell = np.loadtxt("./output_scale_embedding4/feature_tcell_human.txt",delimiter="\t")
mouse_tcell = np.loadtxt("./output_scale_embedding4/feature_tcell_mouse.txt",delimiter="\t")
latent_space = np.concatenate([human_bcell,mouse_bcell,human_tcell,mouse_tcell],axis=0)
ts = TSNE(n_components=2).fit_transform(latent_space)
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
a1 = plt.scatter(ts[:95,0], ts[:95,1], marker='o', color='#ff7f00',s=3)
a2 = plt.scatter(ts[95:194,0], ts[95:194,1], marker='o', color='#4daf4a',s=3)
a3 = plt.scatter(ts[194:264,0], ts[194:264,1], marker='o', color='#f781bf',s=3)
a4 = plt.scatter(ts[264:,0], ts[264:,1], marker='o', color='#a65628',s=3)
p_1 = mpatches.Patch(color='#ff7f00', label='Human-GM12878')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Bcell')
p_3 = mpatches.Patch(color='#f781bf', label='Human-Tcell')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-Tcell')
plt.legend(handles=[p_1, p_2,p_3,p_4],fontsize=9)
plt.show()

---------------------------------------------------------------SCVI-----------------------------------------------------------------------------
import numpy as np
import scvi
import scanpy as sc
import pandas as pd

human_bcell= pd.read_csv("./GE_human_gm100cells.txt",sep=",", index_col = 0)
mouse_bcell= pd.read_csv("./GE_mouse_bcell100cells.txt",sep=",", index_col=0)
mouse_proximal  = pd.read_csv("./GE_mouse_tcell100cells.txt",sep=",",index_col = 0)
human_hek = pd.read_csv("./GE_human_tcell100cells.txt",sep=",", index_col = 0)
human_bcell = pd.DataFrame(human_bcell)
human_hek = pd.DataFrame(human_hek)
mouse_bcell = pd.DataFrame(mouse_bcell)
mouse_proximal = pd.DataFrame(mouse_proximal)

d = pd.concat([human_bcell,mouse_bcell,human_hek,mouse_proximal],axis=1)
d = d.transpose()

d.to_csv("./combined_allcells_gene_enrichment_transposed.txt",sep=",",index=True,header=True)

data = scvi.data.read_csv("./combined_allcells_gene_enrichment_transposed.txt")
sc.pp.filter_genes(data, min_counts=300)

data.layers["counts"] = data.X.copy()
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
data.raw = data

#sc.pp.highly_variable_genes(data,n_top_genes=2000,subset=True,layer="counts",flavor="seurat_v3")

scvi.data.setup_anndata(data)
model = scvi.model.SCVI(data)

model.train()

latent = model.get_latent_representation()

import matplotlib.pyplot
from sklearn.manifold import TSNE
latent_space = np.loadtxt("./scvi_embedding4.txt",delimiter=",")
ts = TSNE(n_components=2).fit_transform(latent_space)
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
a1 = plt.scatter(ts[:100,0], ts[:100,1], marker='o', color='#f781bf',s=3)
a2 = plt.scatter(ts[100:200,0], ts[100:200,1], marker='o', color='#a65628',s=3)
a3 = plt.scatter(ts[200:300,0], ts[200:300,1], marker='o', color='#4daf4a',s=3)
a4 = plt.scatter(ts[300:,0], ts[300:,1], marker='o', color='#ff7f00',s=3)
p_1 = mpatches.Patch(color='#ff7f00', label='Human-GM12878')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Bcell')
p_3 = mpatches.Patch(color='#f781bf', label='Human-Tcell')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-Tcell')
plt.legend(handles=[p_1, p_2,p_3,p_4],fontsize=9)
plt.show()

---------------------------------------------------------------------------SCANORAMA---------------------------------------------------
import numpy as np
import pandas as pd
import scanorama

human_bcell= pd.read_csv("./GE_human_gm100cells.txt",sep=",", index_col = 0)
mouse_bcell= pd.read_csv("./GE_mouse_bcell100cells.txt",sep=",", index_col=0)
mouse_proximal  = pd.read_csv("./GE_mouse_tcell100cells.txt",sep=",",index_col = 0)
human_hek = pd.read_csv("./GE_human_tcell100cells.txt",sep=",", index_col = 0)

gene_list = [np.array(human_bcell.index),np.array(mouse_bcell.index),np.array(human_hek.index),np.array(mouse_proximal.index)]

human_bcell = np.transpose(np.array(human_bcell))
mouse_bcell = np.transpose(np.array(mouse_bcell))
human_hek = np.transpose(np.array(human_hek))
mouse_proximal = np.transpose(np.array(mouse_proximal))
human_bcell = human_bcell.astype('float32')
human_hek = human_hek.astype('float32')
mouse_bcell = mouse_bcell.astype('float32')
mouse_proximal = mouse_proximal.astype('float32')
human_bcell = np.nan_to_num(human_bcell)
human_hek = np.nan_to_num(human_hek)
mouse_bcell = np.nan_to_num(mouse_bcell)
mouse_proximal = np.nan_to_num(mouse_proximal)


dataset = [human_bcell,mouse_bcell,human_hek,mouse_proximal]

integrated, corrected, genes = scanorama.correct(dataset, gene_list, return_dimred=True)

integrated_f = np.concatenate([integrated[0],integrated[1],integrated[2],integrated[3],integrated[4]],axis=0)

np.savetxt("./scanorama.txt",np.transpose(integrated_f),delimiter=",")

import matplotlib.pyplot
latent_space = np.loadtxt("./scanorama_embedding4.txt",delimiter=",")
import numpy as np
from sklearn.manifold import TSNE
ts = TSNE(n_components=2).fit_transform(np.transpose(latent_space))
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
a1 = plt.scatter(ts[:100,0], ts[:100,1], marker='o', color='#ff7f00',s=3)
a2 = plt.scatter(ts[100:200,0], ts[100:200,1], marker='o', color='#4daf4a',s=3)
a3 = plt.scatter(ts[200:300,0], ts[200:300,1], marker='o', color='#f781bf',s=3)
a4 = plt.scatter(ts[300:,0], ts[300:,1], marker='o', color='#a65628',s=3)
p_1 = mpatches.Patch(color='#ff7f00', label='Human-GM12878')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Bcell')
p_3 = mpatches.Patch(color='#f781bf', label='Human-Tcell')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-Tcell')
plt.legend(handles=[p_1, p_2,p_3,p_4],fontsize=9)
plt.show()

-------------------------------------------------------------------MINT---------------------------------------------------------------------------
library(mixOmics)
human_bcell = read.table("./GE_human_gm100cells.txt",sep=",", header=TRUE, row.names = 1)
mouse_bcell= read.table("./GE_mouse_bcell100cells.txt",sep=",", header=TRUE, row.names = 1)
mouse_proximal = read.table("./GE_mouse_tcell100cells.txt",sep=",", header=TRUE, row.names = 1)
human_hek = read.table("./GE_human_tcell100cells.txt",sep=",", header=TRUE, row.names = 1)

colnames(human_bcell)[1:ncol(human_bcell)]  = 'GM12878'
colnames(mouse_bcell)[1:ncol(mouse_bcell)]  = 'Bcell'
colnames(human_hek)[1:ncol(human_hek)]  = 'Tcell'
colnames(mouse_proximal)[1:ncol(mouse_proximal)]  =  'T cell'

combined = cbind(human_bcell,mouse_bcell,human_hek,mouse_proximal)

X = t(combined)
Y = colnames(combined)
Y = as.factor(Y)

colnames(human_bcell)[1:ncol(human_bcell)]  = 'Human'
colnames(mouse_bcell)[1:ncol(mouse_bcell)]  = 'Mouse'
colnames(human_hek)[1:ncol(human_hek)]  = 'Human'
colnames(mouse_proximal)[1:ncol(mouse_proximal)]  = 'Mouse'


combined = cbind(human_bcell,mouse_bcell,human_hek,mouse_proximal)

study = colnames(combined)
study = as.factor(study)

rownames(X) = 1:ncol(combined)

mint.plsda.res = mint.plsda(X = X, Y = Y, study = study, ncomp = 2)
#mint.plsda.res # lists the different functions
color = c('#4daf4a','#ff7f00','#a65628','#f781bf')
plotIndiv(mint.plsda.res, legend = TRUE, title = '', subtitle = '',col=color)
