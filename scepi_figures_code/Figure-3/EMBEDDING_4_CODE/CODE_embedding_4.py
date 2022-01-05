------------------------------ scEpiSearch ---------------------------------------------------------------
import numpy as np
import pandas as pd
from sklearn.cluster import SpectralClustering
from itertools import cycle, islice
import matplotlib.pyplot as plt
from sklearn.neighbors import kneighbors_graph
import itertools
import networkx as nx

adj = np.loadtxt("./scepi_adj_embedding_2.txt",delimiter=",")

adj = np.nan_to_num(adj)
adj = np.abs(adj)                                   # no negative weights
adj = adj - np.diag(np.diag(adj))                   # no self loops
adj = np.triu(adj) + np.transpose(np.triu(adj))

ann = []
l = [np.repeat(1,822),np.repeat(2,1213),np.repeat(3,480),np.repeat(4,1621)]
l = list(itertools.chain(*l))
for i in range(adj.shape[0]):
    ann.append("D"+str(l[i])+" "+"Q"+str(i))
	

G = nx.from_numpy_matrix(np.matrix(adj>8))
# color=colors[clustering.labels_]
colors = np.array(list(islice(cycle(['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']),int(max(l)+1))))
# add black color for outliers (if any)
# colors = np.append(colors, ["#000000"])
color=colors[l]
pos=nx.spring_layout(G)

labels = {}    
for i,node in enumerate(G.nodes()):
    labels[node] = ann[i]

fig, ax = plt.subplots(figsize=(10,10))

nodes = nx.draw_networkx_nodes(G, pos=pos, ax=ax,node_color=color,cmap=plt.get_cmap('jet'),node_size=2)
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

p_1 = mpatches.Patch(color='#ff7f00', label='Human-Ex Neuron')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Forebrain')
p_3 = mpatches.Patch(color='#f781bf', label='Human-HSC')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-HSC')

plt.legend(handles=[p_1, p_2,p_3,p_4],fontsize=20)
# nx.draw_networkx_labels(G, pos, labels, font_size=9)
# nx.draw_networkx_edges(G, pos=pos, ax=ax)


------------------------------------------------------------ SCALE ------------------------------------------------------------------------------------------
#RAN FOLLOWING COMMAND FOR EACH QUERY DATASET TO GET FEATURES IN LATENT SPACE
python SCALE.py --dataset ./DATASET.tsv --n_centroids 1 --batch_size 500 --seed 43 --min_peaks 400 --lr 0.0002 --max_iter 500
######################################################################################

import matplotlib.pyplot
ts_mouse_hsc = np.loadtxt("./scale_data/result/feature_mouse_hsc.txt",delimiter="\t")
ts_human_hsc = np.loadtxt("./scale_data/result/feature_human_hsc.txt",delimiter="\t")
ts_mouse_neuron = np.loadtxt("./scale_data/result/feature_mouse_brain.txt",delimiter="\t")
ts_human_neuron = np.loadtxt("./scale_data/result/feature_human_brain.txt",delimiter="\t")
latent_space = np.concatenate([ts_human_neuron,ts_mouse_neuron,ts_human_hsc,ts_mouse_hsc],axis=0)
ts = TSNE(n_components=2).fit_transform(latent_space)
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
a1 = plt.scatter(ts[:804,0], ts[:804,1], marker='o', color='#ff7f00',s=0.5)
a2 = plt.scatter(ts[804:2017,0], ts[804:2017,1], marker='o', color='#4daf4a',s=0.5)
a3 = plt.scatter(ts[2017:2459,0], ts[2017:2459,1], marker='o', color='#f781bf',s=0.5)
a4 = plt.scatter(ts[2459:,0], ts[2459:,1], marker='o', color='#a65628',s=0.5)
p_1 = mpatches.Patch(color='#ff7f00', label='Human-Ex Neuron')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Forebrain')
p_3 = mpatches.Patch(color='#f781bf', label='Human-HSC')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-HSC')
plt.legend(handles=[p_1, p_2,p_3,p_4],fontsize=7)
plt.show()

------------------------------------------------------------ SCVI ------------------------------------------------------------------------------------------
import numpy as np
import scvi
import scanpy as sc

human_neuron = pd.read_csv(",./GE_GSE97942_human_brain.txt",sep=",", index_col = 0)
mouse_neuron= pd.read_csv("./GE_mouse_forebrain.txt",sep=",", index_col=0)
mouse_hsc = pd.read_csv("./GE_mouse_HSC.txt",sep=",",index_col = 0)
human_hsc = pd.read_csv("./GE_human_HSC.txt",sep=",", index_col = 0)

d = pd.concat([human_neuron,mouse_neuron,mouse_hsc,human_hsc],axis=1)
d = d.transpose()

d.to_csv("./combined_allcells_gene_enrichment_transposed.txt",sep=",",index=True,header=True)

data = scvi.data.read_csv("./combined_allcells_gene_enrichment_transposed.txt")
sc.pp.filter_genes(data, min_counts=3)

data.layers["counts"] = data.X.copy()
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
data.raw = data

#sc.pp.highly_variable_genes(data,n_top_genes=2000,subset=True,layer="counts",flavor="seurat_v3")

scvi.data.setup_anndata(data, layer="counts")
model = scvi.model.SCVI(data)

model.train()

latent = model.get_latent_representation()
np.savetxt("./scvi_latent.txt",delimiter=",")


import matplotlib.pyplot
latent_space = np.loadtxt("./scvi_latent.txt",delimiter=",")
import numpy as np
from sklearn.manifold import TSNE
ts = TSNE(n_components=2).fit_transform(latent_space)
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
a1 = plt.scatter(ts[:822,0], ts[:822,1], marker='o', color='#ff7f00',s=0.5)
a2 = plt.scatter(ts[822:2035,0], ts[822:2035,1], marker='o', color='#4daf4a',s=0.5)
a3 = plt.scatter(ts[2035:3656,0], ts[2035:3656,1], marker='o', color='#f781bf',s=0.5)
a4 = plt.scatter(ts[3656:,0], ts[3656:,1], marker='o', color='#a65628',s=0.5)
p_1 = mpatches.Patch(color='#ff7f00', label='Human-Ex Neuron')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Forebrain')
p_3 = mpatches.Patch(color='#f781bf', label='Human-HSC')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-HSC')
plt.legend(handles=[p_1, p_2,p_3,p_4],fontsize=7)
plt.show()

------------------------------------------------------------ SCANORAMA ------------------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import scanorama

human_neuron = pd.read_csv("./GE_GSE97942_human_brain.txt",sep=",", index_col = 0)
mouse_neuron= pd.read_csv("./GE_mouse_forebrain.txt",sep=",", index_col=0)
mouse_hsc = pd.read_csv("./GE_mouse_HSC.txt",sep=",",index_col = 0)
human_hsc = pd.read_csv("./GE_human_HSC.txt",sep=",", index_col = 0)

gene_list = [np.array(human_neuron.index),np.array(human_neuron.index),np.array(human_neuron.index),np.array(human_neuron.index)]

human_neuron = np.transpose(np.array(human_neuron))
mouse_neuron = np.transpose(np.array(mouse_neuron))
human_hsc = np.transpose(np.array(human_hsc))
mouse_hsc = np.transpose(np.array(mouse_hsc))

dataset = [human_neuron,mouse_neuron,mouse_hsc,human_hsc]

integrated, corrected, genes = scanorama.correct(dataset, gene_list, return_dimred=True)

integrated_f = np.concatenate([integrated[0],integrated[1],integrated[2],integrated[3],integrated[4]],axis=0)

np.savetxt("./scanorama.txt",np.transpose(integrated_f),delimiter=",")

import matplotlib.pyplot
latent_space = np.loadtxt("./scanorama.txt",delimiter=",")
import numpy as np
from sklearn.manifold import TSNE
ts = TSNE(n_components=2).fit_transform(np.transpose(latent_space))
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
a1 = plt.scatter(ts[:822,0], ts[:822,1], marker='o', color='#ff7f00',s=0.5)
a2 = plt.scatter(ts[822:2035,0], ts[822:2035,1], marker='o', color='#4daf4a',s=0.5)
a3 = plt.scatter(ts[2035:3656,0], ts[2035:3656,1], marker='o', color='#f781bf',s=0.5)
a4 = plt.scatter(ts[3656:,0], ts[3656:,1], marker='o', color='#a65628',s=0.5)
p_1 = mpatches.Patch(color='#ff7f00', label='Human-Ex Neuron')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Forebrain')
p_3 = mpatches.Patch(color='#f781bf', label='Human-HSC')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-HSC')
plt.legend(handles=[p_1, p_2,p_3,p_4],fontsize=7)
plt.show()

---------------------------------------------------------- MINT --------------------------------------------------------
library(mixOmics)
human_neuron = read.table("./GE_GSE97942_human_brain.txt",sep=",", header=TRUE, row.names = 1)
mouse_neuron= read.table("./GE_mouse_forebrain.txt",sep=",", header=TRUE, row.names = 1)
mouse_hsc = read.table("./GE_mouse_HSC.txt",sep=",", header=TRUE, row.names = 1)
human_hsc = read.table("./GE_human_HSC.txt",sep=",", header=TRUE, row.names = 1)

colnames(human_neuron)[1:ncol(human_neuron)]  = 'Ex. Neuron'
colnames(mouse_neuron)[1:ncol(mouse_neuron)]  = 'Forebrain'
colnames(human_hsc)[1:ncol(human_hsc)]  = 'Human HSC'
colnames(mouse_hsc)[1:ncol(mouse_hsc)]  =  'Mouse HSC'

combined = cbind(human_neuron,mouse_neuron,human_hsc,mouse_hsc)

X = t(combined)
Y = colnames(combined)
Y = as.factor(Y)

colnames(human_neuron)[1:ncol(human_neuron)]  = 'Human'
colnames(mouse_neuron)[1:ncol(mouse_neuron)]  = 'Mouse'
colnames(human_hsc)[1:ncol(human_hsc)]  = 'Human'
colnames(mouse_hsc)[1:ncol(mouse_hsc)]  = 'Mouse'


combined = cbind(human_neuron,mouse_neuron,human_hsc,mouse_hsc)

study = colnames(combined)
study = as.factor(study)

rownames(X) = 1:ncol(combined)

mint.plsda.res = mint.plsda(X = X, Y = Y, study = study, ncomp = 2)
#mint.plsda.res # lists the different functions
color = c('#ff7f00','#4daf4a','#f781bf','#a65628')
plotIndiv(mint.plsda.res, legend = TRUE, title = '', subtitle = '',col=color)



