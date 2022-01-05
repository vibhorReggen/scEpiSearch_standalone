---------------------------------------------- scEpiSearch ----------------------------------------------------------------------------
import numpy as np
import pandas as pd
from sklearn.cluster import SpectralClustering
from itertools import cycle, islice
import matplotlib.pyplot as plt
from sklearn.neighbors import kneighbors_graph

adj = np.loadtxt("./scepi_adj_embedding_1.txt",delimiter=",")

adj = np.nan_to_num(adj)
adj = np.abs(adj)                                   # no negative weights
adj = adj - np.diag(np.diag(adj))                   # no self loops
adj = np.triu(adj) + np.transpose(np.triu(adj))

import itertools
ann = []
l = [np.repeat(1,30),np.repeat(2,30),np.repeat(3,10),np.repeat(4,12),np.repeat(5,11),np.repeat(6,14),np.repeat(7,10),np.repeat(8,10)]
l = list(itertools.chain(*l))
for i in range(adj.shape[0]):
    ann.append("D"+str(l[i])+" "+"Q"+str(i))
	import networkx as nx
G = nx.from_numpy_matrix(np.matrix(adj>7))
# color=colors[clustering.labels_]
colors = np.array(list(islice(cycle(['#984ea3','#a65628','#f781bf','#ff7f00','#4daf4a','#e41a1c','#999999','#dede00']),int(max(l)+1))))
# add black color for outliers (if any)
# colors = np.append(colors, ["#000000"])
color=colors[l]
pos=nx.spring_layout(G)

labels = {}    
for i,node in enumerate(G.nodes()):
    labels[node] = ann[i]

fig, ax = plt.subplots(figsize=(5,5))

nodes = nx.draw_networkx_nodes(G, pos=pos, ax=ax,node_color=color,cmap=plt.get_cmap('jet'),node_size=30)
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

p_1 = mpatches.Patch(color='#a65628', label='Mouse-HSC')
p_2 = mpatches.Patch(color='#f781bf', label='Human-HSC')
p_3 = mpatches.Patch(color='#ff7f00', label='Human-Neuron')
p_4 = mpatches.Patch(color='#4daf4a', label='Mouse-Neuron')
p_5 = mpatches.Patch(color='#e41a1c', label='Mouse-Bcell')
p_6 = mpatches.Patch(color='#999999', label='Human-GM')
p_7 = mpatches.Patch(color='#dede00', label='Human-GM-Different Batch')
p_8 = mpatches.Patch(color='#984ea3', label='Human-Myoblast')

plt.legend(handles=[p_1, p_2,p_3,p_4,p_5,p_6,p_7,p_8],fontsize=7)
# nx.draw_networkx_labels(G, pos, labels, font_size=9)
# nx.draw_networkx_edges(G, pos=pos, ax=ax)

--------------------------------------------------SCALE -----------------------------------------------------------------------------------------------
#RAN FOLLOWING COMMAND FOR EACH QUERY DATASET TO GET FEATURES IN LATENT SPACE
python SCALE.py --dataset ./DATASET.tsv --n_centroids 1 --batch_size 500 --seed 43 --min_peaks 400 --lr 0.0002 --max_iter 500
######################################################################################

import matplotlib.pyplot
ts_mouse_hsc = np.loadtxt("./output_mouse_hsc/feature.txt",delimiter="\t")
ts_human_hsc = np.loadtxt("./output_human_hsc_SCALE/feature.txt",delimiter="\t")
ts_mouse_neuron = np.loadtxt("./output_mouse_forebrain_SCALE/feature.txt",delimiter="\t")
ts_human_neuron = np.loadtxt("./output_human_neuron/feature.txt",delimiter="\t")
ts_mouse_gm = np.loadtxt("./output_bcell_mouse_SCALE/feature.txt",delimiter="\t")
ts_human_gm = np.loadtxt("./output_human_GM_SCALE/feature.txt",delimiter="\t")
ts_human_gm_gse68103 = np.loadtxt("./output_human_GM_GSE68103_SCALE/feature.txt",delimiter="\t")
ts_human_myoblast = np.loadtxt("./output_human_myoblast_SCALE/feature.txt",delimiter="\t")
latent_space = np.concatenate([ts_human_neuron,ts_mouse_neuron,ts_human_hsc,ts_mouse_hsc,ts_mouse_gm,ts_human_gm,ts_human_gm_gse68103,ts_human_myoblast],axis=0)
ts = TSNE(n_components=2).fit_transform(latent_space)
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
a1 = plt.scatter(ts[:13,0], ts[:13,1], marker='o', color='#ff7f00')
a2 = plt.scatter(ts[13:23,0], ts[13:23,1], marker='o', color='#4daf4a')
a3 = plt.scatter(ts[23:53,0], ts[23:53,1], marker='o', color='#f781bf')
a4 = plt.scatter(ts[53:84,0], ts[53:84,1], marker='o', color='#a65628')
a5 = plt.scatter(ts[84:95,0], ts[84:95,1], marker='o', color='#e41a1c')
a6 = plt.scatter(ts[95:109,0], ts[95:109,1], marker='o', color='#999999')
a7 = plt.scatter(ts[109:115,0], ts[109:115,1], marker='o', color='#dede00')
a8 = plt.scatter(ts[115:,0], ts[115:,1], marker='o', color='#984ea3')
p_1 = mpatches.Patch(color='#ff7f00', label='Human-Ex Neuron')
p_2 = mpatches.Patch(color='#4daf4a', label='Mouse-Neuron')
p_3 = mpatches.Patch(color='#f781bf', label='Human-HSC')
p_4 = mpatches.Patch(color='#a65628', label='Mouse-HSC')
p_5 = mpatches.Patch(color='#e41a1c', label='Mouse-Bcell')
p_6 = mpatches.Patch(color='#999999', label='Human-GM')
p_7 = mpatches.Patch(color='#dede00', label='Human-GM-Different Batch')
p_8 = mpatches.Patch(color='#984ea3', label='Human Myoblast')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(handles=[p_1, p_2,p_3,p_4,p_5,p_6,p_7,p_8], prop={'size': 10},loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

------------------------------------------------------------ SCVI ------------------------------------------------------------------------------------------
import numpy as np
import scvi
import scanpy as sc

human_neuron = pd.read_csv("./human_neuron.txt",sep=",", index_col = 0)
mouse_neuron= pd.read_csv("./mouse_neuron.txt",sep=",", index_col = 0)
mouse_bcell = pd.read_csv("./mouse_bcells.txt",sep=",", index_col = 0)
human_gm = pd.read_csv("./human_gm.txt",sep=",", index_col = 0)
human_gm_gse68103 = pd.read_csv("./human_gm_GSE68103.txt",sep=",",index_col = 0)
human_hsc = pd.read_csv("./human_hsc.txt",sep=",", index_col = 0)
mouse_hsc = pd.read_csv("./mouse_hsc.txt",sep=",", index_col = 0)
human_myoblast = pd.read_csv("./human_myoblast.txt",sep=",",index_col = 0)

d = pd.concat([human_neuron,mouse_neuron,mouse_hsc,human_hsc,mouse_bcell,human_gm,human_myoblast],axis=1)
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
latent_space = np.loadtxt("./gene_enrichment_scores_embedding_data/scvi_latent_Space.txt",delimiter=",")
ts = TSNE(n_components=2).fit_transform(latent_space)
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
a1 = plt.scatter(ts[:30,0], ts[:30,1], marker='o', color='#f781bf')
a2 = plt.scatter(ts[30:61,0], ts[30:61,1], marker='o', color='#a65628')
a3 = plt.scatter(ts[61:71,0], ts[61:71,1], marker='o', color='#4daf4a')
a4 = plt.scatter(ts[71:83,0], ts[71:83,1], marker='o', color='#ff7f00')
a5 = plt.scatter(ts[83:94,0], ts[83:94,1], marker='o', color='#e41a1c')
a6 = plt.scatter(ts[94:108,0], ts[94:108,1], marker='o', color='#999999')
a7 = plt.scatter(ts[108:119,0], ts[108:119,1], marker='o', color='#dede00')
a8 = plt.scatter(ts[119:,0], ts[119:,1], marker='o', color='#984ea3')
p_1 = mpatches.Patch(color='#f781bf', label='Human HSC')
p_2 = mpatches.Patch(color='#a65628', label='Mouse HSC')
p_3 = mpatches.Patch(color='#4daf4a', label='Mouse Neuron')
p_4 = mpatches.Patch(color='#ff7f00', label='Human Neuron')
p_5 = mpatches.Patch(color='#e41a1c', label='Mouse-Bcell')
p_6 = mpatches.Patch(color='#999999', label='Human-GM')
p_7 = mpatches.Patch(color='#dede00', label='Human-GM-Different Batch')
p_8 = mpatches.Patch(color='#984ea3', label='Human Myoblast')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(handles=[p_1, p_2,p_3,p_4,p_5,p_6,p_7,p_8], prop={'size': 10},loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

----------------------------------------------------SCANORAMA ------------------------------------------------------------------------
import numpy as np
import pandas as pd
import scanorama
human_neuron = pd.read_csv("./human_neuron.txt",sep=",", index_col = 0)
mouse_neuron= pd.read_csv("./mouse_neuron.txt",sep=",", index_col = 0)
mouse_bcell = pd.read_csv("./mouse_bcells.txt",sep=",", index_col = 0)
human_gm = pd.read_csv("./human_gm.txt",sep=",", index_col = 0)
human_gm_gse68103 = pd.read_csv("./human_gm_GSE68103.txt",sep=",",index_col = 0)
human_hsc = pd.read_csv("./human_hsc.txt",sep=",", index_col = 0)
mouse_hsc = pd.read_csv("./mouse_hsc.txt",sep=",", index_col = 0)
human_myoblast = pd.read_csv("./human_myoblast.txt",sep=",",index_col = 0)

gene_list = [np.array(human_neuron.index),np.array(human_neuron.index),np.array(human_neuron.index),np.array(human_neuron.index),np.array(human_neuron.index),np.array(human_neuron.index),np.array(human_neuron.index),np.array(human_neuron.index)]

human_neuron = np.transpose(np.array(human_neuron))
mouse_neuron = np.transpose(np.array(mouse_neuron))
mouse_bcell = np.transpose(np.array(mouse_bcell))
human_gm = np.transpose(np.array(human_gm))
human_gm_gse68103 = np.transpose(np.array(human_gm_gse68103))
human_hsc = np.transpose(np.array(human_hsc))
mouse_hsc = np.transpose(np.array(mouse_hsc))
human_myoblast = np.transpose(np.array(human_myoblast))

dataset = [mouse_hsc,human_hsc,mouse_neuron,human_neuron,mouse_bcell,human_gm,human_gm_gse68103,human_myoblast]

integrated, corrected, genes = scanorama.correct(dataset, gene_list, return_dimred=True)

integrated_f = np.concatenate([integrated[0],integrated[1],integrated[2],integrated[3],integrated[4],integrated[5],integrated[6],integrated[7]],axis=0)

np.savetxt("./scanorama.txt",np.transpose(integrated_f),delimiter=",")

import matplotlib.pyplot
latent_space = np.loadtxt("./scanorama.txt",delimiter=",")
ts = TSNE(n_components=2).fit_transform(latent_space)
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
a1 = plt.scatter(ts[:30,0], ts[:30,1], marker='o', color='#f781bf')
a2 = plt.scatter(ts[30:61,0], ts[30:61,1], marker='o', color='#a65628')
a3 = plt.scatter(ts[61:71,0], ts[61:71,1], marker='o', color='#4daf4a')
a4 = plt.scatter(ts[71:83,0], ts[71:83,1], marker='o', color='#ff7f00')
a5 = plt.scatter(ts[83:94,0], ts[83:94,1], marker='o', color='#e41a1c')
a6 = plt.scatter(ts[94:108,0], ts[94:108,1], marker='o', color='#999999')
a7 = plt.scatter(ts[108:119,0], ts[108:119,1], marker='o', color='#dede00')
a8 = plt.scatter(ts[119:,0], ts[119:,1], marker='o', color='#984ea3')
p_1 = mpatches.Patch(color='#f781bf', label='Human HSC')
p_2 = mpatches.Patch(color='#a65628', label='Mouse HSC')
p_3 = mpatches.Patch(color='#4daf4a', label='Mouse Neuron')
p_4 = mpatches.Patch(color='#ff7f00', label='Human Neuron')
p_5 = mpatches.Patch(color='#e41a1c', label='Mouse-Bcell')
p_6 = mpatches.Patch(color='#999999', label='Human-GM')
p_7 = mpatches.Patch(color='#dede00', label='Human-GM-Different Batch')
p_8 = mpatches.Patch(color='#984ea3', label='Human Myoblast')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
plt.legend(handles=[p_1, p_2,p_3,p_4,p_5,p_6,p_7,p_8], prop={'size': 10},loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

---------------------------------------------------------------------------MINT---------------------------------------------------------------------
library(mixOmics)
human_neuron = read.table("./human_neuron.txt",sep=",", header=TRUE, row.names = 1)
mouse_neuron= read.table("./mouse_neuron.txt",sep=",", header=TRUE, row.names = 1)
mouse_bcell = read.table("./mouse_bcells.txt",sep=",", header=TRUE, row.names = 1)
human_gm = read.table("./human_gm.txt",sep=",", header=TRUE, row.names = 1)
human_gm_gse68103 = read.table("./human_gm_GSE68103.txt",sep=",", header=TRUE, row.names = 1)
human_hsc = read.table("./human_hsc.txt",sep=",", header=TRUE, row.names = 1)
mouse_hsc = read.table("./mouse_hsc.txt",sep=",", header=TRUE, row.names = 1)
human_myoblast = read.table("./human_myoblast.txt",sep=",", header=TRUE, row.names = 1)

colnames(human_neuron)[1:ncol(human_neuron)]  = 'Human neuron'
colnames(mouse_neuron)[1:ncol(mouse_neuron)]  = 'Mouse neuron'
colnames(human_myoblast)[1:ncol(human_myoblast)]  = 'Human myoblast'
colnames(human_gm_gse68103)[1:ncol(human_gm_gse68103)]  =  'Human bcell-1'
colnames(human_gm)[1:ncol(human_gm)]  = 'Human bcell-2'
colnames(mouse_bcell)[1:ncol(mouse_bcell)]  = 'Mouse bcell'
colnames(human_hsc)[1:ncol(human_hsc)]  = 'Human HSC'
colnames(mouse_hsc)[1:ncol(mouse_hsc)]  = 'Mouse HSC'

combined = cbind(human_neuron,mouse_neuron,human_myoblast,human_gm_gse68103,human_gm,mouse_bcell,human_hsc,mouse_hsc)

X = t(combined)
Y = colnames(combined)
Y = as.factor(Y)

colnames(human_neuron)[1:ncol(human_neuron)]  = 'Human'
colnames(mouse_neuron)[1:ncol(mouse_neuron)]  = 'Mouse'
colnames(human_myoblast)[1:ncol(human_myoblast)]  = 'Human'
colnames(human_gm_gse68103)[1:ncol(human_gm_gse68103)]  =  'Human'
colnames(human_gm)[1:ncol(human_gm)]  = 'Human'
colnames(mouse_bcell)[1:ncol(mouse_bcell)]  = 'Mouse'
colnames(human_hsc)[1:ncol(human_hsc)]  = 'Human'
colnames(mouse_hsc)[1:ncol(mouse_hsc)]  = 'Mouse'


combined = cbind(human_neuron,mouse_neuron,human_myoblast,human_gm_gse68103,human_gm,mouse_bcell,human_hsc,mouse_hsc)

study = colnames(combined)
study = as.factor(study)

rownames(X) = 1:128

mint.plsda.res = mint.plsda(X = X, Y = Y, study = study, ncomp = 2)
#mint.plsda.res # lists the different functions
color = c('#999999','#dede00','#f781bf','#984ea3','#ff7f00','#e41a1c','#a65628','#4daf4a')


plotIndiv(mint.plsda.res, legend = TRUE,col=color, title = '', subtitle = '')



