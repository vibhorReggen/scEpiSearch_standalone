import multiprocessing, itertools, pickle
from scipy.stats import hypergeom
import numpy as np
import pandas as pd
import numpy_indexed as npi
import matplotlib.pyplot as plt

class gene_enrichment():
	def __init__(self):
		pass

	def apply_par(self, args):
		d, nearest_gene, gene_list = args
		m_gene = d.apply(lambda x : self.gene_enrich(x , nearest_gene, gene_list), axis = 0)
		return m_gene

	def hypergeometric_test(self, foreground , background):
		fgc1 = list(foreground.iloc[:,0])
		fgd1 = list(foreground.iloc[:,1])
		fgc2 = list(foreground.iloc[:,2])
		fgd2 = list(foreground.iloc[:,3])

		for i in range(len(fgd1)):
			if fgd1[i] <= fgd2[i]:
				fgd2[i] = 2000000
			else:
				fgd1[i] = 2000000

		afc = np.append(fgc1, fgc2)
		afd = np.append(fgd1,fgd2)
		pos = np.where(afd < 1000000)
		afc = afc[pos]

		#---------------------------------
		bgc1 = list(background.iloc[:,0])
		bgd1 = list(background.iloc[:,1])
		bgc2 = list(background.iloc[:,2])
		bgd2 = list(background.iloc[:,3])


		for i in range(len(bgd1)):
			if bgd1[i] <= bgd2[i]:
				bgd2[i] = 2000000
			else:
				bgd1[i] = 2000000

		abc = np.append(bgc1, bgc2)
		abd = np.append(bgd1,bgd2)
		pos = np.where(abd < 1000000)
		abc = abc[pos]

		#---------------------------------
		N = int(len(abc)) + 1
		n = int(len(afc)) + 1

		fgGenes =  set(afc)
		fgGenes = list(fgGenes)
		apvals = list()

		for i in range(len(fgGenes)):
			pos = np.where(afc == fgGenes[i])[0]
			pos1 = np.where(abc == fgGenes[i])[0]
			k = len(pos)
			K = len(pos1)
			mn = min(K,n) +1
			pval = 0

			for j in range(k,mn):
				pval = pval + hypergeom.pmf(j, N, n, K)
			apvals.append(pval)

		summary = np.hstack((np.expand_dims(apvals,axis=1),np.expand_dims(fgGenes, axis=1)))
		summaryDF = pd.DataFrame(summary, columns=["P-value", "Genes"])
		return summaryDF

	def gene_enrich(self, x, nearest_gene, gene_list):
		col_ind = int(x.name)
		fore = np.take(nearest_gene , x , axis = 0)
		gene_enrich = self.hypergeometric_test(fore, nearest_gene)
		gene_pval = list(gene_enrich['Genes'])
		pval_total = gene_enrich['P-value'].astype(float)
		pval_total = list(pval_total)
		pval_total.append(1)
		pval_total = -np.log2(pval_total)
		index = pval_total.argsort()[-50:][::-1]
		mark_gene = np.asarray(gene_pval)[index]
		mark_gene_enrichment = np.asarray(pval_total)[index]
		ind = npi.indices(gene_pval ,gene_list, missing = -1)
		pval_corr = np.asarray(pval_total)[ind]
		return pd.Series([mark_gene, pval_corr,mark_gene_enrichment])

	def gene_enrichment_calc(self,epi,gene_list,nearest_gene,query_type):
		gene_enriched = np.zeros([epi.shape[1], len(gene_list)])
		marker_gene = np.empty([epi.shape[1], len(range(50))] ,dtype = 'object')

		Z =10000
		sorted_col_idx = np.argsort(epi, axis=0)[epi.shape[0]-Z::,:]

		workers = 20
		pool = multiprocessing.Pool(processes=workers)
		split_dfs = np.array_split(pd.DataFrame(sorted_col_idx), workers, axis=1)

		res = pd.concat(pool.map(self.apply_par, [(d ,nearest_gene, gene_list) for d in split_dfs]), axis = 1)
		pool.close()
		res_0 = list(itertools.chain.from_iterable(res.iloc[0]))
		res_1 = list(itertools.chain.from_iterable(res.iloc[1]))
		res_2 = list(itertools.chain.from_iterable(res.iloc[2]))
		marker_gene = pd.DataFrame(np.array(res_0).reshape(sorted_col_idx.shape[1],50))
		marker_gene_enrichment = np.array(np.array(res_2).reshape(sorted_col_idx.shape[1], 50))		

		if query_type == 1 or query_type == 3:
			a_file = open("./meta_human/markers_human.pkl", "rb")
			dict_markers = pickle.load(a_file)
			a_file.close()
		else:
			a_file = open("./meta_mouse/markers_mouse.pkl", "rb")
			dict_markers = pickle.load(a_file)
			a_file.close()
		celltype_markers = pd.DataFrame(index = range(marker_gene.shape[0]), columns = range(marker_gene.shape[1]))
        
		for i in range(marker_gene.shape[0]):
			marker_gene.iloc[i,:] = marker_gene.iloc[i,:].str.upper()
			res = marker_gene.iloc[i,:].map(dict_markers)
			celltype_markers.iloc[i,:] = res


		gene_enriched = pd.DataFrame(np.array(res_1).reshape(sorted_col_idx.shape[1], len(gene_list)))
		gene_enriched = np.transpose(np.array(gene_enriched))

		gene_enriched = pd.DataFrame(gene_enriched)
		gene_enriched.index = gene_list

		gene_enriched.to_csv('./enrichment_scores.txt', sep=" ")
		np.savetxt('./marker_gene.txt', marker_gene , delimiter=" ", fmt = '%s')
		celltype_markers.to_csv('./celltype_markers.txt',sep="\t",header=False,index=False)
		np.savetxt('./marker_gene_enrichment.txt',marker_gene_enrichment, delimiter=" ",fmt='%f')

		m_c1 = pd.Series(np.array(marker_gene).ravel()).dropna().value_counts()
		plt.figure(figsize=(10,5))
		plt.bar(m_c1.iloc[:60].index,m_c1.iloc[:60].values, width=0.7,color='g', edgecolor='black')
		plt.xticks(rotation=90,fontsize=8)
		plt.xlim([0,60])
		plt.savefig('./gene_frequency.png')
		return np.array(gene_enriched)
