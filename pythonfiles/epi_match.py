from gene_enrichment import gene_enrichment
import sklearn.preprocessing, multiprocessing, gzip, matplotlib
import matplotlib.pyplot as plt  , seaborn as sns , SimpSOM as sps, matplotlib.pyplot as plt
from collections import OrderedDict
from sklearn import cluster
#from keras.models import model_from_json
#from keras.models import Model
#sgd = keras.optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')
import numpy as np
import pandas as pd
import numpy_indexed as npi
matplotlib.use('Agg')
plt.switch_backend('agg')
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

class epi_match():
	def __init__(self):
		self.gene_enriched_obj = gene_enrichment()
		pass

	def rand_jitter(self, arr):
		stdev = .1*(max(arr)-min(arr))
		return arr + np.random.randn(len(arr)) * stdev

	def epi_matching(self, epi, query_type, gene_list, nearest_gene, top_study, epi_gene, gene_enriched, gene_name_exp_loc):
		pinf = float('+inf')
		ninf = float('-inf')
		gene_enriched[gene_enriched == ninf] = 0.0
		gene_enriched[gene_enriched == pinf] = np.finfo(np.float16).max
        
		where_are_NaNs = np.isnan(gene_enriched)
		gene_enriched[where_are_NaNs] = 0

		if(query_type == 1):
			correlation_matrix = np.zeros([epi.shape[1],742297])
		else:
			correlation_matrix = np.zeros([epi.shape[1],81173])

		# gene_enriched = self.gene_enriched_obj.gene_enrichment_calc(epi,gene_list,nearest_gene)
		# print(np.sum(gene_enriched))
		#
		# N = 2000
		# sorted_col_idx_epi = np.argsort(epi, axis=0)[epi.shape[0]-N::,:]
		# top_epi_gene = pd.DataFrame(sorted_col_idx_epi).apply(lambda x: np.take(epi_gene, indices = x) , axis = 0)
		# gene_name_exp_loc = top_epi_gene.apply(lambda x : npi.indices(gene_list ,x, missing = -1)  , axis = 0)

		#autoencoder output
		#reduced_gene_enriched = self.autoencoder_output(query_type, gene_enriched)
		#print(reduced_gene_enriched.shape)

		#find relevant clusters
		# res = self.clusters_corr_epi(query_type,reduced_gene_enriched)
		# print(res.shape)

		res = self.clusters_corr_epi(query_type,gene_enriched)

		# N = 10
		sorted_clust = np.argsort(res, axis=0)[res.shape[0]-20::,:]
		
		if query_type==1:
			index_value = np.genfromtxt("./storage/scepisearch/human/epi_new/mean_array_subclusterindex.txt",dtype='str')
			print(index_value)
			sorted_clust = index_value[sorted_clust.ravel()].reshape((sorted_clust.shape[0],sorted_clust.shape[1]))
			print(sorted_clust.shape)

		un = np.unique(sorted_clust)
		id_clust = OrderedDict()
		for i,j in enumerate(un):
			id_clust[str(j)]=[]
			id_clust[str(j)][:] = np.where(sorted_clust == j)[1]

		workers= 10
		p = multiprocessing.Pool(processes=workers)
		keys, values= zip(*id_clust.items())

		manager = multiprocessing.Manager()
		final_res = manager.dict()
		lock = manager.Lock()

		for i in range(epi.shape[1]):
			final_res[i,'corr'] = np.array([])
			final_res[i,'index'] = np.array([])
			final_res[i,'pval'] = np.array([])
			final_res[i,'adj_pval'] = np.array([])

		p.map(self.cluster_epi_par, [(val,key,final_res,lock,query_type,gene_enriched) for val,key in zip(values,keys)])
		p.close()

		final_corr = np.zeros([epi.shape[1],top_study], dtype='int')
		pval_epi = np.zeros([epi.shape[1],top_study])
		final_fdr = np.zeros([epi.shape[1],top_study])

		for i in range(epi.shape[1]):
			ind = np.argsort(final_res[i,'corr'])[::-1][:top_study]
			correlation_matrix[i,final_res[i,'index'].astype(int)] = final_res[i,'corr']
			final_corr[i,:] = final_res[i,'index'][ind].astype(int)
			pval_epi[i,:] = final_res[i,'pval'][ind]
			final_fdr[i,:] = final_res[i,'adj_pval'][ind]
		
		ind_corrmat_nonzero = ~np.all(correlation_matrix == 0, axis=0)
		ind_corrmat_nonzero = np.where(ind_corrmat_nonzero)[0]
		correlation_matrix = correlation_matrix[:,~np.all(correlation_matrix == 0, axis=0)]

		pinf = float('+inf')
		ninf = float('-inf')
		correlation_matrix[correlation_matrix == ninf] = 0.0
		correlation_matrix[correlation_matrix == pinf] = np.finfo(np.float16).max

		where_are_NaNs = np.isnan(correlation_matrix)
		correlation_matrix[where_are_NaNs] = 0

		np.savetxt('./fdr_epi.txt', final_fdr , delimiter=" ", fmt='%f')
		np.savetxt('./epi.txt', final_corr, fmt='%d', delimiter=" ")
		np.savetxt('./pval_epi.txt', pval_epi, delimiter=" ", fmt='%f')
		#np.savetxt('./autoencoder.txt', reduced_gene_enriched , delimiter=" ", fmt='%f')
		np.savetxt('./enrichment_scores.txt', gene_enriched , delimiter=" ", fmt='%f')

		net_size = int(np.ceil(np.sqrt(5*np.sqrt(epi.shape[1]))))
		net = sps.somNet(net_size, net_size, correlation_matrix, PBC=False)
		net.train(0.01, 1000)
		bmuList1 = net.project(correlation_matrix, show=False, printout=False)
		clust_infor = self.cluster_new(net, correlation_matrix, bmuList1, type= 'DBSCAN',cutoff=1, min_samples=1)
		clust_infor.ix[:,3] = pd.Series(self.rand_jitter(np.array(clust_infor.ix[:,3])))
		clust_infor.ix[:,4] = pd.Series(self.rand_jitter(np.array(clust_infor.ix[:,4])))
		clust_infor.to_csv('./tsne.txt',sep =" ", index = False , header = False)

		correlation_matrix = pd.DataFrame(correlation_matrix)
		if query_type == 1:
			meta_exp = pd.read_csv("./meta_human/metadata_epi_new.csv", header=None,sep="@")
		else:
			meta_exp = pd.read_csv("./meta_mouse/metadata_epi.csv", header=None,sep="@")

		meta_exp = meta_exp.iloc[:,0]
		meta_exp = np.array(meta_exp)
		meta_exp = meta_exp[ind_corrmat_nonzero]
		correlation_matrix.columns = meta_exp
		correlation_matrix.to_csv('./correlation_matrix.txt' , sep=" ")
		ax = sns.clustermap(correlation_matrix, figsize=(8, 4))
		if correlation_matrix.shape[1]>20:
			ax.ax_heatmap.tick_params(axis='x', which='both', length=0,labelsize=0)
		ax.savefig('./heatmap.png')

		return gene_enriched, gene_name_exp_loc
	
	def clusters_corr_null(self,null_model,exp_ref):
		nr, nc = exp_ref.shape
		xvec = ro.FloatVector(exp_ref.transpose().reshape((exp_ref.size)))
		xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)

		nry, ncy = null_model.shape
		xvecy = ro.FloatVector(null_model.transpose().reshape((null_model.size)))
		yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)

		res = ro.r.cor(xr,yr, method="spearman")
		res = np.array(res)
		res = np.transpose(res)
		return res


	def median_calc_null(self, x , corr_mat, exp_ref):
		#ind_top = [x for x in x if x != -1]
		mat_top = np.take(exp_ref , x , axis = 0)
		med_top = np.median(mat_top , axis=0)
		#print(med_top.shape)
		corr_mat[x.name,:] = med_top

	def cluster_epi_par(self, args):
		val,key,final_res,lock,query_type,gene_enriched = args
		if query_type==1:
			f0 = gzip.GzipFile('./storage/scepisearch/human/epi_new/clusters_allgene_new_subclusters/clust_'+str(key)+'.npy.gz', "r")
			exp_ref = np.load(f0)
			null_model = np.load('./storage/scepisearch/human/epi_new/null_model.npy')
			#null_model = np.load(f3)
			clust_infor = pd.read_csv("./storage/scepisearch/human/epi_new/clusters_subclusters.txt",sep=" ")
		else:
			f0 = gzip.GzipFile('./storage/scepisearch/mouse/epi/clusters_allgene/clust_'+str(key)+'.npy.gz', "r")
			exp_ref = np.load(f0)
			f3 = gzip.GzipFile('./storage/scepisearch/mouse/epi/null_model/clust_'+str(key)+'.npy.gz', "r")
			corr_mat = np.load(f3)
			clust_infor = pd.read_csv("./storage/scepisearch/mouse/epi/clusters.txt",sep=" ")

		if query_type == 1:
			pinf = float('+inf')
			ninf = float('-inf')
			null_model[null_model == ninf] = 0.0
			null_model[null_model == pinf] = np.finfo(np.float16).max
			where_are_NaNs = np.isnan(null_model)
			null_model[where_are_NaNs] = 0
        
			exp_ref[exp_ref == ninf] = 0.0
			exp_ref[exp_ref == pinf] = np.finfo(np.float16).max
			where_are_NaNs = np.isnan(exp_ref)
			exp_ref[where_are_NaNs] = 0
			corr_mat = self.clusters_corr_null(null_model,exp_ref)		
			ind_ref = list(clust_infor.loc[clust_infor['cluster_no'] == str(key), 'id'])
		else:
			ind_ref = list(clust_infor.loc[clust_infor['cluster_no'] == int(key), 'id'])

		#gene_name_exp_loc_query = np.take(gene_name_exp_loc, list(val), axis=1)
		#
		#res = gene_name_exp_loc_query.apply(lambda x : self.median_calc(x , exp_ref) ,axis = 0)
		#res = np.array(res)

		query = np.take(gene_enriched, list(val), axis=1)

		nr, nc = query.shape
		xvec = ro.FloatVector(query.transpose().reshape((query.size)))
		xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)

		nry, ncy = exp_ref.shape
		xvecy = ro.FloatVector(exp_ref.transpose().reshape((exp_ref.size)))
		yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)

		res = ro.r.cor(xr,yr, method="spearman")
		res = np.array(res)
		res = np.transpose(res)

		pval = []
		sorted_10_idx = np.argsort(res, axis=0)[res.shape[0]-50::,:]
		sorted_raveled = sorted_10_idx.ravel()
		col_idx = np.arange(res.shape[1])
		val_top = res[sorted_10_idx, col_idx]
		val_top_raveled = val_top.ravel()
		for a, b in zip(sorted_raveled, val_top_raveled):
			pval.append(sum(corr_mat[:,a] >= b))
		pval = np.reshape(pval , (sorted_10_idx.shape[0] , sorted_10_idx.shape[1]))

		pval = pval/float(1000)
		pval = pval + 0.00001

		#pd.DataFrame(sorted_col).apply(lambda x : self.median_calc_null(x , corr_mat ,exp_ref) ,axis = 0)
		#corr_mat = np.array(corr_mat)

		#score_ranking = np.argsort(pval, axis=1)[:,pval.shape[1]-50::] 
		#top_ref = np.unique(score_ranking)
		#res = res[top_ref,:]

		p_adjust_epi = stats.p_adjust(FloatVector(pval.ravel()), method = 'BH')
		p_adjust_epi = np.array(p_adjust_epi).reshape(pval.shape[0], pval.shape[1])

		result = np.array(ind_ref)[sorted_10_idx.ravel()]
		result = result.reshape(sorted_10_idx.shape[0],sorted_10_idx.shape[1])

		with lock:
			for i,j in enumerate(val):
				final_res[j,'corr'] = np.append(final_res[j,'corr'], val_top[:,i])
				final_res[j,'index'] = np.append(final_res[j,'index'], result[:,i])
				final_res[j,'pval'] = np.append(final_res[j,'pval'], pval[:,i])
				final_res[j,'adj_pval'] = np.append(final_res[j,'adj_pval'],p_adjust_epi[:,i])

	def median_calc(self, x , exp_ref):
		ind_top = [x for x in x if x != -1]
		mat_top = np.take(exp_ref , ind_top , axis = 0)
		med_top = np.median(mat_top , axis=0)
		return med_top

	def clusters_corr_epi(self,query_type,reduced_gene_enriched):
		if query_type==1:
			mean_array_reduced = np.load("./storage/scepisearch/human/epi_new/mean_array.npy")
		else:
			mean_array_reduced = np.load("./storage/scepisearch/mouse/epi/mean_array.npy")
		nr, nc = reduced_gene_enriched.shape
		xvec = ro.FloatVector(reduced_gene_enriched.transpose().reshape((reduced_gene_enriched.size)))
		xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)

		nry, ncy = mean_array_reduced.shape
		xvecy = ro.FloatVector(mean_array_reduced.transpose().reshape((mean_array_reduced.size)))
		yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)

		res = ro.r.cor(xr,yr, method="spearman")
		res = np.array(res)
		res = np.transpose(res)
		return res

	def autoencoder_output(self, query_type,gene_enriched):
		gene_enriched = pd.DataFrame(np.transpose(gene_enriched))
		if(query_type==1):
			json_file = open('./storage/scepisearch/human/epi/model_human.json', 'r')
			loaded_model_json = json_file.read()
			json_file.close()
			loaded_model = model_from_json(loaded_model_json)
			loaded_model.load_weights("./storage/scepisearch/human/epi/model_human.h5")
		else:
			json_file = open('./storage/scepisearch/mouse/epi/model_mouse.json', 'r')
			loaded_model_json = json_file.read()
			json_file.close()
			loaded_model = model_from_json(loaded_model_json)
			loaded_model.load_weights("./storage/scepisearch/mouse/epi/model_mouse.h5")
		loaded_model.compile(loss = 'mean_squared_error', optimizer=sgd)
		encoder = Model(inputs=loaded_model.input, outputs=loaded_model.get_layer("dense_5").output)

		highest_non_inf = gene_enriched.max().loc[lambda v: v<np.Inf].max()
		gene_enriched_new = gene_enriched.replace(np.Inf, highest_non_inf + 100)
		gene_enriched_new = sklearn.preprocessing.scale(gene_enriched_new, axis = 0)
		reduced_gene_enriched = encoder.predict(gene_enriched_new)

		return np.transpose(np.array(reduced_gene_enriched))

	def cluster_new(self, net, array, bmuList, type='DBSCAN', cutoff=1, min_samples=5):
		cluster_info = pd.DataFrame(index=range(array.shape[0]), columns=['id','color','cluster_no','xc','yc'])

		if type=='DBSCAN':
			cl = cluster.DBSCAN(eps=cutoff, min_samples=min_samples).fit(bmuList)

		j=0
		randCl = []
		for k in range(len(np.unique(cl.labels_))):
			randCl.append("#%06x" % np.random.randint(0, 0xFFFFFF))
		for i in range(len(cl.labels_)):
			cluster_info.iloc[j,0] = i
			cluster_info.iloc[j, 1] = randCl[cl.labels_[i]]
			cluster_info.iloc[j,2] = cl.labels_[i]
			cluster_info.iloc[j,3] = bmuList[i][0]
			cluster_info.iloc[j,4] = net.netHeight-bmuList[i][1]
			j = j+1

		return cluster_info
