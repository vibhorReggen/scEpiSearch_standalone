import pandas as pd
import numpy as np
import sklearn.preprocessing, multiprocessing, gzip
#sgd = keras.optimizers.SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')
import numpy_indexed as npi
from collections import OrderedDict
class exp_match():
	def __init__(self):
		pass

	def exp_matching(self, epi, gene_enriched, query_type, top_study, gene_name_exp_loc,cluster,anno_unanno):
		#autoencoder output
		#reduced_gene_enriched = self.autoencoder_output(query_type, gene_enriched)
		# print(reduced_gene_enriched.shape)

		#find relevant clusters
		res = self.clusters_corr_exp(query_type,gene_enriched)
		# print(res.shape)
	
		if query_type==1:
			index_value = np.genfromtxt("./storage/scepisearch/human/exp_new/clusters_celltypewise/subclusters/mean_array_subclusterindex.txt",dtype='str')
			#print(index_value)
		else:
			index_value = np.genfromtxt("./storage/scepisearch/mouse/exp_new/clusters_celltype/mean_array_subclusterindex.txt",dtype='str')
    
		#N = 10
		N = int(cluster)
		if anno_unanno==1 and query_type==1:
			N = 30
      
		sorted_clust = np.argsort(res, axis=0)[res.shape[0]-N::,:]
		#     print(sorted_clust)
		sorted_clust = index_value[sorted_clust.ravel()].reshape((sorted_clust.shape[0],sorted_clust.shape[1]))
		print(sorted_clust.shape)

		un = np.unique(sorted_clust)
		id_clust = OrderedDict()
		for i,j in enumerate(un):
			id_clust[str(j)]=[]
			id_clust[str(j)][:] = np.where(sorted_clust == j)[1]

		workers= 20
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

		p.map(self.cluster_exp_par, [(val,key,final_res,lock,query_type,gene_name_exp_loc,anno_unanno) for val,key in zip(values,keys)])
		p.close()

		final_corr = np.zeros([epi.shape[1],top_study], dtype='int')
		pval_epi = np.zeros([epi.shape[1],top_study])
		final_fdr = np.zeros([epi.shape[1],top_study])

		for i in range(epi.shape[1]):
			ind = np.argsort(final_res[i,'corr'])[::-1][:top_study]
			final_corr[i,:] = final_res[i,'index'][ind].astype(int)
			pval_epi[i,:] = final_res[i,'pval'][ind]
			final_fdr[i,:] = final_res[i,'adj_pval'][ind]

		np.savetxt('./fdr_exp.txt', final_fdr , delimiter=" ", fmt='%f')
		np.savetxt('./exp.txt', final_corr, fmt='%d', delimiter=" ")
		np.savetxt('./top_clusters.txt', sorted_clust, fmt='%s', delimiter=" ")
		np.savetxt('./pval_exp.txt', pval_epi, delimiter=" ", fmt='%f')
		#np.savetxt('./median.txt',final_corr,delimiter=" ",fmt='%f')

	def autoencoder_output(self, query_type, gene_enriched):
		gene_enriched = pd.DataFrame(np.transpose(gene_enriched))
		if(query_type==1):
			json_file = open('./storage/scepisearch/human/exp/model_human.json', 'r')
			loaded_model_json = json_file.read()
			json_file.close()
			loaded_model = model_from_json(loaded_model_json)
			loaded_model.load_weights("./storage/scepisearch/human/exp/model_human.h5")
		else:
			json_file = open('./storage/scepisearch/mouse/exp/model_mouse.json', 'r')
			loaded_model_json = json_file.read()
			json_file.close()
			loaded_model = model_from_json(loaded_model_json)
			loaded_model.load_weights("./storage/scepisearch/mouse/exp/model_mouse.h5")
		loaded_model.compile(loss = 'mean_squared_error', optimizer=sgd)
		encoder = Model(inputs=loaded_model.input, outputs=loaded_model.get_layer("dense_5").output)

		highest_non_inf = gene_enriched.max().loc[lambda v: v<np.Inf].max()
		gene_enriched_new = gene_enriched.replace(np.Inf, highest_non_inf + 100)
		gene_enriched_new = sklearn.preprocessing.scale(gene_enriched_new, axis = 0)
		reduced_gene_enriched = encoder.predict(gene_enriched_new)

		return np.transpose(np.array(reduced_gene_enriched))

	def clusters_corr_exp(self,query_type,gene_enriched):
		if query_type==1:
			mean_array = np.load("./storage/scepisearch/human/exp_new/clusters_celltypewise/subclusters/mean_array.npy")
		else:
			mean_array = np.load("./storage/scepisearch/mouse/exp_new/clusters_celltype/mean_array.npy")
			
			
		where_are_NaNs = np.isnan(mean_array)
		mean_array[where_are_NaNs] = 0
    
		where_are_NaNs = np.isnan(gene_enriched)
		gene_enriched[where_are_NaNs] = 0


		nr, nc = gene_enriched.shape
		xvec = ro.FloatVector(gene_enriched.transpose().reshape((gene_enriched.size)))
		xr = ro.r.matrix(xvec, nrow=nr, ncol=nc)

		nry, ncy = mean_array.shape
		xvecy = ro.FloatVector(mean_array.transpose().reshape((mean_array.size)))
		yr = ro.r.matrix(xvecy, nrow=nry, ncol=ncy)

		res = ro.r.cor(xr,yr, method="spearman")
		res = np.array(res)
		res = np.transpose(res)
		return res

	def median_calc(self, x , exp_ref):
		ind_top = [x for x in x if x != -1]
		mat_top = np.take(exp_ref , ind_top , axis = 0)
		med_top = np.median(mat_top , axis=0)
		return med_top
		
	def median_calc_null(self, x , corr_mat, exp_ref):
		#ind_top = [x for x in x if x != -1]
		mat_top = np.take(exp_ref , x , axis = 0)
		med_top = np.median(mat_top , axis=0)
		#print(med_top.shape)
		corr_mat[x.name,:] = med_top

	def cluster_exp_par(self, args):
		val,key,final_res,lock,query_type,gene_name_exp_loc,anno_unanno = args
		if(query_type == 1):
			f0 = gzip.GzipFile('./storage/scepisearch/human/exp_new/clusters_celltypewise/subclusters/clust_'+str(key)+'.npy.gz', "r")
			exp_ref = np.load(f0)
			#f3 = gzip.GzipFile('/home/cellsearch/cellatlassearch_shreya/epigenome_search/human/exp/null_model/clust_'+str(key)+'.npy.gz', "r")
			#corr_mat = np.load(f3)
			sorted_col = np.load("./storage/scepisearch/human/exp_new/clusters_celltypewise/null_idx_enrichment.npy")
			clust_infor = pd.read_csv("./storage/scepisearch/human/exp_new/clusters_celltypewise/subclusters/clusters_final.txt",sep=" ",dtype='str')
			unanno = np.loadtxt("./meta_human/exp_unknown_human.txt", delimiter=",")
		else:
			f0 = gzip.GzipFile('./storage/scepisearch/mouse/exp_new/clusters_celltype/subclusters_twentymillion/clust_'+str(key)+'.npy.gz', "r")
			exp_ref = np.load(f0)
			#f3 = gzip.GzipFile('./searchProject/storage/scepisearch/mouse_new_exp/exp/null_model/clust_'+str(key)+'.npy.gz', "r")
			#corr_mat = np.load(f3)
			sorted_col = np.loadtxt("./storage/scepisearch/mouse/exp_new/clusters_celltype/gene_name_exp_loc_mouse.txt",delimiter=",")
			sorted_col=sorted_col.astype(int)
			clust_infor = pd.read_csv("./storage/scepisearch/mouse/exp_new/clusters_celltype/clusters_final.txt",sep=" ",dtype='str')
			unanno = np.loadtxt("./meta_mouse/exp_unknown_mouse.txt", delimiter=",")
		ind_ref = np.array(clust_infor.loc[clust_infor['cluster_no'] == str(key), 'id'])
		if anno_unanno==1:
			exp_ref = exp_ref[:,~np.isin(ind_ref.astype(int),unanno.astype(int))]
			ind_ref = ind_ref[~np.isin(ind_ref.astype(int),unanno.astype(int))]

		if exp_ref.shape[1]==0:
			return
		
		where_are_NaNs = np.isnan(exp_ref)
		exp_ref[where_are_NaNs] = 0

		query = np.take(gene_name_exp_loc, list(val), axis=1)
		res = query.apply(lambda x : self.median_calc(x , exp_ref) ,axis = 0)
		res = np.array(res)
		if query_type == 1:
			corr_mat = np.zeros([21159,exp_ref.shape[1]])
		else:
			corr_mat = np.zeros([23595,exp_ref.shape[1]])
		pd.DataFrame(sorted_col).apply(lambda x : self.median_calc_null(x , corr_mat ,exp_ref) ,axis = 0)
		corr_mat = np.array(corr_mat)

		#score_ranking = np.argsort(corr_mat, axis=1)[:,corr_mat.shape[1]-50::] 
		#top_ref = np.unique(score_ranking)
		#res = res[top_ref,:]

		pval = []
		sorted_10_idx = np.argsort(res, axis=0)[res.shape[0]-10::,:]
		sorted_raveled = sorted_10_idx.ravel()
		col_idx = np.arange(res.shape[1])
		val_top = res[sorted_10_idx, col_idx]
		val_top_raveled = val_top.ravel()
		for a, b in zip(sorted_raveled, val_top_raveled):
			pval.append(sum(corr_mat[:,a] >= b))
		pval = np.reshape(pval , (sorted_10_idx.shape[0] , sorted_10_idx.shape[1]))

		pval = pval/float(1000)
		pval = pval + 0.00001

		score_ranking = np.argsort(pval, axis=1)[:,pval.shape[1]-50::] 
		top_ref = np.unique(score_ranking)
		res = res[top_ref,:]

		#pval = np.random.uniform(low=0.01, high=0.05, size=(pval.shape[0],pval.shape[1]))

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
