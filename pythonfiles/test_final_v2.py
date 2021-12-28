#!/usr/bin/python
import numpy as np , pandas as pd
from foreground_gen import foreground_gen
from exp_match import exp_match
from epi_match import epi_match
from gene_enrichment import gene_enrichment
from multiprocessing import Process
import numpy_indexed as npi
from itertools import cycle, islice
# from exp_match_corr import exp_match_corr

class test_final_v2():
	def __init__(self,chr_file,query_file,top_study,query_type,acc_fast,active_poised,cluster,form_type,anno_unanno):
		self.exp_obj = exp_match()
		self.epi_obj = epi_match()
		# self.exp_obj_corr = exp_match_corr()
		self.gene_enriched_obj = gene_enrichment()
		self.nearest_gene_obj = foreground_gen()
		self.process_query(chr_file, query_file, top_study, query_type, acc_fast,active_poised,cluster,form_type,anno_unanno)

	def get_digit(self,x):
		return(int(re.search(r'\d+', x).group()))

	def process_query(self, chr_file, epi_path, top_study, query_type, acc_fast,active_poised,cluster,form_type,anno_unanno):
				if(query_type==1 or query_type == 3):
					sc_gene_path = './storage/scepisearch/human/genes_21159.txt'
				else:
					sc_gene_path = './storage/scepisearch/mouse/gene_mouse.csv'

				gene_list = list()
				with open(sc_gene_path, 'r') as fl:
					for l in fl.readlines():
						gene_list.append(l.rstrip('\n'))

				nearest_gene,epi = self.nearest_gene_obj.nearest_gene_accurate(query_type,chr_file,acc_fast,epi_path)
				if not nearest_gene.empty:
					acc_score_path = './acc_score.csv'
					acc_score = np.loadtxt(acc_score_path)
					# epi = np.loadtxt(epi_path, delimiter=",")
					acc_score[acc_score == 0] = 1
					acc_fast=2

					epi = np.array(epi)/(acc_score[:, np.newaxis])

					start_nearest = list(nearest_gene.ix[:,1])
					end_nearest = list(nearest_gene.ix[:,3])
					genename_start = list(nearest_gene.ix[:,0])
					genename_end = list(nearest_gene.ix[:,2])

					#gene names from nearest genes from either of the two gene locations
					epi_gene = [0] * len(start_nearest)
					for k in range(len(start_nearest)):
						if start_nearest[k] <= end_nearest[k]:
							epi_gene[k] = genename_start[k]
						else:
							epi_gene[k] = genename_end[k]

					gene_enriched = self.gene_enriched_obj.gene_enrichment_calc(epi,gene_list,nearest_gene,query_type)
					if int(active_poised)==1:
						N = 500
					else:
						N = 1000
					
					if query_type == 1 or query_type == 2:
						sorted_col_idx_epi = np.argsort(gene_enriched, axis=0)[gene_enriched.shape[0]-N::,:]
						top_epi_gene_epi = pd.DataFrame(sorted_col_idx_epi).apply(lambda x: np.take(gene_list, indices = x) , axis = 0)
						gene_name_exp_loc_epi = top_epi_gene_epi.apply(lambda x : npi.indices(gene_list ,x, missing = -1)  , axis = 0)

						sorted_col_idx_exp = np.argsort(epi, axis=0)[epi.shape[0]-N::,:]
						top_epi_gene_exp = pd.DataFrame(sorted_col_idx_exp).apply(lambda x: np.take(epi_gene, indices = x) , axis = 0)
						
						enhancers = top_epi_gene_exp
						enhancers = enhancers.transpose()
						enhancers.to_csv("./enhancers.txt",sep=" ",header=None, index = None)
						
						gene_name_exp_loc_exp = top_epi_gene_exp.apply(lambda x : npi.indices(gene_list ,x, missing = -1)  , axis = 0)

						# self.epi_obj.epi_matching(epi,query_type,gene_list,nearest_gene,top_study,epi_gene,gene_enriched,gene_name_exp_loc)
						# self.exp_obj.exp_matching(epi, gene_enriched , query_type, top_study, gene_name_exp_loc)
						
						if form_type=='search':
							self.epi_obj.epi_matching(epi,query_type,gene_list,nearest_gene,top_study,epi_gene,gene_enriched,gene_name_exp_loc_epi)
							self.exp_obj.exp_matching(epi,gene_enriched,query_type,top_study,gene_name_exp_loc_exp,cluster,anno_unanno)
							#p1 = Process(target=self.epi_obj.epi_matching, args=(epi,query_type,gene_list,nearest_gene,top_study,epi_gene,gene_enriched,gene_name_exp_loc_epi))
							#p2 = Process(target=self.exp_obj.exp_matching, args=(epi,gene_enriched,query_type,top_study,gene_name_exp_loc_exp,cluster,anno_unanno))
							#p1.start()
							#p2.start()
							#p1.join()
							#p2.join()
						else:
							self.exp_obj.exp_matching(epi,gene_enriched,query_type,top_study,gene_name_exp_loc_exp,cluster,anno_unanno)

						print('cell search done')
					else:
						gene_mouse = list()
						with open("./storage/scepisearch/mouse/gene_mouse.csv", 'r') as fl:
							for l in fl.readlines():
								gene_mouse.append(l.rstrip('\n'))
						human_gene_lower = [x.lower() for x in gene_list]
						mouse_gene_lower = [x.lower() for x in gene_mouse]
						intersect_gene = list(set(human_gene_lower) & set(mouse_gene_lower))
						gene_enriched_mouse = np.zeros((23595,epi.shape[1])) 
						
						for i in range(len(human_gene_lower)):
							if human_gene_lower[i] in mouse_gene_lower:
								ind = mouse_gene_lower.index(human_gene_lower[i])
								gene_enriched_mouse[ind,:] = gene_enriched[i,:]
							else:
								continue
						
						epi_gene = [x.lower() for x in epi_gene]
						sorted_col_idx_exp = np.argsort(epi, axis=0)[epi.shape[0]-1000::,:]
						top_epi_gene_exp = pd.DataFrame(sorted_col_idx_exp).apply(lambda x: np.take(epi_gene, indices = x) , axis = 0)
						gene_name_exp_loc = top_epi_gene_exp.apply(lambda x : npi.indices(mouse_gene_lower ,x, missing = -1)  , axis = 0)						

						enhancers = top_epi_gene_exp
						enhancers = enhancers.transpose()
						enhancers.to_csv("./enhancers.txt",sep=" ",header=None, index = None)


						self.exp_obj.exp_matching(epi,gene_enriched_mouse,query_type,top_study,gene_name_exp_loc,cluster,anno_unanno)
						self.epi_obj.epi_matching(epi,query_type,mouse_gene_lower,nearest_gene,top_study,epi_gene,gene_enriched_mouse,gene_name_exp_loc)
						print('cell search done')

