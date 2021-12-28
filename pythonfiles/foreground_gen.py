import numpy as np
import pandas as pd
from collections import OrderedDict
import subprocess

class foreground_gen():
	def __init__(self):
		pass

	def reg_gen_dist(self, chr_start , chr_end , ref_start , ref_end, strand):
		rs = chr_start + 1
		re = chr_end
		gs = ref_start + 1
		ge = ref_end
		if strand == '+':
			l = [gs - rs , gs - re]
			l.sort()
			return l[0]
		else:
			l = [ge - rs , ge - re]
			l.sort()
			return l[0]

	def is_intronic(self,chr_start , chr_end , ref_start , ref_end):
		rs = chr_start + 1
		re = chr_end
		gs = ref_start + 1
		ge = ref_end
		if (rs > ge or gs>re):
			return 1
		else:
			return 0

	def foreground_calc(self, d_chr , d ,nearest_gene):
		larg_dist = 1e10
		for i in d_chr.keys():
			for val in d_chr[i]:
				cdL = - larg_dist   # left
				cdR =  larg_dist    # right
				cdI = larg_dist     # intronic
				cgL, cgR, cgI = "","",""
				for val2 in d[i]:
					if(self.is_intronic(val['Start'] , val['End'] , val2['txStart'] , val2['txEnd'])):
						dist = abs(self.reg_gen_dist(val['Start'] , val['End'] , val2['txStart'] , val2['txEnd'], val2['strand']))
						if (dist < cdI):
							cdI = dist
							cgI = val2['name2']
					else:
						dist = self.reg_gen_dist(val['Start'] , val['End'] , val2['txStart'] , val2['txEnd'], val2['strand'])
						if (dist <= 0):
							if (dist > cdL):
								cdL = dist
								cgL = val2['name2']
						elif (dist < cdR):
								cdR = dist
								cgR = val2['name2']
				cdist = []
				cdist.extend((abs(cdL) , cdR , cdI))
				cgene = []
				cgene.extend((cgL , cgR , cgI))
				cd_idx = np.argsort(cdist)
				cd = np.array(cdist)[cd_idx[:2]]
				cg = np.array(cgene)[cd_idx[:2]]
				if( cg[1] != ''):
					tmp_g = cgene[cd_idx[2]]
					cg[1] == cg[0] and (tmp_g != '') and cg[1] == tmp_g and cd[1] == cdist[cd_idx[2]]
				else:
					cg[1] = cg[0]
					cd[1] = cd[0]
				nearest_gene.iloc[val['index'],0] = cg[0]
				nearest_gene.iloc[val['index'],1] = cd[0]
				nearest_gene.iloc[val['index'],2] = cg[1]
				nearest_gene.iloc[val['index'],3] = cd[1]
		return nearest_gene
	
	def remove_zeropeaks(self,epi,chr_file):
		chr = pd.read_csv(chr_file, sep ='\t', header = None)
		ind = np.where(~np.asarray(epi).any(axis=1))[0]
		ind1 = set(range(epi.shape[0])) - set(ind)
		ind1 = sorted(ind1)
		epi1 = np.take(np.array(epi), list(ind1), axis = 0)
		chr1 = np.take(np.array(chr) , list(ind1) , axis =0)
		pd.DataFrame(chr1).to_csv(chr_file,sep= "\t", header = False , index = False)
		return epi1

	def nearest_gene_accurate(self, query_type, chr_file, acc_fast, query_file):
		epi = np.loadtxt(query_file,delimiter=",")
		#epi = self.remove_zeropeaks(epi,chr_file)
		
		if (query_type == 1 or query_type==3):
			ref = pd.read_csv('./storage/scepisearch/human/refseq-hg19.txt' , sep = '\t')
			ref.loc[:,'chrom'] = (ref['chrom'].str.split("_", expand=True)).iloc[: , 0]
			chr = pd.read_csv(chr_file, sep ='\t', header = None)

			d = OrderedDict()
			for i in ref['chrom'].unique():
				d[i] = [{'name' : ref['name'][j] , 'strand' : ref['strand'][j] , 'txStart' : ref['txStart'][j] ,'txEnd' : ref['txEnd'][j] ,'exonCount' : ref['exonCount'][j] ,'name2': ref['name2'][j]} for j in ref[ref['chrom']==i].index]

			cmd = ['Rscript ./storage/scepisearch/human/accessibility_score_faster/global_score.R '+chr_file+' ./acc_score.csv ./foreground.csv']
			process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		else:
			ref = pd.read_csv('./storage/scepisearch/mouse/refGene.txt' , sep = '\t')
			ref.loc[:,'chrom'] = (ref['chrom'].str.split("_", expand=True)).iloc[: , 0]
			chr = pd.read_csv(chr_file, sep ='\t', header = None)

			d = OrderedDict()
			for i in ref['chrom'].unique():
				d[i] = [{'name' : ref['name'][j] , 'strand' : ref['strand'][j] , 'txStart' : ref['txStart'][j] ,'txEnd' : ref['txEnd'][j] ,'exonCount' : ref['exonCount'][j] ,'name2': ref['name2'][j]} for j in ref[ref['chrom']==i].index]

			cmd = ['Rscript ./storage/scepisearch/mouse/accessibility_score_faster/global_score.R '+chr_file+' ./acc_score.csv ./foreground.csv']
			process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		#print(cmd)
		#print(process)
		if(process == 0):
			nearest_gene_epi_path = './foreground.csv'
			nearest_gene = pd.read_csv(nearest_gene_epi_path, sep="\t", header=None)
			if (acc_fast==2):
				acc_score = np.loadtxt("./acc_score.csv")
				epi = epi[(nearest_gene!=0).all(axis=1)]
				overlap = epi.shape[0]/nearest_gene.shape[0]
				overlap = np.array([overlap])
				acc_score = acc_score[(nearest_gene!=0).all(axis=1)]
				nearest_gene = nearest_gene[(nearest_gene!=0).all(axis=1)]
				np.savetxt("./acc_score.csv",acc_score)
				np.savetxt("./overlap.txt",overlap,fmt='%f')
				return nearest_gene, epi
			else:
				ind_bool = (nearest_gene!=0).all(axis=1)
				ind_zero = np.where(ind_bool)[0]
				overlap = len(ind_zero)/nearest_gene.shape[0]
				overlap = np.array([overlap])
				d_chr = OrderedDict()
				for i in chr.iloc[:,0].unique():
					d_chr[i] = [{'Start' : chr.iloc[j,1] ,'End' : chr.iloc[j,2] ,'index' : j} for j in chr[chr.iloc[:,0]==i].index if j not in ind_zero]
				nearest_gene = self.foreground_calc(d_chr , d , nearest_gene)
				nearest_gene.to_csv("./foreground.csv",header=False,sep="\t")
				np.savetxt("./overlap.txt",overlap,fmt='%f')
				return nearest_gene, epi
		else:
			nearest_gene = pd.DataFrame({})
			return nearest_gene,epi
