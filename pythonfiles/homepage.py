from __future__ import print_function
from operator import itemgetter
import itertools, subprocess, sys, os, csv, random, glob, re, gzip, time, matplotlib,  multiprocessing , sklearn.preprocessing
import numpy as np , pandas as pd , matplotlib.pyplot as plt , numpy_indexed as npi , seaborn as sns , SimpSOM as sps, matplotlib.pyplot as plt
from sklearn import cluster
from scipy.stats import hypergeom
matplotlib.use('Agg')
plt.switch_backend('agg')
from collections import OrderedDict
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle, gc, linecache
from test_final_v2 import test_final_v2
from pdf_report_generator import generate_pdf
import sys, uuid, math, random
from FITSPhase1 import FITSPhase1
from FITSPhase2 import FITSPhase2
from itertools import repeat , cycle, islice

# for Python3
from tkinter import *
from tkinter.ttk import *
from tkinter.filedialog import askopenfilename
from tkinter import messagebox

from resultpage import ResultPage
from resultpage_old import ResultPage_old
from resultpage_emb import ResultPage_emb


##global parameters
class HomePage(Frame):
	"""
		GUI Application for Transcriptome Search Standalone Application
	"""
	def __init__(self, master):
		"""
			Constructor
		"""
		self.master = master
		Frame.__init__(self, self.master)
		self.grid()
		# logo
		img = PhotoImage(file="logo.png")
		label = Label(self.master, image=img)
		label.image = img # keep a reference!
		label.grid(row=0, column=0, columnspan=3, padx=50, pady=10)
		# Notebook Tabs
		self.note = Notebook(self.master)
		self.homeTab = Frame(self.note)
		self.ImpTab = Frame(self.note)
		self.EmbTab = Frame(self.note)
		self.contactUsTab = Frame(self.note)
		self.note.add(self.homeTab, text = "SEARCH FORM")
		self.note.add(self.ImpTab, text = "IMPUTATION")
		self.note.add(self.EmbTab, text = "EMBEDDING")
		self.note.add(self.contactUsTab, text = "CONTACT US")
		self.note.grid(row=1, column=0, columnspan=3, pady=10, padx=50)
		self.master.protocol("WM_DELETE_WINDOW", self.close_window)

	def close_window(self):
		try:
			if tkMessageBox.askokcancel("Caution", "Are you sure you want to Exit?", parent=self.master):
				print ("Quitting .. ")
				self.master.destroy()
		except Exception as e:
			if messagebox.askokcancel("Caution", "Are you sure you want to Exit? ", parent=self.master):
				print ("Quitting .. ")
				self.master.destroy()

	def openfile(self, choice):
		if choice == 1:
			self.query_file = askopenfilename(initialdir=os.getcwd(),
			title="Select Query file",
			filetypes=[("csv Files", "*.txt")])
		if choice == 2:
			self.chr_file = askopenfilename(initialdir=os.getcwd(),
			title="Select Chromosome Peak file",
			filetypes=[("bed files", "*.bed")])

	def openfile_emb(self, choice, query_file, chr_file,x):
		if choice == 1:
			#self.query_file[-1]= askopenfilename(initialdir=os.getcwd(),title="Select Query file",filetypes=[("csv Files", "*.txt")])
			self.query_file.append(askopenfilename(initialdir=os.getcwd(),title="Select Query file",filetypes=[("csv Files", "*.txt")]))

		if choice == 2:
			#self.chr_file[-1] = askopenfilename(initialdir=os.getcwd(),title="Select Chromosome Peak file",filetypes=[("bed files", "*.bed")])
			self.chr_file.append(askopenfilename(initialdir=os.getcwd(),title="Select Chromosome Peak file",filetypes=[("bed files", "*.bed")]))


	def verify_arguments(self, query_cells):
		"""TODO : Verify each argument """
		pass

	def output_segmentation(self,query_cells):
		print ("Appending phenotype and p values in results .. ")
		if int(self.var.get()) == 1:
			metadata_epi = './meta_human/metadata_epi_new.csv'
			metadata_exp = './meta_human/metadata_exp_with_tissue.csv'
		elif int(self.var.get()) == 2:
			metadata_epi = './meta_mouse/metadata_epi.csv'
			metadata_exp = './meta_mouse/metadata_exp.csv'
		else:
			metadata_exp = './meta_mouse/metadata_exp.csv'
		if int(self.var.get()) == 1 or int(self.var.get()) == 2:
			output_file_exp = open('./exp.txt')
			output_file_epi = open('./epi.txt')
			file_pval_exp = open('./pval_exp.txt')
			file_pval_epi = open('./pval_epi.txt')
			file_fdr_epi = open('./fdr_epi.txt')
			file_fdr_exp = open('./fdr_exp.txt')
			count = 0
			for line1,line2,line3,line4,line5,line6 in zip(output_file_exp,file_pval_exp,file_fdr_exp,output_file_epi,file_pval_epi,file_fdr_epi):
				matrix1 = list()
				matrix2 = list()
				line1 = line1.split(' ')
				line2 = line2.split(' ')
				line3 = line3.split(' ')
				line4 = line4.split(' ')
				line5 = line5.split(' ')
				line6 = line6.split(' ')
				line1[-1] = line1[-1].strip()
				line2[-1] = line2[-1].strip()
				line3[-1] = line3[-1].strip()
				line4[-1] = line4[-1].strip()
				line5[-1] = line5[-1].strip()
				line6[-1] = line6[-1].strip()
				temp_filename_exp = "./storage/query_exp_" + str(count) + ".txt"
				temp_filename_epi = "./storage/query_epi_" + str(count) + ".txt"
				with open(temp_filename_exp, "w") as f1 , open(temp_filename_epi, "w") as f2:
					for i,j,k,l,m,n in zip(line1, line2, line3, line4, line5, line6):
						line = linecache.getline(metadata_exp, int(i)+1)
						line = line.strip()
						if int(self.var.get()) == 1:
							index = line.find('@')
							line = line[index+1:]
						line = line+"@"+j+"@"+k
						f1.write(line+'\n')
						line_epi = linecache.getline(metadata_epi, int(l)+1)
						line_epi = line_epi.strip()
						line_epi = line_epi+"@"+m+"@"+n
						f2.write(line_epi+'\n')
				count = count + 1
		else:
			output_file_exp = open('./exp.txt')
			file_pval_exp = open('./pval_exp.txt')
			file_fdr_exp = open('./fdr_exp.txt')
			count = 0
			for line1,line2,line3 in zip(output_file_exp,file_pval_exp,file_fdr_exp):
				matrix1 = list()
				matrix2 = list()
				line1 = line1.split(' ')
				line2 = line2.split(' ')
				line3 = line3.split(' ')
				line1[-1] = line1[-1].strip()
				line2[-1] = line2[-1].strip()
				line3[-1] = line3[-1].strip()
				temp_filename_exp = "./storage/query_exp_" + str(count) + ".txt"
				temp_filename_epi = "./storage/query_epi_" + str(count) + ".txt"
				with open(temp_filename_exp, "w") as f1:
					for i,j,k in zip(line1, line2, line3):
						line = linecache.getline(metadata_exp, int(i)+1)
						line = line.strip()
						line = line+"@"+j+"@"+k
						f1.write(line+'\n')
				count = count + 1
		filename = generate_pdf(query_cells, int(self.var.get()))
		print ("All Query Result files processed Successfully .. ")
		print ("")

	def display_results(self, query_cells, results_per_study, query_type):
		# self.master = Tk()
		# self.master.geometry("%dx%d+0+0" % (1366, 768))
		res_page = Toplevel(self.master)
		res_page.title("Results")
		# res_page.resizable(0,0)
		w, h = res_page.winfo_screenwidth(), res_page.winfo_screenheight()
		res_page.geometry("%dx%d+0+0" % (1700, 900))
		res_page.resizable(0,0)
		if query_type == 1 or query_type == 2:
			res_page_app = ResultPage(res_page, query_cells, results_per_study,query_type)
		else:
			res_page_app = ResultPage_old(res_page, query_cells, results_per_study, query_type)
		res_page.mainloop()
		exit()

	def display_results_emb(self, query_cells, results_per_study, query_type,sample_each_data):
		# self.master = Tk()
		# self.master.geometry("%dx%d+0+0" % (1366, 768))
		res_page = Toplevel(self.master)
		res_page.title("Results")
		# res_page.resizable(0,0)
		w, h = res_page.winfo_screenwidth(), res_page.winfo_screenheight()
		res_page.geometry("%dx%d+0+0" % (1500, 900))
		res_page.resizable(0,0)
		res_page_app = ResultPage_emb(res_page, query_cells, results_per_study, query_type,sample_each_data)
		res_page.mainloop()
		exit()


	def demo(self):
		self.chr_file = "chromosome_human.bed"
		self.query_file = "./k562.txt"		
		self.execute()
		

	def execute(self):
		form_type='search'
		"""Main command to execute the script"""
		print("Input wrapper check....")
		cmd = ['python3 ./pythonfiles/wrapper.py '+self.chr_file+' '+self.query_file]
		process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		process = 0
		if process == 0:
			print ("EXECUTION STARTED")
			print ("")
			f_query = open(self.query_file)
			reader = list(csv.reader(f_query,delimiter=","))
			query_cells= len(reader[0])
			top_study = int(self.rps.get())
			cluster = int(self.cls.get())
			#print(int(self.anno_unanno.get()))

			#query_cells= 10
			#top_study = 5
			#cluster = 5
			
			#print(self.chr_file)
			#print(self.query_file)
			# cmd = ['python3 ./pythonfiles/main.py '+self.chr_file+' '+self.query_file+' '+str(self.rps.get())+' '+str(self.var.get())+' '+str(self.var1.get())]
			# process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			test_final_v2(self.chr_file,self.query_file,int(self.rps.get()),int(self.var.get()),int(self.var1.get()),int(self.active_poised.get()),int(self.cls.get()),form_type,int(self.anno_unanno.get()))
			process=0
			if process ==0:
				print("Cell search complete!!!")
				self.output_segmentation(query_cells)
				print ("")
				print ("EXECUTION COMPLETE")
				self.display_results(query_cells,int(self.rps.get()), int(self.var.get()))
			else:
				print("Cell search error")
		else:
			print("INPUT WRAPPER ERROR!!!")

	def execute_background(self):
		form_type='search'
		"""Main command to execute the script"""
		print("Input wrapper check....")
		cmd = ['python3 ./pythonfiles/wrapper.py '+self.chr_file+' '+self.query_file]
		process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		process = 0
		if process == 0:
			print ("EXECUTION STARTED")
			print ("")
			f_query = open(self.query_file)
			reader = list(csv.reader(f_query,delimiter=","))
			query_cells= len(reader[0])
			top_study = int(self.rps.get())
			cluster = int(self.cls.get())
			
			#print(self.chr_file)
			#print(self.query_file)
			cmd = ['nohup python3 ./pythonfiles/main.py '+self.chr_file+' '+self.query_file+' '+str(self.rps.get())+' '+str(self.var.get())+' '+str(self.var1.get())+' '+str(self.active_poised.get())+' '+str(self.cls.get())+' '+str(form_type)+' '+str(self.anno_unanno.get())+' >log.out &']
			process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			pid = process.pid
			print("Please Note process ID of process to track:")
			print(pid)
			#test_final_v2(self.chr_file,self.query_file,int(self.rps.get()),int(self.var.get()),int(self.var1.get()),int(self.active_poised.get()),int(self.cls.get()),form_type,int(self.anno_unanno.get()))
			#print(cmd)
			if process ==0:
				print("Cell search complete!!!")
				self.output_segmentation(query_cells)
				print ("")
				print ("EXECUTION COMPLETE")
				self.display_results(query_cells,int(self.rps.get()), int(self.var.get()))
			else:
				print("Cell search error")
		else:
			print("INPUT WRAPPER ERROR!!!")


	def execute_imp(self):
		"""Main command to execute the script"""
		print ("EXECUTION STARTED")
		print ("")
		f_query = open(self.query_file)
		reader = list(csv.reader(f_query,delimiter=","))
		query_cells= len(reader[0])

		trees = int(self.rps_imp.get())
		topk = int(self.var2_imp.get())
		colwise = int(self.var1_imp.get())
		maxClusters = 8
		msc = 8
		aRFSR = [60,100]
		num = str(uuid.uuid1())
		name2save = 'FITS_OUTPUT'
		maxAllowedLevel = trees
		
		dataX1 = np.loadtxt(self.query_file, delimiter=",")
		dataX1 = np.transpose(dataX1)
		average = dataX1.mean(axis=1)
		dataX1 = dataX1/average[:, np.newaxis]+0.000001
		FITSPhase1(dataX1,maxClusters,msc,aRFSR,maxAllowedLevel,name2save,num)
		FITSPhase2(dataX1, topk, colwise, name2save)
		process=0
		if process ==0:
			print("Imputation complete!!!")
			print ("")
			print ("EXECUTION COMPLETE")
			messagebox.showinfo(title="EXECUTION COMPLETE", message="Final Imputed output exit in Project folder with name FINAL_OUTPUT_(feature/sample).txt!")
		else:
			print("Imputation error")


	def execute_form2(self):
		#EMBEDDING DATASET FORM - 2
		labelframe3 = LabelFrame(self.EmbTab, text="Embedding FORM")
		labelframe3.grid(row=2, column=2, pady = 5, padx=10)

		no_data = int(self.data.get())
		self.vars = []
		self.query_file = []
		self.chr_file = []
		for x in range(1, no_data+1):
			self.vars.append(IntVar())
			#self.query_file[x= StringVar()
			#self.chr_file[x] = StringVar()
			Button(labelframe3, text="Choose Count "+str(x), command=lambda:self.openfile_emb(1,self.query_file,self.chr_file,x)).grid(row=11+x, column=1, padx=5, pady=5, sticky=E)
			Button(labelframe3, text="Choose Peak "+str(x), command=lambda:self.openfile_emb(2,self.query_file,self.chr_file,x)).grid(row=11+x, column=2, padx=5, pady=5, sticky=E)
			R1 = Radiobutton(labelframe3, text="Human Dataset (hg19)", variable=self.vars[-1], value=1)
			R1.grid(row=11+x, column=3)
			R2 = Radiobutton(labelframe3, text="Mouse Dataset (mm9)", variable=self.vars[-1], value=2)
			R2.grid(row=11+x, column=4)
			x = x + 1
		Button(labelframe3, text="SUBMIT - FINAL", command=self.execute_emb).grid(row=11+x+1,column=2, padx=5, pady=5, sticky=E)

		
	def execute_emb(self):
		form_type='emb'
		#print(self.query_file)
		#print(self.chr_file)
		dataset = int(self.data.get())
		final_query_cells = 0
		final_species = 3
		exp = np.array([], dtype=np.int64).reshape(0,int(self.rps.get()))
		pval_exp = np.array([], dtype=np.int64).reshape(0,int(self.rps.get()))
		fdr_exp = np.array([], dtype=np.int64).reshape(0,int(self.rps.get()))
		marker_gene = np.array([], dtype=np.int64).reshape(0,50)
		celltype_markers = pd.DataFrame()
		top_clusters = np.array([], dtype=np.str).reshape(int(self.cls.get()),0)
		#pathway_scores = pd.DataFrame()
		enhancers = pd.DataFrame()
		overlap = ""	
		sample_each_data = []	
		count=1 
		for i in range(0,dataset):
			print ("EXECUTION STARTED")
			print ("")
			f_query = open(self.query_file[i])
			reader = list(csv.reader(f_query,delimiter=","))
			query_cells= len(reader[0])
			sample_each_data.extend(list(repeat(count,query_cells)))
			final_query_cells = final_query_cells + query_cells
			top_study = int(self.rps.get())
			cluster = int(self.cls.get())
			if int(self.vars[i].get()) == 1:
				species = 3
			else:
				species = 2
			#test_final_v2(self.chr_file[i],self.query_file[i],int(self.rps.get()),species,int(self.var1.get()),int(self.active_poised.get()),int(self.cls.get()),form_type,int(self.anno_unanno_emb.get()))
			#cmd = ['Rscript ./pythonfiles/unipath.R '+self.chr_file[i]+' '+self.query_file[i]+' '+str(self.vars[i].get())+'']
			#print(cmd)
			#process = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			#print(process)
			#if process == 0:
			#	pathway_scores_1 = pd.read_csv('./pathway_scores.txt',sep=' ',index_col=0)
			exp_1 = np.loadtxt('./exp.txt', delimiter=" ")
			pval_exp_1 = np.loadtxt('./pval_exp.txt', delimiter=" ")
			fdr_exp_1 = np.loadtxt('./fdr_exp.txt', delimiter=" ")
			marker_gene_1 = np.loadtxt('./marker_gene.txt',delimiter=" ",dtype='str')
			celltype_markers_1 = pd.read_csv('./celltype_markers.txt',sep='\t',header=None)
			top_clusters_1 = np.loadtxt('./top_clusters.txt', delimiter=" ", dtype='str')
			enhancers_1 = pd.read_csv('./enhancers.txt',sep='\t',header=None)
			overlap_1 = np.loadtxt("./overlap.txt",dtype='str')

			exp = np.vstack([exp,exp_1])
			pval_exp = np.vstack([pval_exp,pval_exp_1])	
			fdr_exp = np.vstack([fdr_exp,fdr_exp_1])	
			marker_gene = np.vstack([marker_gene,marker_gene_1])
			top_clusters = np.hstack([top_clusters,top_clusters_1])
			enhancers = pd.concat([enhancers,enhancers_1])
			celltype_markers = pd.concat([celltype_markers,celltype_markers_1])
			#pathway_scores.index = pathway_scores_1.index
			#pathway_scores = pd.concat([pathway_scores_1,pathway_scores_1])
			overlap = overlap+","+str(overlap_1)	
			count=count+1	
		np.savetxt('./top_clusters.txt', top_clusters , delimiter=" ", fmt='%s')
		np.savetxt('./marker_gene.txt', marker_gene , delimiter=" ", fmt = '%s')
		celltype_markers.to_csv('./celltype_markers.txt',sep="\t",header=False,index=False)		
		np.savetxt('./fdr_exp.txt', fdr_exp , delimiter=" ", fmt='%f')
		np.savetxt('./exp.txt', exp, fmt='%d', delimiter=" ")
		np.savetxt('./pval_exp.txt', pval_exp, delimiter=" ", fmt='%f')	
		overlap = np.array([overlap])	
		np.savetxt("./overlap.txt",overlap,fmt='%s')
		enhancers.to_csv("./enhancers.txt",sep=" ",header=None, index = None)
		#pathway_scores.to_csv("./pathway_scores.txt",sep=" ",header=None)

		#

		# #GRAPH
		# top_clusters = np.loadtxt('./top_clusters.txt',delimiter=" ",dtype='str')
		# #adj = np.zeros((top_clusters.shape[1], top_clusters.shape[1]))
		# clus = top_clusters.ravel()
		# res = [t.replace('_', '.') for t in clus] 
		# res = np.array(res).reshape(top_clusters.shape[0],top_clusters.shape[1])
		# adj = kneighbors_graph(np.transpose(res), 5, mode='distance', include_self=True)
		# adj = adj.toarray()
		# adj = np.nan_to_num(adj)
		# adj = np.abs(adj)                                   # no negative weights
		# adj = adj - np.diag(np.diag(adj))                   # no self loops
		# adj = np.triu(adj) + np.transpose(np.triu(adj))
		# # for i in range(top_clusters.shape[1]):
    			# # for j in range(i+1, top_clusters.shape[1]):
        			# # adj[i, j] = len(np.intersect1d(top_clusters[:,i],top_clusters[:,j]))
		# clustering = SpectralClustering(n_clusters=2,assign_labels="discretize",random_state=0,affinity='precomputed').fit(adj)
		# ann = []
		# for i in range(final_query_cells):
    			# ann.append("D"+str(sample_each_data[i])+" "+"Q"+str(i))
		# colors = np.array(list(islice(cycle(['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']),int(max(clustering.labels_) + 1))))
		# # add black color for outliers (if any)
		# colors = np.append(colors, ["#000000"])
		# node_names = ann
		# g = igraph.Graph.Adjacency((adj > 5).astype(bool).tolist(),mode="UNDIRECTED")
		# # Add edge weights and node labels.
		# g.es['weight'] = adj[adj.nonzero()]
		# g.vs['label'] = node_names  # or a.index/a.columns
		# color=colors[clustering.labels_]
		# g.vs["color"] = color
		# igraph.plot(g,"./graph.png",vertex_size=40,vertex_shape="circle",vertex_label_size=9)
		
		print("Cell search complete!!!")
		self.output_segmentation_embedding(final_query_cells)
		print ("")
		print ("EXECUTION COMPLETE")
		self.display_results_emb(final_query_cells,int(self.rps.get()), final_species,sample_each_data)
		
	def output_segmentation_embedding(self,query_cells):
		print ("Appending phenotype and p values in results .. ")
		metadata_exp = './meta_mouse/metadata_exp.csv'
		output_file_exp = open('./exp.txt')
		file_pval_exp = open('./pval_exp.txt')
		file_fdr_exp = open('./fdr_exp.txt')
		count = 0
		for line1,line2,line3 in zip(output_file_exp,file_pval_exp,file_fdr_exp):
			matrix1 = list()
			matrix2 = list()
			line1 = line1.split(' ')
			line2 = line2.split(' ')
			line3 = line3.split(' ')
			line1[-1] = line1[-1].strip()
			line2[-1] = line2[-1].strip()
			line3[-1] = line3[-1].strip()
			temp_filename_exp = "./storage/query_exp_" + str(count) + ".txt"
			temp_filename_epi = "./storage/query_epi_" + str(count) + ".txt"
			with open(temp_filename_exp, "w") as f1:
				for i,j,k in zip(line1, line2, line3):
					line = linecache.getline(metadata_exp, int(i)+1)
					line = line.strip()
					line = line+"@"+j+"@"+k
					f1.write(line+'\n')
			count = count + 1
		filename = generate_pdf(query_cells, 3)
		print ("All Query Result files processed Successfully .. ")
		print ("")
		

	def AddContents(self):
		#### HOME TAB
		## ABOUT Section
		labelframe1 = LabelFrame(self.homeTab, text="ABOUT")
		labelframe1.grid(row=1, column=0, columnspan = 2, pady=5, padx=10)
		Label(labelframe1, text=pickle.load(open('pythonfiles/about.pickle', 'rb')), \
			wraplength = 400).grid(row=0, column=0)
		labelframe5 = LabelFrame(self.homeTab, text="REFERENCE DATASET")
		labelframe5.grid(row=2, column=0, columnspan = 1, pady=5, padx=5)
		Label(labelframe5, text=pickle.load(open('pythonfiles/about_new.pickle', 'rb')), \
			wraplength = 800).grid(row=2, column=0)
		## SEARCH Section
		labelframe2 = LabelFrame(self.homeTab, text="SEARCH FORM")
		labelframe2.grid(row=1, column=2, pady = 5, padx=10)
		Label(labelframe2, text="Search against Reference Data", \
			font="underline").grid(row=0, column=0, columnspan=3, sticky=W, pady=10)

		self.var = IntVar()
		R1 = Radiobutton(labelframe2, text="Human Dataset (hg19)", variable=self.var, value=1)
		R1.grid(row=1, column=0)
		R2 = Radiobutton(labelframe2, text="Mouse Dataset (mm9)", variable=self.var, value=2)	## No command added
		R2.grid(row=1, column=1)
		R3 = Radiobutton(labelframe2, text="Cross-Species (Human (hg19) -> Mouse)", variable=self.var, value=3)	## No command added
		R3.grid(row=1, column=2)

		Label(labelframe2, text="Select Active/Poised State", \
			font="underline").grid(row=2, column=0, columnspan=3, sticky=W, pady=10)
		self.active_poised = IntVar()
		R1 = Radiobutton(labelframe2, text="Active Genes", variable=self.active_poised, value=1)
		R1.grid(row=3, column=0)
		R2 = Radiobutton(labelframe2, text="Poised Genes", variable=self.active_poised, value=2)	## No command added
		R2.grid(row=3, column=1)

		Label(labelframe2, text="Upload Query (count) File").grid(row=7, column=0, sticky=W)
		Button(labelframe2, text="Choose File .. ", command=lambda: \
			self.openfile(1)).grid(row=7, column=1, padx=5, pady=5, sticky=E)
		Label(labelframe2, text="Upload Chromosome Peak File").grid(row=8, column=0, sticky=W)
		Button(labelframe2, text="Choose File .. ", command=lambda: \
			self.openfile(2)).grid(row=8, column=1, padx=5, pady=5, sticky=E)
		Label(labelframe2, text="Results per study").grid(row=9, column=0, sticky=W)
		self.rps = Entry(labelframe2)

		Label(labelframe2, text="No of Top clusters to search").grid(row=10, column=0, sticky=W)
		self.cls = Entry(labelframe2)
		self.cls.grid(row=10, column=1, padx=5, pady=5)
		self.cls.delete(0, END)
		self.cls.insert(5, "20")

		self.var1 = IntVar()
		R3 = Radiobutton(labelframe2, text="Accurate", variable=self.var1, value=1)
		R3.grid(row=11, column=0)
		R4 = Radiobutton(labelframe2, text="Faster", variable=self.var1, value=2)	## No command added
		R4.grid(row=11, column=1)
		self.rps.grid(row=9, column=1, padx=5, pady=5)
		self.rps.delete(0, END)
		self.rps.insert(0, "5")

		
		self.anno_unanno = IntVar()
		Label(labelframe2, text="Select Annotated/Both Annotated & Unannotated Ref. Cells", \
			font="underline").grid(row=12, column=0, columnspan=3, sticky=W, pady=10)
		R1 = Radiobutton(labelframe2, text="Annotated Cells", variable=self.anno_unanno, value=1)
		R1.grid(row=13, column=0)
		R2 = Radiobutton(labelframe2, text="Both Annotated & Unannotated Ref. Cells", variable=self.anno_unanno, value=2)	## No command added
		R2.grid(row=13, column=1)
		
		Button(labelframe2, text="DEMO", command=self.demo).grid(row=14,
			column=0, padx=5, pady=5, sticky=W)
		Button(labelframe2, text="SUBMIT", command=self.execute).grid(row=14,
			column=1, padx=5, pady=5, sticky=E)
		Button(labelframe2, text="SUBMIT IN BACKGROUND", command=self.execute_background).grid(row=14,
			column=2, padx=5, pady=5, sticky=E)
		self.var.set(1)
		self.active_poised.set(2)
		self.var1.set(2)
		self.anno_unanno.set(2)
		# self.update_form()

		#### IMPUTATION TAB
		## ABOUT Section
		labelframe1 = LabelFrame(self.ImpTab, text="ABOUT")
		labelframe1.grid(row=1, column=0, columnspan = 2, pady=5, padx=10)
		Label(labelframe1, text=pickle.load(open('pythonfiles/about_imp.pickle', 'rb')), \
			wraplength = 400).grid(row=0, column=0)
		## SEARCH Section
		labelframe2 = LabelFrame(self.ImpTab, text="Imputation FORM")
		labelframe2.grid(row=1, column=2, pady = 5, padx=10)

		Label(labelframe2, text="Upload Data (count) File").grid(row=7, column=0, sticky=W)
		Button(labelframe2, text="Choose File .. ", command=lambda: \
			self.openfile(1)).grid(row=7, column=1, padx=5, pady=5, sticky=E)

		Label(labelframe2, text="Maximum Trees").grid(row=8, column=0, sticky=W)
		self.rps_imp = Entry(labelframe2)

		self.rps_imp.grid(row=8, column=1, padx=5, pady=5)
		self.rps_imp.delete(0, END)
		self.rps_imp.insert(0, "4")

		Label(labelframe2, text="Phase 2(FeatureWise/SampleWise)").grid(row=9, column=0, sticky=W)
		self.var1_imp = Entry(labelframe2)

		self.var1_imp.grid(row=9, column=1, padx=5, pady=5)
		self.var1_imp.delete(0, END)
		self.var1_imp.insert(0, "1")

		Label(labelframe2, text="Phase 2(correlated matrix feature/sample value)").grid(row=10, column=0, sticky=W)
		self.var2_imp = Entry(labelframe2)

		self.var2_imp.grid(row=10, column=1, padx=5, pady=5)
		self.var2_imp.delete(0, END)
		self.var2_imp.insert(0, "1")

		Button(labelframe2, text="SUBMIT", command=self.execute_imp).grid(row=12,
			column=2, padx=5, pady=5, sticky=E)

		
		#### EMBEDDING TAB
		## ABOUT Section
		labelframe1 = LabelFrame(self.EmbTab, text="ABOUT")
		labelframe1.grid(row=1, column=0, columnspan = 2, pady=5, padx=10)
		Label(labelframe1, text=pickle.load(open('pythonfiles/about_emb.pickle', 'rb')), \
			wraplength = 400).grid(row=0, column=0)
		## SEARCH Section
		labelframe2 = LabelFrame(self.EmbTab, text="Embedding FORM")
		labelframe2.grid(row=1, column=2, pady = 5, padx=10)

		Label(labelframe2, text="No of Dataset").grid(row=2, column=0, sticky=W)
		self.data = Entry(labelframe2)
		self.data.grid(row=2, column=1, padx=5, pady=5)
		self.data.delete(0, END)
		self.data.insert(0, "4")

		
		Label(labelframe2, text="Select Active/Poised State", \
			font="underline").grid(row=3, column=0, columnspan=3, sticky=W, pady=10)
		self.active_poised = IntVar()
		R1 = Radiobutton(labelframe2, text="Active Genes", variable=self.active_poised, value=1)
		R1.grid(row=4, column=0)
		R2 = Radiobutton(labelframe2, text="Poised Genes", variable=self.active_poised, value=2)	## No command added
		R2.grid(row=4, column=1)

		Label(labelframe2, text="Results per study").grid(row=5, column=0, sticky=W)
		self.rps = Entry(labelframe2)
		self.rps.grid(row=5, column=1, padx=5, pady=5)
		self.rps.delete(0, END)
		self.rps.insert(0, "5")

		
		Label(labelframe2, text="No of Top clusters to search").grid(row=6, column=0, sticky=W)
		self.cls = Entry(labelframe2)
		self.cls.grid(row=6, column=1, padx=5, pady=5)
		self.cls.delete(0, END)
		self.cls.insert(0, "20")

		
		R3 = Radiobutton(labelframe2, text="Accurate", variable=self.var1, value=1)
		R3.grid(row=7, column=0)
		R4 = Radiobutton(labelframe2, text="Faster", variable=self.var1, value=2)	## No command added
		R4.grid(row=7, column=1)

		Label(labelframe2, text="Select Annotated/Both Annotated & Unannotated Ref. Cells", \
			font="underline").grid(row=8, column=0, columnspan=3, sticky=W, pady=10)
		self.anno_unanno_emb = IntVar()
		R1 = Radiobutton(labelframe2, text="Annotated Cells", variable=self.anno_unanno_emb, value=1)
		R1.grid(row=9, column=0)
		R2 = Radiobutton(labelframe2, text="Both Annotated & Unannotated Ref. Cells", variable=self.anno_unanno_emb, value=2)	## No command added
		R2.grid(row=9, column=1)


		Button(labelframe2, text="SUBMIT - 1", command=self.execute_form2).grid(row=14,
			column=2, padx=5, pady=5, sticky=E)



		#### COONTACT US TAB
		## Contact Us Section
		labelframe3 = LabelFrame(self.contactUsTab, text="CONTACT US")
		labelframe3.grid(row=1, column=0, padx=10, pady=30)
		Label(labelframe3, text=pickle.load(open('pythonfiles/guide_info.pickle', 'rb')), \
			font=("Helvetica", 12)).grid(row=0, column=0)
		## Developed By Section
		labelframe4 = LabelFrame(self.contactUsTab, text="DEVELOPED BY")
		labelframe4.grid(row=1, column=1, padx=10, pady=30)
		Label(labelframe4, text=pickle.load(open('pythonfiles/author_info.pickle', 'rb')), \
			font=("Helvetica", 12)).grid(row=0, column=0)
