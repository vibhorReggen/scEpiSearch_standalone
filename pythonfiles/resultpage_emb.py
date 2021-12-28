## results page
import matplotlib , sys, subprocess, csv, webbrowser, os, pickle, re
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
from matplotlib.figure import Figure
import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import networkx as nx
from itertools import repeat , cycle, islice
from sklearn.cluster import SpectralClustering
from sklearn.neighbors import kneighbors_graph
from pdf_report_generator import generate_pdf

try:
    # for Python2
    from Tkinter import *
    from ttk import *
    import tkMessageBox
except ImportError:
    # for Python3
    from tkinter import *
    from tkinter.ttk import *
    from tkinter import messagebox



class ResultPage_emb(Frame):
	def __init__(self, master, query_cells, results_per_study, query_type,sample_each_data):
		self.master = master
		self.results_per_study = results_per_study
		self.query_cells = query_cells
		self.query_type = query_type
		self.sample_each_data = sample_each_data
		Frame.__init__(self, self.master)
		self.grid()
		# logo
		img = PhotoImage(file="logo.png")
		imglabel = Label(self.master, image=img)
		imglabel.image = img # keep a reference!
		imglabel.grid(row=0, column=1, columnspan=3, sticky='news', padx=50, pady=10)
		self.labelframe1 = LabelFrame(self.master, text="Query Cells")
		self.labelframe1.grid(row=1, column=0, pady=5, padx=20)
		self.labelframe2 = LabelFrame(self.master, text="Phenotype")
		self.labelframe2.grid(row=1, column=1, columnspan=2, pady=5, padx=20)
		#self.labelframe3 = LabelFrame(self.master, text="Phenotype")
		#self.labelframe3.grid(row=1, column=5, columnspan=2, pady=5, padx=20)
		self.addcontents_frame1(query_cells)
		self.addcontents_frame2()
		#self.addcontents_frame3()
		self.update_information_exp('0') ## display cell 0's information
		#self.update_information_epi('0')
		#Button(self.master, text="Close", command=self.close_window).grid(\
		#	row=2, column=0)
		Button(self.master, text="View Report", command=self.generate_report).grid(\
			row=2, column=1, sticky=E, padx=5)
		Button(self.master, text="Show WordCloud", command=self.display_wordcloud).grid(\
			row=2, column=2, sticky=W, padx=5)
		Button(self.master, text="Show Enhancers", command=self.display_enhancers).grid(\
			row=3, column=1, sticky=W, padx=5)
		Button(self.master, text="Show Graph", command=self.display_graph).grid(\
			row=3, column=2, sticky=W, padx=5)

		#self.master.protocol("WM_DELETE_WINDOW", self.close_window_main)

	def generate_report(self):
		filename = generate_pdf(self.query_cells,self.query_type)
		#try:
		#	tkMessageBox.showinfo("PDF Generated", "File : " + filename , parent=self.master)
		#except Exception as e:
		#	messagebox.showinfo("PDF Generateddd", "File : " + filename , parent=self.master)
		subprocess.Popen(['evince', filename])

	def close_window_main(self):
		try:
			if tkMessageBox.askokcancel("Caution", "Exit? \
				(Make sure you have downloaded the report)", parent=self.master):
				self.clear_files()
				self.master.destroy()
		except Exception as e:
			if messagebox.askokcancel("Caution", "Exit? \
				(Make sure you have downloaded the report)", parent=self.master):
				self.clear_files()
				self.master.destroy()

	def close_window(self):
		try:
			if tkMessageBox.askokcancel("Caution", "Exit?", parent=self.master):
				self.clear_files()
				self.master.destroy()
		except Exception as e:
			if messagebox.askokcancel("Caution", "Exit?", parent=self.master):
				self.clear_files()
				self.master.destroy()


	def clear_files(self):
		pass
		# os.remove('./acc_score.csv')
		# os.remove('./exp.txt')
		# os.remove('./epi.txt')
		# os.remove('./pval_exp.txt')
		# os.remove('./pval_epi.txt')
		# os.remove('./fdr_exp.txt')
		# os.remove('./fdr_epi.txt')
		# os.remove('./foreground.csv')
		# os.remove('./tsne.txt')
		# for i in range(int(self.query_cells)):
		# 	os.remove('data/query_exp_' + str(i) + '.txt')
		# 	os.remove('data/query_epi_' + str(i) + '.txt')
		# os.remove('./heatmap.png')

	def exec_1(self,i):
		self.update_information_exp(i)


	def addcontents_frame1(self, query_cells):
		self.canvas_labelframe1 = Canvas(self.labelframe1, width=150, height=300)
		self.frame_labelframe1 = Frame(self.canvas_labelframe1)
		self.vsb_framelabel1 = Scrollbar(self.labelframe1, orient='vertical', \
			command=self.canvas_labelframe1.yview)
		self.canvas_labelframe1.configure(yscrollcommand=self.vsb_framelabel1.set)
		self.vsb_framelabel1.pack(side='right', fill='y')
		self.canvas_labelframe1.create_window((0,0), window=self.frame_labelframe1, anchor='nw')
		self.canvas_labelframe1.pack(side="left", fill="both", expand=True)
		self.frame_labelframe1.bind('<Configure>', self.onFrame1Configure)
		count=0
		for i in range(int(query_cells)):
			exec("self.query_button_" + str(i) + \
				" = Button(self.frame_labelframe1, text='Data " +str(self.sample_each_data[i])+"- Query " + \
				str(i) + "', command= lambda self = self , query_number = "+str(i)+" : self.exec_1("+str(i)+"))")

			exec("self.query_button_" + str(i) + ".grid(row=" + str(i) + \
				", column=0, padx=20, pady=5)")


	def addcontents_frame2(self):
		self.canvas_labelframe2 = Canvas(self.labelframe2, width=750, height=200)
		self.frame_labelframe2 = Frame(self.canvas_labelframe2)
		self.vsb_framelabel2 = Scrollbar(self.labelframe2, orient='vertical', \
			command=self.canvas_labelframe2.yview)
		self.canvas_labelframe2.configure(yscrollcommand=self.vsb_framelabel2.set)
		self.vsb_framelabel2.pack(side='right', fill='y')
		self.canvas_labelframe2.create_window((0,0), window=self.frame_labelframe2, anchor='nw')
		self.canvas_labelframe2.pack(side="left", fill="both", expand=True)
		self.frame_labelframe2.bind('<Configure>', self.onFrame2Configure)
		Label(self.frame_labelframe2, text="Cell ID", font=("Helvetica", 12))\
		.grid(row=0, column=0, pady=10, padx=5)
		Label(self.frame_labelframe2, text="Experiment ID", font=("Helvetica", 12))\
		.grid(row=0, column=1, pady=10, padx=5)
		Label(self.frame_labelframe2, text="Phenotype", font=("Helvetica", 12))\
		.grid(row=0, column=2, pady=10, padx=5)
		Label(self.frame_labelframe2, text="P value", font=("Helvetica", 12))\
		.grid(row=0, column=3, pady=10, padx=5)
		Label(self.frame_labelframe2, text="Adjusted P Value", font=("Helvetica", 12))\
		.grid(row=0, column=4, pady=10, padx=5)

	def onFrame1Configure(self, event):
		self.canvas_labelframe1.configure(scrollregion=self.canvas_labelframe1.bbox("all"))

	def onFrame2Configure(self, event):
		self.canvas_labelframe2.configure(scrollregion=self.canvas_labelframe2.bbox("all"))

	def update_information_exp(self, query_number):
		self.canvas_labelframe2.destroy()
		self.vsb_framelabel2.destroy()
		self.addcontents_frame2()
		self.labelframe2.config(text=' Nearest Neighbours in RNA-Seq : Query Cell ' + str(query_number) + ' ')
		r=1
		temp_filename_exp= 'storage/query_exp_' + str(query_number) + '.txt'
		with open(temp_filename_exp) as file1:
			reader1 = csv.reader(file1, delimiter="@")
			count = 0
			iterator = 0
			for col in reader1:
				if count == 0:
					prev_study = col[1]
				if count > 0 and col[1] == prev_study:
					iterator = iterator + 1
				elif count > 0:
					iterator = 0
				count = count + 1
				prev_study = col[1]
				if iterator < self.results_per_study:
					for c in range(len(col)):
						if c == 2:
							phenotype_text = col[c]
							Label(self.frame_labelframe2, text=phenotype_text, \
								wraplength = 400).grid(row=r, column=c, pady = 5, sticky=W)
						else:
							self.l = Label(self.frame_labelframe2, text=col[c],\
								foreground="blue", cursor="hand2")
							self.l.grid(row=r, column=c, pady = 5)
							# self.l.bind("<Button-1>", \
							# 	lambda event, self=self, study_id = col[c]: \
							# 	self.studyid_callback(event, study_id))
					r = r+1

	def studyid_callback(self, event, study_id):
		# https://trace.ddbj.nig.ac.jp/DRASearch/study?acc=
		url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=' + str(study_id)
		# webbrowser.open_new(url)
		pid = '"' + str(study_id) + '"'
		# print (project_abstracts[pid])
		top = Toplevel(self.master)
		top.title(pid + ' ABSTRACT')
		Message(top, text=project_abstracts[pid]).pack()

	def onFrame4Configure(self, event):
		self.canvas_labelframe4.configure(scrollregion=self.canvas_labelframe4.bbox("all"))

	def exec_2(self,i):
		self.Genes(i)

	def Genes(self,j):
		self.chat2 = ttk.Treeview(self.newFrame3, height=12,columns=("Genes","Marker of Celltype"))
		self.chat2.heading('#0', text='GENES')
		self.chat2.column('#0', width=150,minwidth=50)

		self.chat2.heading('#1', text='Marker of Celltype')
		self.chat2.column('#1', width=700,minwidth=300)

		p2=np.loadtxt("marker_gene.txt",delimiter=" ",dtype='str')
		p2=np.array(p2)

		p3=np.loadtxt("celltype_markers.txt",delimiter="\t",dtype='str')
		p3=np.array(p3)
		
		self.chat2.grid(row=0,column=0)

		for i in range(len(p2)):
			if j==i:
				data1=p2[j].tolist()
				data2=p3[j].tolist()
				for q,j in zip(data1,data2):
					self.chat2.insert('','end',text = q,valu=[j])
	
	def addcontents_frame4(self, query_cells):
		self.canvas_labelframe4 = Canvas(self.labelframe4, width=150, height=500)
		self.frame_labelframe4 = Frame(self.canvas_labelframe4)
		self.vsb_framelabel4 = Scrollbar(self.labelframe4, orient='vertical', \
			command=self.canvas_labelframe1.yview)
		self.canvas_labelframe4.configure(yscrollcommand=self.vsb_framelabel4.set)
		self.vsb_framelabel4.pack(side='right', fill='y')
		self.canvas_labelframe4.create_window((0,0), window=self.frame_labelframe4, anchor='nw')
		self.canvas_labelframe4.pack(side="left", fill="both", expand=True)
		self.frame_labelframe4.bind('<Configure>', self.onFrame4Configure)
		for i in range(int(query_cells)):
			exec("self.query_button_" + str(i) + \
				" = Button(self.frame_labelframe4, text='Query Cell " + \
				str(i) + "', command= lambda self = self , query_number = "+str(i)+" : self.exec_2("+str(i)+"))")

			exec("self.query_button_" + str(i) + ".grid(row=" + str(i) + \
				", column=0, padx=20, pady=5)")

	def display_wordcloud(self):
		novi = Toplevel()
		canvas = Canvas(novi, width = 1400, height = 700)
		canvas.pack(expand = YES, fill = BOTH)
		gif1 = PhotoImage(file = './wordcloud.png')
	    #image not visual
		obj1 = canvas.create_image(50, 50, image = gif1, anchor = "nw")
		obj2 = canvas.create_text(20, 20, text='Click on Image anywhere to view top 50 enriched Genes for every Query', anchor = "nw")
	    #assigned the gif1 to the canvas object
		canvas.gif1 = gif1
		canvas.tag_bind(obj1, '<ButtonPress-1>', self.onObjectClick)  		

	def onObjectClick(self, event):                  
		win = Toplevel()
		win.wm_title("Enriched Genes")
		w, h = win.winfo_screenwidth(), win.winfo_screenheight()
		win.geometry("%dx%d+0+0" % (1400, 700))
		win.resizable(0,0)
		print(self.query_cells)
		Frame.__init__(self, win)
		self.grid()
		self.labelframe4 = LabelFrame(win, text="Query Cells")
		self.labelframe4.grid(row=1, column=0, pady=5, padx=20)
		self.addcontents_frame4(self.query_cells)

		self.newFrame3=Frame(win,height=80,width=100)
		self.newFrame3.grid(row=1,column=1,padx=0, pady=20,sticky="E")
		self.newFrame3.grid_rowconfigure(1, weight=1)
		self.newFrame3.grid_columnconfigure(0, weight=3)

		win.protocol("WM_DELETE_WINDOW", self.close_window)

	
	def onFrame5Configure(self, event):
		self.canvas_labelframe5.configure(scrollregion=self.canvas_labelframe5.bbox("all"))

	def exec_3(self,i):
		self.Genes_enhancers(i)

	def Genes_enhancers(self,j):
		self.chat2 = ttk.Treeview(self.newFrame2, height=12,columns=("Genes"))
		self.chat2.heading('#0', text='GENES')
		self.chat2.column('#0', width=200,minwidth=50)

		p2=np.loadtxt("enhancers.txt",delimiter=" ",dtype='str')
		p2=np.array(p2)
		
		self.chat2.grid(row=0,column=0)

		for i in range(len(p2)):
			if j==i:
				data1=p2[j].tolist()
				data1 = np.unique(data1)
				for q in zip(data1):
					self.chat2.insert('','end',text = q)
	
	
	def addcontents_frame5(self, query_cells):
		self.canvas_labelframe5 = Canvas(self.labelframe5, width=150, height=500)
		self.frame_labelframe5 = Frame(self.canvas_labelframe5)
		self.vsb_framelabel5 = Scrollbar(self.labelframe5, orient='vertical', \
			command=self.canvas_labelframe5.yview)
		self.canvas_labelframe5.configure(yscrollcommand=self.vsb_framelabel5.set)
		self.vsb_framelabel5.pack(side='right', fill='y')
		self.canvas_labelframe5.create_window((0,0), window=self.frame_labelframe5, anchor='nw')
		self.canvas_labelframe5.pack(side="left", fill="both", expand=True)
		self.frame_labelframe5.bind('<Configure>', self.onFrame5Configure)
		for i in range(int(query_cells)):
			exec("self.query_button_" + str(i) + \
				" = Button(self.frame_labelframe5, text='Data "+str(self.sample_each_data[i])+" - Query " + \
				str(i) + "', command= lambda self = self , query_number = "+str(i)+" : self.exec_3("+str(i)+"))")

			exec("self.query_button_" + str(i) + ".grid(row=" + str(i) + \
				", column=0, padx=20, pady=5)")


	def display_enhancers(self):
		win1 = Toplevel()
		win1.wm_title("Enhancers Genes")
		w, h = win1.winfo_screenwidth(), win1.winfo_screenheight()
		win1.geometry("%dx%d+0+0" % (1400, 700))
		win1.resizable(0,0)
		print(self.query_cells)
		Frame.__init__(self, win1)
		self.grid()
		self.labelframe5 = LabelFrame(win1, text="Query Cells")
		self.labelframe5.grid(row=1, column=0, pady=5, padx=20)
		self.addcontents_frame5(self.query_cells)

		self.newFrame2=Frame(win1,height=80,width=100)
		self.newFrame2.grid(row=1,column=1,padx=0, pady=20,sticky="E")
		self.newFrame2.grid_rowconfigure(1, weight=1)
		self.newFrame2.grid_columnconfigure(0, weight=3)

		win1.protocol("WM_DELETE_WINDOW", self.close_window)

	def update_graph(self):
		print(int(self.cls_sp.get()))
		try: 
			self.canvas.get_tk_widget().pack_forget()
		except AttributeError: 
			pass 
		clustering = SpectralClustering(n_clusters=int(self.cls_sp.get()),assign_labels="discretize",random_state=0,affinity='precomputed').fit(self.adj)
		ann = []
		for i in range(self.query_cells):
    				ann.append("D"+str(self.sample_each_data[i])+" "+"Q"+str(i))

		colors = np.array(list(islice(cycle(['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']),int(max(clustering.labels_) + 1))))
		# add black color for outliers (if any)
		colors = np.append(colors, ["#000000"])
		G = nx.from_numpy_matrix(np.matrix(self.adj))
		color=colors[clustering.labels_]
		pos=nx.spring_layout(G)	
		nx.set_node_attributes(G, self.data_exp)

		labels = {}    
		for i,node in enumerate(G.nodes()):
			labels[node] = ann[i]

		fig, ax = plt.subplots(figsize=(5,5))

		nodes = nx.draw_networkx_nodes(G, pos=pos, ax=ax,node_color=color,alpha=0.8,cmap=plt.get_cmap('jet'),node_size=50)
		#nx.draw_networkx_labels(G, pos, labels, font_size=9)
		#nx.draw_networkx_edges(G, pos=pos, ax=ax)	

		annot = ax.annotate("", xy=(0,0), xytext=(5,5),textcoords="offset points",bbox=dict(boxstyle="round", fc="w"),arrowprops=dict(arrowstyle="->"))
		annot.set_visible(False)

		def update_annot(ind):
			nod = ind["ind"][0]
			node = list(pos.keys())[nod]
			xy = pos[node]
			annot.xy = xy
			node_attr = {'node': node}
			node_attr.update(G.nodes[node])
			text = '\n'.join('{k}: {v}' for k, v in node_attr.items())
			annot.set_text("ID : "+str(ann[nod])+"\nHits:"+str(G.nodes[node]["matches"]))
		
		def hover(event):
			vis = annot.get_visible()
			if event.inaxes == ax:
				cont, ind = nodes.contains(event)
				if cont:
					update_annot(ind)
					annot.set_visible(True)
					fig.canvas.draw_idle()
				else:
					if vis:
						annot.set_visible(False)
						fig.canvas.draw_idle()

		self.canvas = FigureCanvasTkAgg(fig, master=self.novi1)
		self.canvas.draw()
		self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

		#toolbar = NavigationToolbar2Tk(self.canvas, self.novi1 )
		#toolbar.update()
		#self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)	

		fig.canvas.mpl_connect("motion_notify_event", hover)
		

	def display_graph(self):
		self.novi1 = Toplevel()
		self.novi1.geometry("700x700") 		

		#Metadata Dict
		self.data_exp = {}
		for i in range(int(self.query_cells)):	
			self.data_exp[i] = {}
			temp_filename = './storage/query_exp_'+ str(i) + '.txt'
			with open(temp_filename) as file:
				reader = list(csv.reader(file,delimiter="@"))
				reader = np.array(reader)[:,2]
			self.data_exp[i]["matches"] = list(reader)

		#ADJACENCY / GRAPH 
		top_clusters = np.loadtxt('./top_clusters.txt',delimiter=" ",dtype='str')
		self.adj = np.zeros((top_clusters.shape[1], top_clusters.shape[1]))
		for i in range(top_clusters.shape[1]):
			for j in range(i+1, top_clusters.shape[1]):
				self.adj[i, j] = len(np.intersect1d(top_clusters[:,i],top_clusters[:,j]))
		#clus = top_clusters.ravel()
		#res = [t.replace('_', '.') for t in clus] 
		#res = np.array(res).reshape(top_clusters.shape[0],top_clusters.shape[1])
		#adj = kneighbors_graph(np.transpose(res), 5, mode='distance', include_self=True)
		#adj = adj.toarray()
		#adj = np.nan_to_num(adj)
		#adj = np.abs(adj)                                   # no negative weights
		#adj = adj - np.diag(np.diag(adj))                   # no self loops
		#adj = np.triu(adj) + np.transpose(np.triu(adj))
		clustering = SpectralClustering(n_clusters=2,assign_labels="discretize",random_state=0,affinity='precomputed').fit(self.adj)
		ann = []
		for i in range(self.query_cells):
    			ann.append("D"+str(self.sample_each_data[i])+" "+"Q"+str(i))

		colors = np.array(list(islice(cycle(['#377eb8', '#ff7f00', '#4daf4a','#f781bf', '#a65628', '#984ea3','#999999', '#e41a1c', '#dede00']),int(max(clustering.labels_) + 1))))
		# add black color for outliers (if any)
		colors = np.append(colors, ["#000000"])
		G = nx.from_numpy_matrix(np.matrix(self.adj))
		color=colors[clustering.labels_]
		pos=nx.spring_layout(G)	
		nx.set_node_attributes(G, self.data_exp)

		labels = {}    
		for i,node in enumerate(G.nodes()):
			labels[node] = ann[i]

		fig, ax = plt.subplots(figsize=(5,5))

		nodes = nx.draw_networkx_nodes(G, pos=pos, ax=ax,node_color=color,alpha=0.8,cmap=plt.get_cmap('jet'),node_size=50)
		#nx.draw_networkx_labels(G, pos, labels, font_size=9)
		#nx.draw_networkx_edges(G, pos=pos, ax=ax)	

		annot = ax.annotate("", xy=(0,0), xytext=(0.5,0.5),textcoords="offset points",bbox=dict(boxstyle="round", fc="0.8"),arrowprops=dict(arrowstyle="->"))
		annot.set_visible(False)
		
		def update_annot(ind):
			nod = ind["ind"][0]
			node = list(pos.keys())[nod]
			xy = pos[node]
			annot.xy = xy
			node_attr = {'node': node}
			node_attr.update(G.nodes[node])
			text = '\n'.join('{k}: {v}' for k, v in node_attr.items())
			annot.set_text("ID : "+str(ann[nod])+"\nHits:"+str(G.nodes[node]["matches"]))
		
		def hover(event):
			vis = annot.get_visible()
			if event.inaxes == ax:
				cont, ind = nodes.contains(event)
				if cont:
					update_annot(ind)
					annot.set_visible(True)
					fig.canvas.draw_idle()
				else:
					if vis:
						annot.set_visible(False)
						fig.canvas.draw_idle()

		self.canvas = FigureCanvasTkAgg(fig, master=self.novi1)
		self.canvas.draw()
		self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

		toolbar = NavigationToolbar2Tk(self.canvas, self.novi1 )
		toolbar.update()
		self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)	

		fig.canvas.mpl_connect("motion_notify_event", hover)

		labelframe=LabelFrame(self.novi1, text="No of Spectral Clusters")
		labelframe.pack(side=BOTTOM)

		#cls=IntVar(None)
		self.cls_sp=Entry(labelframe,width=10)
		self.cls_sp.pack(side="left")
		self.cls_sp.insert(0,"2")

		button = tk.Button(labelframe, text='SUBMIT', command=self.update_graph)
		button.pack(side="right")

		#WORD CLOUD
		canvas = Canvas(self.novi1 , width = 300, height = 300)
		canvas.pack(expand = YES, fill = BOTH)
		gif1 = PhotoImage(file = './wordcloud.png')
	    	#image not visual
		obj1 = canvas.create_image(50, 50, image = gif1, anchor = "nw")
	    	#assigned the gif1 to the canvas object
		canvas.gif1 = gif1

		#GENES
		self.chat2 = ttk.Treeview(self.newFrame3, height=12,columns=("Genes"))
		self.chat2.heading('#0', text='GENES')
		self.chat2.column('#0', width=150,minwidth=50)

		p2=np.loadtxt("marker_gene.txt",delimiter=" ",dtype='str')
		p2=np.array(p2)
		
		self.chat2.grid(row=0,column=0)

		data1=p2[0].tolist()
		for q,j in data1:
			self.chat2.insert('','end',text = q,valu=[j])

		
		#canvas1 = Canvas(novi1, width = 1400, height = 700)
		#canvas1.pack(expand = YES, fill = BOTH)
		#gif1 = PhotoImage(file = './graph.png')
		#obj1 = canvas1.create_image(50, 50, image = gif1, anchor = "nw")
		#canvas1.gif1 = gif1


