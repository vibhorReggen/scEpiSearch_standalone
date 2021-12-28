## results page
import subprocess, csv, webbrowser, os, pickle, re, matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from pdf_report_generator import generate_pdf
matplotlib.use('Agg')
plt.switch_backend('agg')
# for Python3
from tkinter import *
from tkinter.ttk import *
from tkinter import messagebox
import tkinter as tk
from tkinter import ttk
import PIL.Image
import PIL
#from PIL import ImageTk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


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



class ResultPage_old(Frame):
	def __init__(self, master, query_cells, results_per_study, query_type):
		self.master = master
		self.results_per_study = results_per_study
		self.query_cells = query_cells
		self.query_type = query_type
		Frame.__init__(self, self.master)
		self.grid()
		# logo
		
		p1=pd.read_csv("tsne.txt" ,delimiter=' ')
		p1=np.array(p1)
		self.overlap = np.loadtxt("./overlap.txt",dtype='str')


		x=p1[:,3]
		y=p1[:,4]
		color=np.array(p1[:,1])
		color=color.tolist()
		cmap = plt.cm.Spectral
		norm = plt.Normalize(vmin=4, vmax=5)

		lab=p1[:,0]
		fig = Figure(figsize=(5,3))
		a=fig.add_subplot(111)
		a.scatter(x,y,c=color,picker=1)
		
		img = PhotoImage(file="logo.png")
		imglabel = Label(self.master, image=img)
		imglabel.image = img # keep a reference!
		imglabel.grid(row=0, column=1, sticky='nw', padx=2, pady=0)

		self.labelframe1 = LabelFrame(self.master, text="Query Cells")
		self.labelframe1.grid(row=1, column=0, pady=1, padx=2,sticky="nw")

		self.labelframe2 = LabelFrame(self.master, text="Phenotype")
		self.labelframe2.grid(row=1, column=1, columnspan=1, pady=1, padx=85)
		#self.labelframe3 = LabelFrame(self.master, text="Phenotype")
		#self.labelframe3.grid(row=2, column=1, columnspan=1, pady=1, padx=85)

		self.frame1=Frame(self.master)
		self.frame1.grid(row=0,column=0,padx=2,pady=2)
		self.newFrame=Frame(self.master,height=70,width=600,border=1)
		self.newFrame.grid(row=0,column=1,padx=20,pady=0,sticky="WS")
		canvas = FigureCanvasTkAgg(fig,self.frame1)
		c1=fig.canvas.mpl_connect('pick_event',self.onpick)
		canvas.get_tk_widget().pack()
		canvas.draw()		

		self.addcontents_frame1(query_cells)
		self.addcontents_frame2()
		#self.addcontents_frame3()
		self.update_information_exp('0') ## display cell 0's information
		#self.update_information_epi('0')
		#Button(self.master, text="Close", command=self.close_window).grid(\
		#	row=2, column=0)
		#Button(self.master, text="View Report", command=self.generate_report).grid(\
		#	row=1, column=0, sticky=E, padx=5)
		
		#Button(self.master, text="Show WordCloud", command=self.display_wordcloud).grid(\
		#	row=2, column=2, sticky=W, padx=5)
		#Button(self.master, text="Show Enhancers", command=self.display_enhancers).grid(\
		#	row=3, column=1, sticky=W, padx=5)

		self.b2=Button(self.master, text="View Report", command=self.generate_report)
		self.b2.grid(row=1, column=0, padx=1)
		self.b2.place(x=260,y=320)
		
		self.b4=Button(self.master, text="Show WordCloud", command=self.display_wordcloud)
		self.b4.grid(row=1, column=0, padx=1)
		self.b4.place(x=260,y=400)

		self.b5=Button(self.master, text="Show Gene Frequency Plot", command=self.display_gene_frequency_bar)
		self.b5.grid(row=1, column=0, padx=1)
		self.b5.place(x=260,y=460)

		self.b8=Button(self.master, text="Show Gene Enrichment Table", command=self.gene_enrichment_table)
		self.b8.grid(row=1, column=0, padx=1)
		self.b8.place(x=260,y=510)


		self.b6 = Button(self.master, text="Fraction overlap with accessibility peak list :")
		self.b6.grid(row=1, column=0, padx=1)
		self.b6.place(x=20,y=600)

		self.b7 = Button(self.master, text=self.overlap)
		self.b7.grid(row=1, column=0, padx=1)
		self.b7.place(x=360,y=600)

		#self.master.protocol("WM_DELETE_WINDOW", self.close_window_main)

	def onpick(self,event):
		xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
		ind = event.ind

		p1=pd.read_csv("tsne.txt" ,delimiter=' ')
		p1=np.array(p1)
		self.chat1 = ttk.Treeview(self.newFrame, height=3,columns=('EXPRE_STUDY_ID','EXPR_CELL_ID','EXPR_CELL_TYPE'))
		self.chat1.heading('#0', text='QUERY_NO')
		self.chat1.heading('#1',text='EXPR_STUDY_ID')
		self.chat1.heading('#2',text='EXPR_CELL_ID')
		self.chat1.heading('#3',text='EXPR_CELL_TYPE')
		self.chat1.column('#0', width=150,minwidth=100)
		self.chat1.column('#1', width=150,minwidth=100)
		self.chat1.column('#2', width=150,minwidth=100)
		self.chat1.column('#3', width=150,minwidth=100)
		for i in ind:
			data=p1[i,:]
			cluster=str(data[2])
			clus=[]
			for j in range(len(p1)):
				cl=str(p1[j][2])
				if cl==cluster and i!=j:
					clus.append(p1[j][0])

			clus.append(p1[i,0])

			for i in range(len(clus)):
				cl1=str(clus[i])
				d2=np.loadtxt("storage/query_exp_"+cl1+".txt",delimiter="@",dtype='str')
				values1=d2[0,:]
				query="QUERY"+cl1
				STUDY_ID_EXP=values1[1]
				CELL_ID_EXP=values1[0]
				CELL_TYPE_EXP=values1[2]
				self.chat1.insert('','end',text=query,values=(STUDY_ID_EXP,CELL_ID_EXP,CELL_TYPE_EXP))

			self.xscrollbar = Scrollbar(self.newFrame,orient=HORIZONTAL)
			self.xscrollbar.grid(row=1,column=0,sticky = (N,S,W,E))
			self.chat1.config(xscrollcommand=self.xscrollbar.set)
			self.xscrollbar.config(command=self.chat1.xview)
			self.chat1.grid(row=7, column=0)
			if len(p1)>20:
			   #self.update_information_epi(ind[0])
			   self.update_information_exp(ind[0])
			   #self.Genes(ind[0])


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
		#os.remove('./acc_score.csv')
		#os.remove('./celltype_markers.txt')
		#os.remove('./correlation_matrix.txt')
		#os.remove('./enhancers.txt')
		#os.remove('./enrichment_scores.txt')
		#os.remove('./exp.txt')
		#os.remove('./epi.txt')
		#os.remove('./pval_exp.txt')
		#os.remove('./pval_epi.txt')
		#os.remove('./fdr_exp.txt')
		#os.remove('./fdr_epi.txt')
		#os.remove('./foreground.csv')
		#os.remove('./marker_gene.txt')
		#os.remove('./overlap.txt')
		#os.remove('./tsne.txt')
		#os.remove('./wordcloud.png')
		#for i in range(int(self.query_cells)):
		#   os.remove('storage/query_exp_' + str(i) + '.txt')
		#   os.remove('storage/query_epi_' + str(i) + '.txt')
		#os.remove('./heatmap.jpg')
		#os.remove('./heatmap.png')
		
	def exec_1(self,i):
		self.update_information_exp(i)


	def addcontents_frame1(self, query_cells):
		self.canvas_labelframe1 = Canvas(self.labelframe1, width=150, height=200)
		self.frame_labelframe1 = Frame(self.canvas_labelframe1)
		self.vsb_framelabel1 = Scrollbar(self.labelframe1, orient='vertical', \
			command=self.canvas_labelframe1.yview)
		self.canvas_labelframe1.configure(yscrollcommand=self.vsb_framelabel1.set)
		self.vsb_framelabel1.pack(side='right', fill='y')
		self.canvas_labelframe1.create_window((0,0), window=self.frame_labelframe1, anchor='nw')
		self.canvas_labelframe1.pack(side="left", fill="both", expand=True)
		self.frame_labelframe1.bind('<Configure>', self.onFrame1Configure)
		for i in range(int(query_cells)):
			exec("self.query_button_" + str(i) + \
				" = Button(self.frame_labelframe1, text='Query Cell " + \
				str(i) + "', command= lambda self = self , query_number = "+str(i)+" : self.exec_1("+str(i)+"))")

			exec("self.query_button_" + str(i) + ".grid(row=" + str(i) + \
				", column=0, padx=20, pady=5)")


	def addcontents_frame2(self):
		self.canvas_labelframe2 = Canvas(self.labelframe2, width=750, height=500)
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

	def gene_enrichment_table(self):
		win2 = Toplevel()
		win2.wm_title("Enriched Genes")
		w, h = win2.winfo_screenwidth(), win2.winfo_screenheight()
		win2.geometry("%dx%d+0+0" % (1400, 700))
		win2.resizable(0,0)
		print(self.query_cells)
		Frame.__init__(self, win2)
		self.grid()
		self.labelframe4 = LabelFrame(win2, text="Query Cells")
		self.labelframe4.grid(row=1, column=0, pady=5, padx=20)
		self.addcontents_frame4(self.query_cells)

		self.newFrame3=Frame(win2,height=80,width=100)
		self.newFrame3.grid(row=1,column=1,padx=0, pady=20,sticky="E")
		self.newFrame3.grid_rowconfigure(1, weight=1)
		self.newFrame3.grid_columnconfigure(0, weight=3)
		
		win2.protocol("WM_DELETE_WINDOW", self.close_window)


	def Genes(self,j):
		self.chat2 = ttk.Treeview(self.newFrame3, height=12,columns=("Genes","Marker of Celltype","Enrichment Scores"))
		self.chat2.heading('#0', text='GENES')
		self.chat2.column('#0', width=150,minwidth=50)

		self.chat2.heading('#1', text='Marker of Celltype')
		self.chat2.column('#1', width=350,minwidth=250)

		self.chat2.heading('#2', text='Enrichment Scores')
		self.chat2.column('#2', width=150,minwidth=50)	

		p2=np.loadtxt("marker_gene.txt",delimiter=" ",dtype='str')
		p2=np.array(p2)

		p3=np.loadtxt("celltype_markers.txt",delimiter="\t",dtype='str')
		p3=np.array(p3)

		p4=np.loadtxt("marker_gene_enrichment.txt",delimiter=" ",dtype='str')
		p4=np.array(p4)
		
		self.chat2.grid(row=0,column=0)

		for i in range(len(p2)):
			if j==i:
				data1=p2[j].tolist()
				data2=p3[j].tolist()
				data3=p4[j].tolist()
				for q,j,p in zip(data1,data2,data3):
					self.chat2.insert('','end',text = q,values=([j],[p]))
	
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
				" = Button(self.frame_labelframe5, text='Query Cell " + \
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


	def display_gene_frequency_bar(self):
		novi_f = Toplevel()
		canvas_f = tk.Canvas(novi_f, width = 1400, height = 700)
		canvas_f.pack(expand = YES, fill = BOTH)
		gif_f = PhotoImage(file = './gene_frequency.png')
		#image not visual
		canvas_f.gif_f = gif_f
		canvas_f.create_image(60, 60, image = canvas_f.gif_f, anchor = NW)
		canvas_f.create_text(20, 20, text='Frequency of Gene Bar Plot.', anchor = "nw")




		
		













