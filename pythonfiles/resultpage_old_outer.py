## results page
import subprocess, csv, webbrowser, os, pickle, re, matplotlib
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

class ResultPage_old(Frame):

	def __init__(self, master, query_cells, results_per_study):
		self.master = master
		self.results_per_study = results_per_study
		self.query_cells = query_cells
		Frame.__init__(self, self.master)
		self.grid()
		# logo
		p1=pd.read_csv("tsne.txt" ,delimiter=' ')
		p1=np.array(p1)


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
		self.labelframe1.grid(row=1, column=0,padx=2,pady=1,sticky="W")


		self.newFrame1=Frame(self.master,height=80,width=100)
		self.newFrame1.grid(row=1,column=0,padx=0, pady=20,sticky="E")
		self.newFrame1.grid_rowconfigure(1, weight=1)
		self.newFrame1.grid_columnconfigure(0, weight=1)
		# self.labelframe1.place(x=13, y=305)
		# self.newFrame1.place(x=260,y=290)
		self.labelframe2 = LabelFrame(self.master, text="Phenotype")
		self.labelframe2.grid(row=1, column=1, columnspan=1, pady=1, padx=85)
		# self.labelframe2.place(x=430,y=303)
		self.labelframe3 = LabelFrame(self.master, text="Phenotype")
		self.labelframe3.grid(row=2, column=1, columnspan=1, pady=1, padx=85)
		# self.labelframe3.place(x=430,y=520)
		self.frame1=Frame(self.master)
		self.frame1.grid(row=0,column=0,padx=2,pady=2)
		self.newFrame=Frame(self.master,height=70,width=600,border=1)
		self.newFrame.grid(row=0,column=1,padx=20,pady=0,sticky="WS")
		# self.newFrame.place(x=w-1150,y=h-730)
		canvas = FigureCanvasTkAgg(fig,self.frame1)
		c1=fig.canvas.mpl_connect('pick_event',self.onpick)
		canvas.get_tk_widget().pack()
		canvas.draw()
		self.addcontents_frame1(query_cells)
		self.addcontents_frame2()
		self.addcontents_frame3()
		self.update_information_exp('0') ## display cell 0's information
		self.update_information_epi('0')
		self.b1=Button(self.master, text="Close", command=self.close_window)
		self.b1.grid(row=2,column=0)
		# self.b1.place(x=13,y=600)

		self.b2=Button(self.master, text="Download Report", command=self.generate_report)
		self.b2.grid(row=2, column=0, sticky=E, padx=2)
		# self.b2.place(x=13,y=640)

		self.b3=Button(self.master, text="Show Heat Map", command=self.display_heatmap)
		self.b3.grid(row=2, column=0, sticky=W, padx=2)
		# self.b3.place(x=140,y=640)
		self.master.protocol("WM_DELETE_WINDOW", self.close_window)

	def onpick(self,event):
		xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
		ind = event.ind

		p1=pd.read_csv("tsne.txt" ,delimiter=' ')
		p1=np.array(p1)
		self.chat1 = ttk.Treeview(self.newFrame, height=3,columns=('EPIG_STUDY_ID', 'EPIG_CELL_ID','EPIG_CELL_TYPE','EXPRE_STUDY_ID','EXPR_CELL_ID','EXPR_CELL_TYPE'))
		self.chat1.heading('#0', text='QUERY_NO')
		self.chat1.heading('#1', text='EPIG_STUDY_ID')
		self.chat1.heading('#2', text='EPIG_CELL_ID')
		self.chat1.heading('#3',text='EPIG_CELL_TYPE')
		self.chat1.heading('#4',text='EXPR_STUDY_ID')
		self.chat1.heading('#5',text='EXPR_CELL_ID')
		self.chat1.heading('#6',text='EXPR_CELL_TYPE')
		self.chat1.column('#0', width=114,minwidth=100)
		self.chat1.column('#1', width=114,minwidth=100)
		self.chat1.column('#2', width=114,minwidth=100)
		self.chat1.column('#3', width=114,minwidth=100)
		self.chat1.column('#4', width=114,minwidth=100)
		self.chat1.column('#5', width=114,minwidth=100)
		self.chat1.column('#6', width=114,minwidth=100)
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
				d1=np.loadtxt("storage/query_epi_"+cl1+".txt",delimiter="@",dtype='str')
				d2=np.loadtxt("storage/query_exp_"+cl1+".txt",delimiter="@",dtype='str')
				values=d1[0,:]
				values1=d2[0,:]
				query="QUERY"+cl1
				STUDY_ID=values[1]
				CELL_ID=values[0]
				CELL_TYPE=values[2]
				STUDY_ID_EXP=values1[1]
				CELL_ID_EXP=values1[0]
				CELL_TYPE_EXP=values1[2]
				self.chat1.insert('','end',text=query,values=(STUDY_ID,CELL_ID,CELL_TYPE,STUDY_ID_EXP,CELL_ID_EXP,CELL_TYPE_EXP))

			self.xscrollbar = Scrollbar(self.newFrame,orient=HORIZONTAL)
			self.xscrollbar.grid(row=1,column=0,sticky = (N,S,W,E))
			self.chat1.config(xscrollcommand=self.xscrollbar.set)
			self.xscrollbar.config(command=self.chat1.xview)
			self.chat1.grid(row=7, column=0)
			if len(p1)>20:
			   self.update_information_epi(ind[0])
			   self.update_information_exp(ind[0])
			   self.Genes(ind[0])

	def generate_report(self):
		filename = generate_pdf(self.query_cells)
		try:
			tkMessageBox.showinfo("PDF Generated", "File : " + filename , parent=self.master)
		except Exception as e:
			messagebox.showinfo("PDF Generateddd", "File : " + filename , parent=self.master)
		subprocess.Popen(['evince', filename])

	def close_window(self):
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

	def clear_files(self):
		pass
		#os.remove('./acc_score.csv')
		#os.remove('./exp.txt')
		#os.remove('./epi.txt')
		#os.remove('./pval_exp.txt')
		#os.remove('./pval_epi.txt')
		#os.remove('./fdr_exp.txt')
		#os.remove('./fdr_epi.txt')
		#os.remove('./foreground.csv')
		#for i in range(int(self.query_cells)):
		#   os.remove('storage/query_exp_' + str(i) + '.txt')
		#   os.remove('storage/query_epi_' + str(i) + '.txt')
		#os.remove('./heatmap.gif')

	def exec_1(self,i):
		self.update_information_epi(i)
		self.update_information_exp(i)
		self.Genes(i)

	def Genes(self,j):
		self.chat2 = ttk.Treeview(self.newFrame1, height=12,columns=())
		self.chat2.heading('#0', text='GENES')
		self.chat2.column('#0', width=150,minwidth=50)
		p2=np.loadtxt("marker_gene.txt",delimiter=" ",dtype='str')
		p2=np.array(p2)
		for i in range(len(p2)):
			if j==i:
			   data=p2[j].tolist()
			   # print(len(data))
			   for q in data:
				   self.chat2.insert('','end',text=q,values=())
		self.chat2.grid(row=0,column=0)

	def addcontents_frame1(self, query_cells):
		self.canvas_labelframe1 = Canvas(self.labelframe1, width=150, height=250)
		self.frame_labelframe1 = Frame(self.canvas_labelframe1)
		self.vsb_framelabel1 = Scrollbar(self.labelframe1, orient='vertical', \
			command=self.canvas_labelframe1.yview)
		self.canvas_labelframe1.configure(yscrollcommand=self.vsb_framelabel1.set)
		self.vsb_framelabel1.pack(side='right', fill='y')
		self.canvas_labelframe1.create_window((0,0), window=self.frame_labelframe1, anchor='nw')
		self.canvas_labelframe1.pack(side="left", fill="both", expand=True)
		self.frame_labelframe1.bind('<Configure>', self.onFrame1Configure)
		p1=pd.read_csv("tsne.txt" ,delimiter=' ')
		if len(p1)<=20:
		   for i in range(int(query_cells)):
			   exec("self.query_button_" + str(i) + \
					 " = Button(self.frame_labelframe1, text='Query Cell " + \
					 str(i) + "', command= lambda self = self , query_number = "+str(i)+" : self.exec_1("+str(i)+"))")

			   exec("self.query_button_" + str(i) + ".grid(row=" + str(i) + \
					 ", column=0, padx=2, pady=2)")

	def addcontents_frame2(self):
		self.canvas_labelframe2 = Canvas(self.labelframe2, width=790, height=190)
		self.frame_labelframe2 = Frame(self.canvas_labelframe2)
		self.vsb_framelabel2 = Scrollbar(self.labelframe2, orient='vertical', \
			command=self.canvas_labelframe2.yview)
		self.canvas_labelframe2.configure(yscrollcommand=self.vsb_framelabel2.set)
		self.vsb_label=Scrollbar(self.labelframe2,orient='horizontal',\
			command=self.canvas_labelframe2.xview)
		self.canvas_labelframe2.configure(xscrollcommand=self.vsb_label.set)
		self.vsb_framelabel2.pack(side='right', fill='y')
		self.vsb_label.pack(side='bottom',fill='x')
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

	def addcontents_frame3(self):
		self.canvas_labelframe3 = Canvas(self.labelframe3, width=790, height=190)
		self.frame_labelframe3 = Frame(self.canvas_labelframe3)
		self.vsb_framelabel3 = Scrollbar(self.labelframe3, orient='vertical', \
			command=self.canvas_labelframe3.yview)
		self.canvas_labelframe3.configure(yscrollcommand=self.vsb_framelabel3.set)
		self.vsb_framelabel = Scrollbar(self.labelframe3, orient='horizontal', \
			command=self.canvas_labelframe3.xview)
		self.canvas_labelframe3.configure(xscrollcommand=self.vsb_framelabel.set)
		self.vsb_framelabel3.pack(side='right', fill='y')
		self.vsb_framelabel.pack(side='bottom', fill='x')
		self.canvas_labelframe3.create_window((0,0), window=self.frame_labelframe3, anchor='nw')
		self.canvas_labelframe3.pack(side="left", fill="both", expand=True)
		self.frame_labelframe3.bind('<Configure>', self.onFrame3Configure)
		Label(self.frame_labelframe3, text="Cell ID", font=("Helvetica", 12))\
		.grid(row=0, column=0, pady=10, padx=5)
		Label(self.frame_labelframe3, text="Experiment ID", font=("Helvetica", 12))\
		.grid(row=0, column=1, pady=10, padx=5)
		Label(self.frame_labelframe3, text="Phenotype", font=("Helvetica", 12))\
		.grid(row=0, column=2, pady=10, padx=5)
		Label(self.frame_labelframe3, text="P value", font=("Helvetica", 12))\
		.grid(row=0, column=3, pady=10, padx=5)
		Label(self.frame_labelframe3, text="Adjusted P Value", font=("Helvetica", 12))\
		.grid(row=0, column=4, pady=10, padx=5)

	def onFrame1Configure(self, event):
		self.canvas_labelframe1.configure(scrollregion=self.canvas_labelframe1.bbox("all"))

	def onFrame2Configure(self, event):
		self.canvas_labelframe2.configure(scrollregion=self.canvas_labelframe2.bbox("all"))

	def onFrame3Configure(self, event):
		self.canvas_labelframe3.configure(scrollregion=self.canvas_labelframe3.bbox("all"))

	def update_information_exp(self, query_number):
		self.canvas_labelframe2.destroy()
		self.vsb_framelabel2.destroy()
		self.vsb_label.destroy()
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
					r = r+1

	def update_information_epi(self, query_number):
		self.canvas_labelframe3.destroy()
		self.vsb_framelabel3.destroy()
		self.vsb_framelabel.destroy()
		self.addcontents_frame3()
		self.labelframe3.config(text=' Nearest Neighbours in Single Cell Epigenome : Query Cell ' + str(query_number) + ' ')
		r=1
		temp_filename_epi= 'storage/query_epi_' + str(query_number) + '.txt'
		with open(temp_filename_epi) as file1:
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
							Label(self.frame_labelframe3, text=phenotype_text, \
								wraplength = 400).grid(row=r, column=c, pady = 5, sticky=W)
						else:
							self.l = Label(self.frame_labelframe3, text=col[c],\
								foreground="blue", cursor="hand2")
							self.l.grid(row=r, column=c, pady = 5)
							self.l.bind("<Button-1>", \
								lambda event, self=self, study_id = col[c]: \
								self.studyid_callback(event, study_id))
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

	def display_heatmap(self):
		novi = Toplevel()
		canvas = tk.Canvas(novi, width = 700, height = 700)
		canvas.pack(expand = YES, fill = BOTH)
		gif1 = PhotoImage(file = './heatmap.png')
		#gif1 = PIL.Image.open('./heatmap.png')
		#canvas.image = ImageTk.PhotoImage(gif1)
		#image not visual
		canvas.gif1 = gif1
		canvas.create_image(20, 20, image = canvas.gif1, anchor = NW)
		#assigned the gif1 to the canvas object
		#canvas.gif1 = gif1
