## pdf report_generator
# from reportlab.pdfgen import canvas
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
import PIL
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import Table, TableStyle
from reportlab.lib.units import cm, inch
from reportlab.lib import colors
from reportlab.lib.units import inch
import matplotlib.pyplot as plt
import numpy as np
from wordcloud import WordCloud, STOPWORDS 
import matplotlib.pyplot as plt 
import csv, time

filename = 'report/results_' + str(time.ctime()) + '.pdf'

def generate_pdf(query_cells,query_type):
	print(query_type)
	doc = SimpleDocTemplate(filename,
                        rightMargin=72,leftMargin=72,
                        topMargin=72,bottomMargin=18)
	Story=[]

	logo = 'logo.jpg'
	im = Image(logo, 6*inch, 1.5*inch)
	Story.append(im)

	full_name = "ScEpiSearch Report"
	address_parts = ["IIIT Delhi", "New Delhi"]
	formatted_time = time.ctime()

	styles=getSampleStyleSheet()
	styles.add(ParagraphStyle(name='Justify', alignment=TA_JUSTIFY))

	ptext = '<font size=12>%s</font>' % formatted_time
	Story.append(Paragraph(ptext, styles["Normal"]))

	Story.append(Spacer(1, 12))

	# Create return address_parts
	ptext = '<font size=12>%s</font>' % full_name
	Story.append(Paragraph(ptext, styles["Normal"]))

	Story.append(Spacer(1, 12))

	for part in address_parts:
	    ptext = '<font size=12>%s</font>' % part.strip()
	    Story.append(Paragraph(ptext, styles["Normal"]))

	Story.append(Spacer(1, 12))
	Story.append(Spacer(1, 12))

	all_cells = [(0, 0), (-1, -1)]
	header = [(0, 0), (-1, 0)]
	column0 = [(0, 0), (0, -1)]
	column1 = [(1, 0), (1, -1)]
	column2 = [(2, 0), (2, -1)]
	column3 = [(3, 0), (3, -1)]
	column4 = [(4, 0), (4, -1)]
	table_style = TableStyle([('VALIGN', all_cells[0], all_cells[1], 'TOP'),
		('FONTSIZE', all_cells[0], all_cells[1], 6),
	    ('LINEBELOW', header[0], header[1], 1, colors.black),
	    ('ALIGN', column0[0], column0[1], 'LEFT'),
	    ('ALIGN', column1[0], column1[1], 'LEFT'),
	    ('ALIGN', column2[0], column2[1], 'LEFT'),
	     ('ALIGN', column3[0], column3[1], 'LEFT'),
    	('ALIGN', column4[0], column4[1], 'LEFT'),
	])

	# PDF Table - Column Widths
	# colWidths = [2 * cm,2* cm, 10 * cm,2.5 * cm,2.5 * cm]

	ptext = '<font size=8>Matching in RNA-Seq Data : </font>'
	Story.append(Paragraph(ptext, styles["Normal"]))

	cells_all = list()
	
	if query_type == 1 or query_type ==2:
		for i in range(int(query_cells)):
			ptext = '<font size=8>Query Cell Number : %s</font>' % str(i)
			Story.append(Paragraph(ptext, styles["Normal"]))
			data_exp = list()
			data_exp.append(['CellId', 'StudyId', 'Phenotype', 'P value', 'Adjusted P value'])
			temp_filename = './storage/query_exp_'+ str(i) + '.txt'
			with open(temp_filename) as file:
				reader = list(csv.reader(file,delimiter="@"))
			for i in range(len(reader)):
				cells_all.append(reader[i][2])	
			reader = data_exp + reader

			for index, row in enumerate(reader):
				for col, val in enumerate(row):
					# if col == 4 or col==0 or index == 0:
						# reader[index][col] = val.strip("'[]")
					# else:
					ptext = '<font size=6>%s</font>' % val
					reader[index][col] = Paragraph(ptext, styles['Normal'])

			t = Table(reader)
			t.setStyle(table_style)
			Story.append(t)
			Story.append(Spacer(1, 12))
			Story.append(Spacer(1, 12))

		Story.append(Spacer(1, 12))
		Story.append(Spacer(1, 12))

		ptext = '<font size=8>Matching in Single Cell Epigenome Data: </font>'
		Story.append(Paragraph(ptext, styles["Normal"]))

		for i in range(int(query_cells)):
			ptext = '<font size=8>Query Cell Number : %s</font>' % str(i)
			Story.append(Paragraph(ptext, styles["Normal"]))
			data_exp = list()
			data_exp.append(['CellId', 'StudyId', 'Phenotype', 'P value', 'Adjusted P value'])
			temp_filename = './storage/query_epi_'+ str(i) + '.txt'
			with open(temp_filename) as file:
				reader = list(csv.reader(file,delimiter="@"))
			for i in range(1,len(reader)):
				cells_all.append(reader[i][2])
			reader = data_exp + reader

			for index, row in enumerate(reader):
				for col, val in enumerate(row):
					# if col != 4 or index == 0:
						# reader[index][col] = val.strip("'[]()")
					# else:
					ptext = '<font size=6>%s</font>' % val
					reader[index][col] = Paragraph(ptext, styles['Normal'])

			t = Table(reader)
			t.setStyle(table_style)
			Story.append(t)
			Story.append(Spacer(1, 12))
			Story.append(Spacer(1, 12))

		heatmap_image = './heatmap.png'
		im = PIL.Image.open(heatmap_image)
		bg = PIL.Image.new("RGB", im.size, (255,255,255))
		bg.paste(im,im)
		bg.save("./heatmap.jpg")
		im = Image("./heatmap.jpg", 6*inch, 1.5*inch)
		Story.append(im)
	else:
		for i in range(int(query_cells)):
			ptext = '<font size=8>Query Cell Number : %s</font>' % str(i)
			Story.append(Paragraph(ptext, styles["Normal"]))
			data_exp = list()
			data_exp.append(['CellId', 'StudyId', 'Phenotype', 'P value', 'Adjusted P value'])
			temp_filename = './storage/query_exp_'+ str(i) + '.txt'
			with open(temp_filename) as file:
				reader = list(csv.reader(file,delimiter="@"))
			for i in range(len(reader)):
				cells_all.append(reader[i][2])	
			reader = data_exp + reader

			for index, row in enumerate(reader):
				for col, val in enumerate(row):
					# if col == 4 or col==0 or index == 0:
						# reader[index][col] = val.strip("'[]")
					# else:
					ptext = '<font size=6>%s</font>' % val
					reader[index][col] = Paragraph(ptext, styles['Normal'])

			t = Table(reader)
			t.setStyle(table_style)
			Story.append(t)
			Story.append(Spacer(1, 12))
			Story.append(Spacer(1, 12))

		Story.append(Spacer(1, 12))
		Story.append(Spacer(1, 12))

		ptext = '<font size=8>Matching in Single Cell Epigenome Data: </font>'
		Story.append(Paragraph(ptext, styles["Normal"]))
		
	doc.build(Story)
	comment_words = '' 
	stopwords = set(STOPWORDS) 
	cells_all = np.array(cells_all)
	#print(cells_all)
	cells_all = cells_all.astype(str)
	comment_words += " ".join(cells_all)+" "
	wordcloud = WordCloud(width = 300, height = 300, background_color ='white', stopwords = stopwords, min_font_size = 10).generate(comment_words) 
  
	# plot the WordCloud image                        
	plt.figure(figsize = (4, 4), facecolor = None) 
	plt.imshow(wordcloud) 
	plt.axis("off") 
	plt.tight_layout(pad = 0) 
	#plt.show()
	plt.savefig('./wordcloud.png')


	 
	return filename
