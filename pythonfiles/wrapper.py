import csv
import sys
import os
import zipfile
import gzip
import shutil
import glob
def check_data_validity_chr(file):
    with open(file, "r") as csvfile:
        try:
            dialect = csv.Sniffer().sniff(csvfile.read(), delimiters = "\t")
            return 1
        except:
            return 0

def check_data_validity_count(file):
    with open(file, "r") as csvfile:
        try:
            dialect = csv.Sniffer().sniff(csvfile.read(), delimiters = ",")
            return 1
        except:
            print("Wrong Delimiter")
            return 0



chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM','chrUn']
chr_path = sys.argv[1]
count_path = sys.argv[2]

chr_path_f = chr_path
count_path_f = count_path

new = check_data_validity_chr(chr_path_f)
new1 = check_data_validity_count(count_path_f)
n2 = 0
n3 = 0
n4 = False
line1 = sum(1 for line in open(count_path_f))
line2 = sum(1 for line in open(chr_path_f))
if (line1 == line2):
	n2 = 1

with open(chr_path_f) as f:
	reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
	first_row = next(reader)
	num_cols = len(first_row)
if(num_cols == 3):
	n3 = 1

l1 = list()
with open(chr_path_f) as f:
	for line in f:
		l1.append(line.split("\t")[0])

n4 = set(l1) <= set(chr_list)
if(new == 1 and new1 == 1 and n2 == 1 and n3 == 1 and n4 == True):
	sys.exit(0)
else:
	sys.exit(1)
