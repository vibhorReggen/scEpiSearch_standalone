import sys
from test_final_v2 import test_final_v2
import numpy as np

def main(argv):
    chr_file = argv[0]
    query_file = argv[1]
    top_study = int(argv[2])
    query_type = int(argv[3])
    acc_fast = int(argv[4])
    active_poised = int(argv[5])
    cluster = int(argv[6])
    form_type = str(argv[7])
    anno_unanno = int(argv[8])
    test_final_v2(chr_file,query_file,top_study,query_type,acc_fast,active_poised,cluster,form_type,anno_unanno)

if __name__ == "__main__":
	main(sys.argv[1:])
