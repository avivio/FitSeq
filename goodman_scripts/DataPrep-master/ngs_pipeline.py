__author__ = 'dell7'


BIN_DIR = '/home/labs/pilpel/avivro/ngs_pipeline/data/Project_avivro/'

def go_over_bins(bin_dir=BIN_DIR):


    for root, dirs, files in os.walk(fq_dir, followlinks=True):
            print root
            print dirs
            pring files
