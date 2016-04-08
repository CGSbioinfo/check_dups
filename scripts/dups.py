#!/usr/bin/env python
import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
from joblib import Parallel, delayed
import multiprocessing
import subprocess
sys.path.insert(0,'/usr/local/bin/')
import functions
import argparse

def read_sample_names(x):
    sampleNames = []
    sample_names_file = open(x,'r')
    for line in sample_names_file:
        sampleNames.append(line.strip())
    return(sampleNames)

def markdups_picard(i):
	bamfiles = os.listdir(in_dir)
	bamfiles = [bamfiles[y] for y, x in enumerate(bamfiles) if re.findall(i, x)]
	bamfiles = [bamfiles[y] for y, x in enumerate(bamfiles) if re.findall(r'sortedByCoord\.out\.bam$', x)][0]
	print bamfiles
	os.system("java -jar ~/tools/picard-tools-1.127/picard.jar MarkDuplicates " +
		" INPUT=" + in_dir + "/" + bamfiles +
		" OUTPUT=" + out_dir + "/" + i + "_marked_dups.bam"
		" M=" + out_dir + "/" + i + "dup_metrics.txt"
		" REMOVE_DUPLICATES=FALSE ASSUME_SORTED=TRUE")

def extract_dups(i):
	os.system("samtools view -bf 0x400 " + out_dir + '/' + i + "_marked_dups.bam > " + out_dir + "/" + i + "_dups.bam")
	os.system("samtools index " + out_dir + "/" + i + "_dups.bam")

#########################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Counting Reads')
	parser.add_argument('--in_dir', help='Path to folder containing bam files. Default=alignedReads/', default='alignedReads/')
	parser.add_argument('--out_dir', help='Path to output folder. Default=mark_dups/', default='mark_dups/')
	parser.add_argument('--sample_names', help='Text file with sample names to process. Default=sample_names_dups.txt', default='sample_names_dups.txt')

	args=parser.parse_args()
	in_dir=args.in_dir
	sample_names_file=args.sample_names
	out_dir=args.out_dir
	functions.make_sure_path_exists(out_dir)

	path=os.getcwd()
	print path
	os.chdir(path)
	
	# Read sample names text file
	sampleNames = read_sample_names(sample_names_file)
	

	# Mark duplicates 
	#Parallel(n_jobs=2)(delayed(markdups_picard)(i) for i in sampleNames)

	# Create a bam with just duplicates
	Parallel(n_jobs=3)(delayed(extract_dups)(i) for i in sampleNames)