#!/usr/bin/env python
import os
import gzip
import subprocess
import sys, time, random, re ,requests, logging, glob, tabix
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager
import multiprocessing

logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter=logging.Formatter("%(asctime)s - %(message)s")

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)

sOutputDir="/mnt/towel/BRCA/VCF/Analysis/SKATO/SNUvsTCGAre_withoutAsian"
sCorrectionDir="/mnt/towel/BRCA/VCF/Analysis/SKATO/SNUvsTCGAre_withoutAsian"


#	

def RunSKAT(sFile):

	
	
	os.system("Rscript /mnt/towel/BRCA/Code/R/SKAT_O_Test.R -G "+sOutputDir+"/GroupID.txt -I "+sFile+" -O "+sOutputDir+"/Result")

	








if __name__=="__main__":
	StartTime=(time.ctime())

	number_of_cpus=40
	manager=Manager()

	p = Pool(number_of_cpus)
	sFilelist=glob.glob(sOutputDir+"/Data/*.txt")
	Pool_map = p.map(RunSKAT, sFilelist)
	

	print(StartTime)
	print((time.ctime()))
	
