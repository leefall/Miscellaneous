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


def Multipletest(sDir,sOutputDir):
	try:
		os.mkdir(sOutputDir)
	except:
		pass
	sHeader="SNUvsTCGA"
	os.system("cat "+sDir+"/Result/Result_*.txt > "+sOutputDir+"/"+sHeader+"_SKATO.txt")
	os.system("Rscript /mnt/towel/Ophthalmology/Code/R/Statistics/MultipletestCorrectionRare.R -I "+sOutputDir+"/"+sHeader+"_SKATO.txt -O "+sOutputDir)







if __name__=="__main__":
	StartTime=(time.ctime())

	number_of_cpus=30
	manager=Manager()

	p = Pool(number_of_cpus)
	sFilelist=glob.glob(sOutputDir+"/Data/*.txt")
#	

	
	
	Multipletest(sOutputDir,sCorrectionDir)
	print(StartTime)
	print((time.ctime()))
	
	
