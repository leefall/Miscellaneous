#!/usr/bin/env python
import os
import gzip
import subprocess
import sys, time, random, re ,requests, logging, glob, tabix
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager
import multiprocessing

def Cosmic_Filter(value, queue):
	print(value)
	os.system("Rscript "+value)
	return
	


def ExtractID(sFile):
	fp=gzip.open(sFile,"rt")
	for sLine in fp.readlines():
		if sLine[0:2]=="#C":
			sLine=sLine.strip()
			sRealID=sLine.split("\t")[-1]
			return sRealID



def gVCFtoVCF(sFile):
	
	sID=sFile.split(".")[0]
	
	sRealID=sID.split("_")[0]
#	sRealID=ExtractID("/mnt/towel/leefall2/UKBiobank/VCF/"+sID+".vcf.gz")
	subprocess.call('''
	
	
	
	
	
	
	vt decompose '''+sID+'''.vcf.gz | /storage/home/leefall2/clara/ANNOVAR/annovar/convert2annovar.pl -format vcf4 --withzyg - > /mnt/alpha/leefall2/UKBB/ANNOVAR/'''+sRealID+'''.txt
	
	
	''',shell=True)

if __name__=="__main__":
	StartTime=(time.ctime())
	os.chdir("/mnt/alpha/leefall2/UKBB/vcf")
	#os.chdir("/storage/home/leefall2/clara/New_Protected_mutation/KIRC/Test_step/KIRC_CoxphfLOHMulti_PFI_R")
	sFilelist=glob.glob("*.vcf.gz")
	#sFilelist=sFilelist[0:1]
	number_of_cpus=30
	manager=Manager()

	p = Pool(number_of_cpus)
	
	Pool_map = p.map(gVCFtoVCF, sFilelist)
	
	
	print(StartTime)
	print((time.ctime()))
	

