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
	
	
def gVCFtoVCF(sFile):
	sID=sFile.split("/")[-1]
	sID=sID.split(".")[0]
	
	fp=gzip.open(sFile,"rt")
	
	fout=open("/mnt/alpha/leefall2/UKBB/vcf/"+sID+".vcf","w")
	
	dNotGT={"0/0":1,"0/.":1,"./0":1,"./.":1}
	
	for sLine in fp.readlines():
		if sLine[0]=="#":
			#if "contig" in sLine:
#				pass
	#		else:
			fout.write(sLine)
		else:
			t=sLine.split("\t")
			sGene=t[-1]
			sGenotype=sGene.split(":")[0]
			if sGenotype in dNotGT:
				pass
			else:
				t=sLine.split("\t")
				t[4]=t[4].replace(",<*>","")
				fout.write("{0}".format("\t".join(t)))
	fout.close()
	subprocess.call('''
	bgzip /mnt/alpha/leefall2/UKBB/vcf/'''+sID+'''.vcf
	tabix -p vcf /mnt/alpha/leefall2/UKBB/vcf/'''+sID+'''.vcf.gz
	
	
	''',shell=True)
	
	
#
#	subprocess.call('''
#	
#	
#	gzip -dc '''+sFile+'''| awk '{if($5!="<NON_REF>"||$) print}' > /mnt/towel/leefall2/UKBiobank/vcf/'''+sID+'''.vcf
#	''',shell=True)
	
#	subprocess.call('''
#	bgzip /mnt/towel/leefall2/UKBiobank/vcf/'''+sID+'''.vcf
#	tabix -p vcf /mnt/towel/leefall2/UKBiobank/vcf/'''+sID+'''.vcf.gz
#	
#	
#	''',shell=True)

if __name__=="__main__":
	StartTime=(time.ctime())

	number_of_cpus=30
	manager=Manager()

	p = Pool(number_of_cpus)
	#sFilelist=glob.glob("/mnt/towel/UKBiobank/FE_vcfs/*.gvcf.gz")#4011848
	sFilelist1=glob.glob("/mnt/Hercules/UKBB/UKB_200K_WES_gVCFs/1/*.gz")
	sFilelist2=glob.glob("/mnt/Hercules/UKBB/UKB_200K_WES_gVCFs/2/*.gz")
	sFilelist3=glob.glob("/mnt/Hercules/UKBB/UKB_200K_WES_gVCFs/3/*.gz")
	sFilelist4=glob.glob("/mnt/Hercules/UKBB/UKB_200K_WES_gVCFs/4/*.gz")
	sFilelist5=glob.glob("/mnt/Hercules/UKBB/UKB_200K_WES_gVCFs/5/*.gz")
	sFilelist6=glob.glob("/mnt/Hercules/UKBB/UKB_200K_WES_gVCFs/6/*.gz")
	sFinalFilelist=sFilelist1+sFilelist2+sFilelist3+sFilelist4+sFilelist5+sFilelist6
	
	
	
	Pool_map = p.map(gVCFtoVCF, sFinalFilelist)
	
	
	
	
	
	
	

	print(StartTime)
	print((time.ctime()))
	

