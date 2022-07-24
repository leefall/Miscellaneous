#!/usr/bin/env python
import gzip
import numpy as np
import os, subprocess, glob, time, tabix
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager


logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter=logging.Formatter("%(asctime)s - %(message)s")

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)


dConversionID=dict()

#fConverse=open("/mnt/towel/Ophthalmology/Code/ParseVariant/AnalysisReadySomatic/New_LabelID.txt")
fConverse=open("/mnt/towel/Ophthalmology/Code/ParseVariant/AnalysisReadySomatic/New_LabelID_PrimEnucMeta_Primary11.txt")
fConverse.readline()
for sLine in fConverse.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	dConversionID[t[-1]]=t[0]





def MAFWriter(sFile):
	
	lSynonymous=["upstream","downstream","intronic","synonymous_SNV","unknown","upstream;downstream"]
	sOut=sFile.split("/")[-1]
	fou=open("/mnt/towel/BRCA/VCF/Analysis/Signature/TNBC_SNU"+"/Signature_"+sOut.replace(".txt",".maf"),"w")
	
	
	fou.write("Hugo\tEntrez\tCenter\tGenome\tChrom\tStart\tEnd\tStrand\tClassification\tVariant_Type\tRef\tAlt1\tAlt2\tdbSNP\tSNP_Val_status\tTumor_Sample\n")


	fp=open(sFile)
	fp.readline()
	#for sKey in dVariantDict.keys():
	for sLine in fp.readlines():
		#print sKey
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[1],t[2],t[7],t[8])
		
		sPassKey=0
		
		
		sID=t[-2]
				
		sRealMolecularEffects=t[4]
		
		
		if t[5] =="SNV":
			sVariantType="SNP"
		elif t[5] =="Deletion":
			sVariantType="DEL"
		elif t[5] =="Insertion":
			sVariantType="INS"
		
		fou.write("{0}\n".format(\
					"\t".join(map(str,[t[0],".","PMSNH","GRCh37",sChr,str(t[2]),str(t[3]),"1",sRealMolecularEffects,sVariantType,sRef,".",sAlt,".",".",sID]))))
			
			
				
		
if __name__=="__main__":
	StartTime=(time.ctime())

	
	dVariantDict=dict()

	MAFWriter("/mnt/towel/BRCA/FinalFreeze/TNBC_All_SNU.txt")
	
	print(StartTime)
	print(time.ctime())	
	
	
	
