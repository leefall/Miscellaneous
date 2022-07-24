#!/usr/bin/env python
import gzip
import numpy as np
import scipy.stats as stats
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




def AssembleData(sCase,sControl,sCaseID,sControlID,sOutputDir):
	
	fpCase=open(sCase,"r")
	fpControl=open(sControl,"r")
	try:
		os.mkdir(sOutputDir)
	except:
		pass
	fout=open(sOutputDir+"/For_Fisher_Input.txt","w")
	
	dExonicVariantDict={'Frame_Shift_Del':1,'In_Frame_Ins':1,'Missense_Mutation':1,'Nonsense_Mutation':1,'Nonsense_Mutation':1,'In_Frame_Del':1,'Frame_Shift_Ins':1,'Frame_Shift_Substitution':1,'In_Frame_Shift_Substitution':1}
	dExonicVariantDict={'Frame_Shift_Del':1,'In_Frame_Ins':1,'Missense_Mutation':1,'Nonsense_Mutation':1,'Nonsense_Mutation':1,'In_Frame_Del':1,'Frame_Shift_Ins':1,'Frame_Shift_Substitution':1,'In_Frame_Shift_Substitution':1,"splicing":1}
	dCaseDict={}
	dControlDict={}
	
	
	
	fTargetIDFile=open(sCaseID)
	dCaseSampleDict=dict()
	
	
	for sLine in fTargetIDFile.readlines():
		sLine=sLine.strip()
		dCaseSampleDict[sLine]=1
	
	
	fControlIDFile=open(sControlID)
	dControlSampleDict=dict()
	
	
	for sLine in fControlIDFile.readlines():
		sLine=sLine.strip()
		dControlSampleDict[sLine]=1
	
	
	
	nCaseSample=len(dCaseSampleDict.keys())
	nControlSample=len(dControlSampleDict.keys())
	fpCase.readline()
	
	lControlIDs=dControlSampleDict.keys()
	
	
	for sLine in fpCase:
		sLine=sLine.strip()
		t=sLine.split("\t")
		
		(sChr,sPosition,sRef,sAlt,sVariantType,sGene,sEffect,sIDs)=(t[0],t[1],t[2],t[3],t[6],t[4],t[5],t[-1])
	
		
		if sEffect in dExonicVariantDict:
			sKey="_".join([sGene,sChr,sPosition,sRef])
			#print(sKey)
			lIDS=sIDs.split(",")
			for sID in lIDS:
				if sKey in dCaseDict:
					dCaseDict[sKey]+=1
				else:
					dCaseDict[sKey]=1
	fpControl.readline()
	for sLine in fpControl:
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,sPosition,sRef,sAlt,sVariantType,sGene)=(t[0],t[1],t[6],t[7],t[4],t[-1])
		#try:
		sChr=sChr.replace("chr","")
		sTCGAID=t[-3]
		
		if sTCGAID in dControlSampleDict:
		
			sKey="_".join([sGene,sChr,sPosition,sRef])
			if sKey in dControlDict:
				dControlDict[sKey]+=1
			else:
				dControlDict[sKey]=1
			
		
	
	
	lAllKeys=list(set(dCaseDict.keys())|set(dControlDict.keys()))
	
	
	nUsedCase=0
	for sKey in lAllKeys:
		sKeyWrite=1
		if not sKey in dCaseDict.keys():
			nCaseMut=0
			nCaseNonMut=nCaseSample
			
		else:
		
			nCaseMut=dCaseDict[sKey]
			nCaseNonMut=nCaseSample-int(nCaseMut)
		if not sKey in dControlDict.keys():
			nControlMut=0
			nControlNonmut=nControlSample
		else:
			nControlMut=dControlDict[sKey]
			nControlNonmut=nControlSample-int(nControlMut)
		if 1:
			
			
			fout.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(sKey,str(nCaseMut),str(nCaseNonMut),str(nControlMut),str(nControlNonmut)))

	
	





def Fisherstest(sDir):
	
	
	
	
	os.chdir(sDir)


	
	fp=open("./For_Fisher_Input.txt","r")
	fout=open("Result_Fisher.txt","w")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sGene,sChr,nPosition,sRef)=(t[0].split("_"))
		nOddsratio,nPvalue = stats.fisher_exact([[int(float(t[1])), int(float(t[2]))], [int(float(t[3])), int(float(t[4]))]])
		if sChr=="X":
			sChr="23"
		fout.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(sGene,sChr,nPosition,sRef,nPvalue,nOddsratio,",".join([str(float(t[1])), str(float(t[2])), str(float(t[3])), str(float(t[4]))])))
		



def Multipletest(sDir):
	sHeader="DEP"
	
	
	
	os.system("Rscript /mnt/towel/BRCA/Code/R/MultipletestCorrectionFisher.R -I "+sDir+"/Result_Fisher.txt -O "+sDir)

	



def consumer_task(q, cosmic_dict):
	while not q.empty():
		value=q.get(True, 0.05)
		
		
		
		
		
		
		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))


def producer_task(q, cosmic_dict):

	lChrlist=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
	for i in lChrlist:
		value=i
		cosmic_dict[value]=None
		logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
		q.put(value)





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()

	number_of_cpus=30


	#Case : TNBC_NonSynonymous_SNU.maf
	#Control : TNBC_TCGA_Vanderbilt.maf
	sCase="/mnt/towel/BRCA/FinalFreeze/StrictFilter_Fudan_SomaticMutation.txt"
	sControl="/storage/home/leefall2/mypro/METABRIC/Intersected_Realsorted_data_mutations_mskcc.maf"
	sCaseID="/mnt/towel/BRCA/FinalFreeze/ID/FUDAN_ID.txt"
	sControlID="/storage/home/leefall2/mypro/METABRIC/TNBC_mut_ID.txt"
	sOutputDir="/mnt/towel/BRCA/VCF/Analysis/Fisher/FUDAN_METABRIC"
	AssembleData(sCase,sControl,sCaseID,sControlID,sOutputDir)
	
	
	
	Fisherstest(sOutputDir)
	Multipletest(sOutputDir)
	
	
	
	
