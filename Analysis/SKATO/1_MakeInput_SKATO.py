#!/usr/bin/env python
import gzip
import numpy as np
#import scipy.stats as stats
from scipy.stats import mannwhitneyu
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




def SKATOInput(sCaseInputFile,sControlInput,sOutputDir,sCaseIDFile,sControlIDFile):
	dExonicVariantDict={'splicing':1,'Frame_Shift_Del':1,'In_Frame_Ins':1,'Missense_Mutation':1,'Nonsense_Mutation':1,'Nonsense_Mutation':1,'In_Frame_Del':1,'Frame_Shift_Ins':1,'Frame_Shift_Substitution':1,'In_Frame_Shift_Substitution':1}
	fpCase=open(sCaseInputFile,"r")
	fpControl=open(sControlInput,"r")
	
	try:
		os.mkdir(sOutputDir)
	except:
		pass
		
	
	try:
		os.mkdir(sOutputDir+"/Data")
	except:
		pass
	
	
	try:
		os.mkdir(sOutputDir+"/Result")
	except:
		pass
	
	dCaseDict={}
	dControlDict={}
	
	lCaseTargetID=[]
	lControlTargetID=[]
	fCaseTargetIDFile=open(sCaseIDFile)
	fControlTargetIDFile=open(sControlIDFile)
	
	
	
	
	
	
	for sLine in fCaseTargetIDFile.readlines():
		sLine=sLine.strip()
		lCaseTargetID.append(sLine)
		
	
	for sLine in fControlTargetIDFile.readlines():
		sLine=sLine.strip()
		lControlTargetID.append(sLine)
	
	ParseCaseDataDict(dCaseDict, fpCase,lCaseTargetID)
	ParseCaseDataDict(dControlDict,fpControl, lControlTargetID)
	
	
	
	
	
	lGenelist=list(set(dCaseDict.keys())|set(dControlDict.keys()))
	
	
	fGroup=open(sOutputDir+"/GroupID.txt","w")
	
	for i in range(len(lCaseTargetID)):
		fGroup.write("{0}\t1\n".format(lCaseTargetID[i]))
	
	for i in range(len(lControlTargetID)):
		fGroup.write("{0}\t0\n".format(lControlTargetID[i]))
	
	for sGene in lGenelist:
		
		
		fDataout=open(sOutputDir+"/Data/"+sGene+".txt","w")
		
		try:
			lCaseKeylist=list(dCaseDict[sGene].keys())
		except:
			lCaseKeylist=[]
		try:
			lControlKeylist=list(dControlDict[sGene].keys())
		except:
			lControlKeylist=[]
		lKeylist=lCaseKeylist+lControlKeylist
		lUnifiedID=lCaseTargetID+lControlTargetID
		
		if not sGene in dCaseDict.keys():
			dCaseDict[sGene]={}
		
		if not sGene in dControlDict.keys():
			dControlDict[sGene]={}
		
		lUniqueKey=list(set(lKeylist))
		
		fDataout.write("{0}\n".format("\t".join(lUnifiedID)))
		for sKey in lUniqueKey:
			
			lWriteGenotypelist=[]
			if sKey in dCaseDict[sGene].keys():
				dKeyCasedict=dCaseDict[sGene][sKey]
			else:
				dKeyCasedict={}
			
			if sKey in dControlDict[sGene].keys():
				dKeyControldict=dControlDict[sGene][sKey]
			else:
				dKeyControldict={}
				
			for sID in lUnifiedID:
				if sID in dKeyCasedict:
					lWriteGenotypelist.append(dKeyCasedict[sID])
				elif sID in dKeyControldict:
					lWriteGenotypelist.append(dKeyControldict[sID])
				else:
					lWriteGenotypelist.append(0)
			
			
			if not 1 in lWriteGenotypelist:
				pass
			else:
				fDataout.write("{1}\t{0}\n".format("\t".join(map(str,lWriteGenotypelist)),sKey))
			
			
			
			
		
		

	



if __name__=="__main__":
	StartTime=(time.ctime())

	
	SKATOInput("/mnt/towel/BRCA/FinalFreeze/Cated_Assemble_Cluster.txt","/mnt/towel/BRCA/FinalFreeze/Nonsynonymous_re_TCGA.txt","/mnt/towel/BRCA/VCF/Analysis/SKATO/SNUvsTCGAre_withoutAsian","/mnt/towel/BRCA/FinalFreeze/ID/Final_SNUID.txt","/mnt/towel/BRCA/FinalFreeze/ID/TCGA_ID_withoutAsian.txt")
	
	
	
	
