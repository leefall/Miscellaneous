#!/usr/bin/env python
import glob
import os
import sys

def getSampleName():
	dUUID_Sample=dict()
	dNormalDict=dict()
	dSoliddict=dict()
	#fp=open("gdc_manifest_20190921_124302.txt","r")
	#fp=open("TCGAID_table.txt","r")
	#fp=open("/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA/Curated_Pairedset.txt")
	fp=open("STAD_TCGAID_table.txt")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sUUID, dSample)=(t[0],t[1])
		#dUUID_Sample[sUUID]=sTCGAID
		sTCGAID=dSample[0:12]
		sSampleStatus=dSample[13:15]
			
		if not sSampleStatus=="01":
			if sTCGAID in dNormalDict.keys():
				dNormalDict[sTCGAID][sSampleStatus]=1
			else:
				dNormalDict[sTCGAID]={"10":0,"11":0}
				dNormalDict[sTCGAID][sSampleStatus]=1
	fp.close()
	
	for sTCGAID in dNormalDict.keys():
		if ((dNormalDict[sTCGAID]["10"]==1) and (dNormalDict[sTCGAID]["11"]==1)):
			#dBlooddict[sTCGAID]=0
			pass
		elif ((dNormalDict[sTCGAID]["10"]==0) and (dNormalDict[sTCGAID]["11"]==1)):
			dSoliddict[sTCGAID]=0
		elif ((dNormalDict[sTCGAID]["10"]==1) and (dNormalDict[sTCGAID]["11"]==0)):
			#dBlooddict[sTCGAID]=0
			pass
		else:
			print(sTCGAID)
	
	#print(dSoliddict)
	#print(len(dSoliddict))
	
#	fp=open("/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA/Curated_Pairedset.txt")
	fp=open("STAD_TCGAID_table.txt")	
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sUUID, dSample)=(t[0],t[1])
		#dUUID_Sample[sUUID]=sTCGAID
		sTCGAID=dSample[0:12]
		sSampleStatus=dSample[13:15]
		
		
		#if not sSampleStatus=="01":
		if sSampleStatus=="11":
			if sTCGAID in dSoliddict.keys():
					#dNormalDict[sTCGAID][sSampleStatus]=1
				dUUID_Sample[sUUID]=dSample
			else:
				pass
		else:
			dUUID_Sample[sUUID]=dSample


		
		
		
		
		
		
		
		
		#if sSampleStatus=="01":
#			sStatus="T"
		#else:
#			sStatus="N"
		
		
		
		
		
	return dUUID_Sample

def MakeSymbol(sDir,dUUID_Sample):
	lBam=glob.glob("*.bam*")
	for i in lBam:
		if not ".bai" in i:
			sBam=i
	dSample=dUUID_Sample[sDir]
	#print(sBam, dSample)
	sTCGAID=dSample[0:12]
	sSampleStatus=dSample[13:15]
	if sSampleStatus=="01":
		sStatus="T"
		sCurrent=os.getcwd()
		os.system("ln -s "+sCurrent+"/"+sBam+" ../Symbol/"+sTCGAID+"-"+sStatus+".bam")
	elif sSampleStatus=="06":
		#sStatus="N"
		print("Metastasis")
		pass
	elif sSampleStatus=="10":
		sStatus="N"
		sCurrent=os.getcwd()
		os.system("ln -s "+sCurrent+"/"+sBam+" ../Symbol/"+sTCGAID+"-"+sStatus+".bam")
	elif sSampleStatus=="11":
		sStatus="N"
		sCurrent=os.getcwd()
		os.system("ln -s "+sCurrent+"/"+sBam+" ../Symbol/"+sTCGAID+"-"+sStatus+".bam")
	else:
		print(sDir)
		print(sTCGAID)
		print("Not Paired Sample")
#		sys.exit()
	#if sTCGAID in dIDDict.keys():
#		dIDDict[sTCGAID]+=1
#	else:
		#dIDDict[sTCGAID]=1
		
	
	
	
	
	
	#if sTCGAID in sNonTNBCBRCA:
#		fout.write(sTCGAID)
		#fout.write("\t")
		#fout.write(dSample)
		#fout.write("\t")
		#fout.write(sDir)
		#fout.write("\n")
	
	
	sCurrent=os.getcwd()
	#os.system("ln -s "+sCurrent+"/"+sBam+" ../"+dSample+".bam")
	
	return dIDDict

if __name__=="__main__":
	os.chdir("/mnt/backup/TCGA_STAD/Paired_mc3")
	try:
		os.mkdir("Symbol")
	except:
		pass
	
	#lDirlist=os.listdir('.')
	#print(lDirlist)
	
	#for top, dirs, files in os.walk("./"):
#		print(top, dirs, files)
	
	
	dUUID_Sample=getSampleName()
	#print(dUUID_Sample)
	#fp=open("/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA/Clinical/NonNA_nonTNBC_BRCA_Cell_Curated_Outcome.txt")
	#fp.readline()
	#sNonTNBCBRCA=set()
	#for sLine in fp.readlines():
#		sBRCAID=sLine.split("\t")[0]
		#sNonTNBCBRCA.add(sBRCAID)
	
	dIDDict=dict()
	
	lDirlist=next(os.walk('.'))[1]
	#fout=open("/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA/FileforProcess.txt","w")
	for sDir in lDirlist:
		#print(sDir)
		
		if sDir in dUUID_Sample:
			#print(sDir)
			
		
		
			os.chdir(sDir)
			MakeSymbol(sDir,dUUID_Sample)
			os.chdir("..")
			dUUID_Sample[sDir]=1
	
	#for sID in dUUID_Sample.keys():
#		if dUUID_Sample[sID]!=1:
			#print(sID)
		
		
		
	



