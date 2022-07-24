#!/usr/bin/env python
import glob
import os
import sys

		
		
		
def MakeSymbol(sDir):
	lBam=glob.glob("*.xml")
	
	for i in lBam:
		
		sBam=i
		sCurrent=os.getcwd()
		sUUID=sCurrent.split("/")[-2]
		sTail=sBam.split("_")[-1]
		#MetaSymbol
		#os.system("ln -s "+sCurrent+"/"+sBam+" ../Paired/"+sTCGAID+"-"+sStatus+".bam")
		#print("ln -s "+sCurrent+"/"+sBam+" ../../../MetaSymbol/"+sUUID+"_"+sTail)
		os.system("ln -s "+sCurrent+"/"+sBam+" ../../../MetaSymbol/"+sUUID+"_"+sTail)

#	else:
		#dIDDict[sTCGAID]=1
		
	
	
	
	
	
	#if sTCGAID in sNonTNBCBRCA:
#		fout.write(sTCGAID)
		#fout.write("\t")
		#fout.write(dSample)
		#fout.write("\t")
		#fout.write(sDir)
		#fout.write("\n")
	
	
	#sCurrent=os.getcwd()
	#os.system("ln -s "+sCurrent+"/"+sBam+" ../"+dSample+".bam")
	
	return dIDDict

if __name__=="__main__":
	os.chdir("/mnt/alpha/leefall2/TCGA_BLCA/Metadata")

	
	#lDirlist=os.listdir('.')
	#print(lDirlist)
	
	#for top, dirs, files in os.walk("./"):
#		print(top, dirs, files)
	
	
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
		
		os.chdir(sDir)
		lZipFiles=glob.glob("*.tar.gz")
		for sFile in lZipFiles:
			os.system("tar -xzvf "+sFile)
		lSubDirlist=next(os.walk('.'))[1]
		
		for sSubDir in lSubDirlist:
			os.chdir(sSubDir)
			MakeSymbol(sSubDir)
			os.chdir("..")
			
		os.chdir("..")
	
	
	
#for sID in dUUID_Sample.keys():
#		if dUUID_Sample[sID]!=1:
			#print(sID)
		
		
		
	



