#!/usr/bin/env python
import gzip
import time
import sys
fSampleID=open("/mnt/alpha/leefall2/TCGA_LGG/mc3/allmc3SamplealiquotBarcode.txt","w")
fNormalSampleID=open("/mnt/alpha/leefall2/TCGA_LGG/mc3/allmc3NormalSamplealiquotBarcode.txt","w")
def UVMID_Parsing():
	
	fp=open("/mnt/alpha/leefall2/TCGA_LGG/NoNAClnicalBarcode.txt")
	sAllDownloadableFile=set()
	
	fp.readline()
	for sLine in fp.readlines():
		sLine=sLine.strip()
		#t=sLine.split("\t")
		#(sID, sER, sPR, sHER2)=(t[0],t[1],t[2],t[3])
		
#		if ((sER=="Negative") and (sPR=="Negative") and (sHER2=="Negative")):
#		sTCGAID=sLine
		sFilename=sLine.split("\t")
		sTCGAID=sFilename[0]
		#print(sTCGAID)
#		if sTCGAID[13:15]=="10":
		#sTCGAID=sTCGAID[0:12]
		sAllDownloadableFile.add(sTCGAID)
		
	return(sAllDownloadableFile)
	

StartTime=(time.ctime())


sUVMID=UVMID_Parsing()

#print(len(sUVMID))
fout=open("/mnt/alpha/leefall2/TCGA_LGG/mc3/MC3LGG.txt","w")
fID=open("/mnt/alpha/leefall2/TCGA_LGG/mc3/allmc3Barcode.txt","w")
fzip=gzip.open("/mnt/towel/TCGA/mc3/mc3.v0.2.8.PUBLIC.maf.gz","rt")
#fout.write(fzip.readline())
sBarcodeset=set()
sSampleset=set()
sNormalSampleset=set()
dNormalSample=dict()

for sLine in fzip.readlines():
	t=sLine.split("\t")
	sTCGATumorID=t[15]
	sTCGANormalID=t[16]
	sTCGAID=sTCGATumorID[0:12]
	sTCGASampleID=sTCGATumorID
	sTCGANormalSampleID=sTCGANormalID
	if sTCGAID[13:15]=="06":
		pass
	elif sTCGAID[13:15]=="02":
		pass
	else:
		if sTCGAID in sUVMID:
			fout.write(sLine)
			sBarcodeset.add(sTCGAID)
			sSampleset.add(sTCGASampleID)
			sNormalSampleset.add(sTCGANormalSampleID)
			if sTCGAID in dNormalSample.keys():
				if sTCGANormalSampleID[13:15]=="10":
					dNormalSample[sTCGAID]["Blood"]=sTCGANormalSampleID
				elif sTCGANormalSampleID[13:15]=="11":
					dNormalSample[sTCGAID]["Solid"]=sTCGANormalSampleID
			else:
				dNormalSample[sTCGAID]=dict()
				if sTCGANormalSampleID[13:15]=="10":
					dNormalSample[sTCGAID]["Blood"]=sTCGANormalSampleID
				elif sTCGANormalSampleID[13:15]=="11":
					dNormalSample[sTCGAID]["Solid"]=sTCGANormalSampleID
				elif sTCGANormalSampleID[13:15]=="12":
					dNormalSample[sTCGAID]["Solid"]=sTCGANormalSampleID
				else:
					print(sTCGAID)
					print(sTCGANormalSampleID)
					pass
					#sys.exit()
	#print(t)
	#break
for sID in sBarcodeset:
	fID.write("{0}\n".format(sID))

for sID in sSampleset:
	if sID[13:15]=="06":
		pass
	elif sID[13:15]=="02":
		pass
	else:
		fSampleID.write("{0}\n".format(sID))

for sTCGAID in dNormalSample:
	#if sID[13:15]=="06":
#		pass
	#else:
		#fNormalSampleID.write("{0}\n".format(sID))
	
	if len(dNormalSample[sTCGAID])==2:
		fNormalSampleID.write("{0}\n".format(dNormalSample[sTCGAID]["Blood"]))
	else:
		for i in dNormalSample[sTCGAID].keys():
			fNormalSampleID.write("{0}\n".format(dNormalSample[sTCGAID][i]))
		
		


print("Start Time")
print(StartTime)
print("End Time")
print(time.ctime())


