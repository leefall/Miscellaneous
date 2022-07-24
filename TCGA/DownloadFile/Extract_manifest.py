#!/usr/bin/env python
dTNBCID=dict()
dNonTNBCFileID=dict()

dCheckID=dict()


fClinic=open("/mnt/alpha/leefall2/TCGA_STAD/NoNASTAD.txt")
fClinic.readline()
for sLine in fClinic.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	dTNBCID[t[0]]=0

print(len(dTNBCID))

fIDFile=open("/mnt/alpha/leefall2/TCGA_STAD/STAD_TCGAID_table.txt")

fIDFile.readline()
n=0
for sLine in fIDFile.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	sSampleID=t[-1]
	sTCGAID=sSampleID[0:12]
	if sTCGAID in dTNBCID:
		TissueCode=sSampleID[13:15]
		if TissueCode=="06":
			pass
		elif TissueCode=="10":
			#dNonTNBCFileID[t[0]]=0
			if sTCGAID in dCheckID.keys():
				dCheckID[sTCGAID]["Normal"]="Yes"
				dTNBCID[sTCGAID]+=1
			else:
				dCheckID[sTCGAID]={"Normal":"No","Tumor":"No"}
				dCheckID[sTCGAID]["Normal"]="Yes"
				dTNBCID[sTCGAID]+=1
		elif TissueCode=="11":
			if sTCGAID in dCheckID.keys():
				dCheckID[sTCGAID]["Normal"]="Yes"
				dTNBCID[sTCGAID]+=1
			else:
				dCheckID[sTCGAID]={"Normal":"No","Tumor":"No"}
				dCheckID[sTCGAID]["Normal"]="Yes"
				dTNBCID[sTCGAID]+=1
		elif TissueCode=="01":
			if sTCGAID in dCheckID.keys():
				dCheckID[sTCGAID]["Tumor"]="Yes"
				dTNBCID[sTCGAID]+=1
			else:
				dCheckID[sTCGAID]={"Normal":"No","Tumor":"No"}
				dCheckID[sTCGAID]["Tumor"]="Yes"
				dTNBCID[sTCGAID]+=1
			
			
fIDFile.close()
fIDFile=open("/mnt/alpha/leefall2/TCGA_STAD/STAD_TCGAID_table.txt")
fIDFile.readline()
for sLine in fIDFile.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	sSampleID=t[-1]
	sTCGAID=sSampleID[0:12]
	if sTCGAID in dCheckID:
		if ((dCheckID[sTCGAID]["Tumor"]=="Yes") & (dCheckID[sTCGAID]["Normal"]=="Yes")):
			dNonTNBCFileID[t[0]]=sTCGAID
			
			
			
		
	
	

			
		#print(sLine)
	n+=1
#print(str(n))
#print(len(dTNBCID))
#print(len(dNonTNBCFileID))

print(len(dNonTNBCFileID))
fFile=open("/mnt/alpha/leefall2/TCGA_STAD/gdc_manifest_20220321_064652.txt")
foutFile=open("/mnt/alpha/leefall2/TCGA_STAD/NonNA_STAD_gdc_manifest_20220321_064652.txt","w")
foutFile.write(fFile.readline())
for sLine in fFile.readlines():
	t=sLine.split("\t")
	if t[0] in dNonTNBCFileID:
		foutFile.write(sLine)
	
fFile.close()
for sID in dTNBCID.keys():
	if dTNBCID[sID]!=2:
		print(sID)
		print(str(dTNBCID[sID]))


fFile=open("/mnt/alpha/leefall2/TCGA_STAD/gdc_manifest_20220321_064652.txt")
foutFile=open("/mnt/alpha/leefall2/TCGA_STAD/Paired_NonNA_STAD_gdc_manifest_20220321_064652.txt","w")
foutFile.write(fFile.readline())
for sLine in fFile.readlines():
	t=sLine.split("\t")
	if t[0] in dNonTNBCFileID:
		if dTNBCID[dNonTNBCFileID[t[0]]]==2:
			foutFile.write(sLine)
	








