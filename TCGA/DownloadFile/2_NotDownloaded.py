#!/usr/bin/env python
import sys
dUUIDSampleID=dict()

sSampleset=set()
sNormalSampleset=set()


dNormalSampleBarcode=dict()
dTumorSampleBarcode=dict()


#fp=open("/mnt/alpha/leefall2/TCGA_OV/mc3/allmc3NormalSampleBarcode.txt")
fp=open("/mnt/alpha/leefall2/TCGA_LGG/mc3/allmc3NormalSamplealiquotBarcode.txt")
for sLine in fp.readlines():
	sLine=sLine.strip()
	dNormalSampleBarcode[sLine]=0
fp.close()

fp=open("/mnt/alpha/leefall2/TCGA_LGG/mc3/allmc3SamplealiquotBarcode.txt")

for sLine in fp.readlines():
	sLine=sLine.strip()
	dTumorSampleBarcode[sLine]=0
fp.close()


#print(len(dNonTNBCFileID))
fFile=open("/mnt/alpha/leefall2/TCGA_LGG/All_NonNAmc3_gdc_manifest_20211213_044749.txt")
#foutFile=open("/mnt/alpha/leefall2/TCGA_HNSC/NonNA_gdc_manifest_20211116_094540.txt","w")
fPaired=open("/mnt/alpha/leefall2/TCGA_LGG/rawbam/Paired_NonNA_LGG_gdc_manifest.txt")
dPaired=dict()
fPaired.readline()
for sLine in fPaired.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	dPaired[t[0]]=1
	
	
fTCGAID=open("/mnt/alpha/leefall2/TCGA_LGG/LGG_TCGAaliquotID_table.txt","r")
fTCGAID.readline()
for sLine in fTCGAID.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	dUUIDSampleID[t[0]]=t[-1]


fout=open("/mnt/alpha/leefall2/TCGA_LGG/NotPaired_NonNA_gdc_manifest_20211116_094540.txt","w")
sFirst=fFile.readline().strip()
fout.write(sFirst)
fout.write("ID\n")
fExtraPaired=open("/mnt/alpha/leefall2/TCGA_LGG/Paired_mc3_NonNA_gdc_manifest_20211116_094540.txt","w")
fExtra=open("/mnt/alpha/leefall2/TCGA_LGG/ForReview_mc3_NonNA_gdc_manifest_20211116_094540.txt","w")
fExtra.write(sFirst)
fExtra.write("ID\n")
dTumor=dict()
dNormal=dict()
for sLine in fFile.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	if not t[0] in dPaired:
		if t[0] in dUUIDSampleID.keys():
			fout.write(sLine)
			fout.write("\t")
			fout.write(dUUIDSampleID[t[0]])
			fout.write("\t")
			fout.write(dUUIDSampleID[t[0]][0:12])
			fout.write("\n")
			sSampleBarcode=dUUIDSampleID[t[0]]
			sTCGAID=sSampleBarcode[0:12]
			#print(sSampleBarcode)
			if sSampleBarcode in dTumorSampleBarcode:
				if sTCGAID in dTumor.keys():
					dTumor[sTCGAID]["Tumor"]+=1
				else:
					dTumor[sTCGAID]=dict()
					dTumor[sTCGAID]["Tumor"]=1
#				fExtra.write(sLine)
#				fExtra.write("\t")
#				fExtra.write(dUUIDSampleID[t[0]])
#				fExtra.write("\t")
#				fExtra.write(dUUIDSampleID[t[0]][0:12])
#				fExtra.write("\n")
			elif sSampleBarcode in dNormalSampleBarcode:
				if sTCGAID in dNormal.keys():
					dNormal[sTCGAID]["Normal"]+=1
				else:
					dNormal[sTCGAID]=dict()
					dNormal[sTCGAID]["Normal"]=1
#				fExtra.write(sLine)
#				fExtra.write("\t")
#				fExtra.write(dUUIDSampleID[t[0]])
#				fExtra.write("\t")
#				fExtra.write(dUUIDSampleID[t[0]][0:12])
#				fExtra.write("\n")
#			
			
			
			
		else:
			print("NotUUID sample")
			print(sLine)
			sys.exit()
		

fFile.close()

fFile=open("/mnt/alpha/leefall2/TCGA_LGG/All_NonNAmc3_gdc_manifest_20211213_044749.txt")
fExtraPaired.write(fFile.readline())
#sFirst=fFile.readline().strip()

for sLine in fFile.readlines():
	sLine=sLine.strip()
	t=sLine.split("\t")
	if not t[0] in dPaired:
		if t[0] in dUUIDSampleID.keys():
			sSampleBarcode=dUUIDSampleID[t[0]]
			sTCGAID=sSampleBarcode[0:12]
			if sSampleBarcode in dTumorSampleBarcode:
				if ((dTumor[sTCGAID]["Tumor"]==1) and (dNormal[sTCGAID]["Normal"]==1)):
					fExtraPaired.write(sLine)
					fExtraPaired.write(dUUIDSampleID[t[0]])
					fExtraPaired.write("\t")
					fExtraPaired.write(dUUIDSampleID[t[0]][0:12])
					fExtraPaired.write("\n")
				else:
					fExtra.write(sLine)
					fExtra.write("\t")
					fExtra.write(dUUIDSampleID[t[0]])
					fExtra.write("\t")
					fExtra.write(dUUIDSampleID[t[0]][0:12])
					fExtra.write("\n")
					
					
			elif sSampleBarcode in dNormalSampleBarcode:
				if ((dTumor[sTCGAID]["Tumor"]==1) and (dNormal[sTCGAID]["Normal"]==1)):
					fExtraPaired.write(sLine)
					fExtraPaired.write(dUUIDSampleID[t[0]])
					fExtraPaired.write("\t")
					fExtraPaired.write(dUUIDSampleID[t[0]][0:12])
					fExtraPaired.write("\n")
				else:
					fExtra.write(sLine)
					fExtra.write("\t")
					fExtra.write(dUUIDSampleID[t[0]])
					fExtra.write("\t")
					fExtra.write(dUUIDSampleID[t[0]][0:12])
					fExtra.write("\n")
					
					
				
				
				
#				if sSampleBarcode in dExtra.keys():
#					dExtra[sSampleBarcode]["Normal"]=1
#				else:
#					dExtra[sSampleBarcode]=dict()
#					dExtra[sSampleBarcode]["Normal"]=1




#foutFile.write(fFile.readline())
#for sLine in fFile.readlines():
#	t=sLine.split("\t")
#	if t[0] in dNonTNBCFileID:
#		foutFile.write(sLine)
#	
#fFile.close()
#for sID in dTNBCID.keys():
#	if dTNBCID[sID]!=2:
#		print(sID)
#		print(str(dTNBCID[sID]))
#
##dUUIDSampleID[sUUID]=sTCGAID
##fFile=open("/mnt/alpha/leefall2/TCGA_UCEC/gdc_manifest_20211028_130109.txt")
#
#			
			
			
			
			#foutFile.write(sLine)
			
			
			
	








