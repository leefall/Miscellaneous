#!/usr/bin/env python
sManifest="Paired_NonNA_LIHC_gdc_manifest_20220205_062511.txt"
fp=open(sManifest)
lManifest=[]
sHeader=fp.readline()

n=0
nTotal=1
for sLine in fp.readlines():
	if n==0:
		fout=open(str(nTotal)+"_"+str(nTotal+40)+"_"+sManifest,"w")
		lManifest.append(str(nTotal)+"_"+str(nTotal+40)+"_"+sManifest)
		fout.write(sHeader)
		fout.write(sLine)
		n+=1
		nTotal+=1
	elif n==40:
		fout.write(sLine)
		fout.close()
		n=0
		nTotal+=1
	else:
		fout.write(sLine)
		n+=1
		nTotal+=1


fManifest=open("Raw_Split_Download.sh","w")

for sFile in lManifest:
	fManifest.write("gdc-client download -t gdc-user-token.2022-01-23T02_46_40+09_00.txt -m "+sFile+" -d /mnt/alpha/leefall2/TCGA_LIHC/rawbam/Paired \n")





