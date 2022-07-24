#!/usr/bin/env python
import glob
import os






l50000vcflist=glob.glob("/mnt/towel/leefall2/UKBiobank/vcf/*.vcf.gz")
dVCFdict=dict()
sVCFset=set()
lHeaderlist=[]
sHeaderset=set()
for sFile in l50000vcflist:
	sName=sFile.split(".")[0]
	sName=sName.split("/")[-1]
	dVCFdict[sName]=0
	sVCFset.add(sName)


#print(len(dVCFdict))


fp=open("/mnt/towel/UKBiobank/Clinical/39581/HES/hesin_diag.txt","r")
fout=open("/mnt/towel/leefall2/UKBiobank/Clinical/Intersected_hesin_diag.txt","w")


sHeaderLine=fp.readline()

fout.write(sHeaderLine)
n=0



for sLine in fp.readlines():
	t=sLine.split("\t")
	if t[0] in dVCFdict:
		fout.write(sLine)


#for sHeader in sHeaderLine.split(","):
#	sHeader=sHeader.split("-")[0]
#	sHeader=sHeader.replace('"','')
##	print(sHeader)
#	if sHeader in dVCFdict:
#		lHeaderlist.append(n)
#	else:
#		pass
#	sHeaderset.add(sHeader)
#	n+=1




#print(lHeaderlist)
#print(dVCFdict)
#print(len(sHeaderset))
#print(len(sVCFset))

#print(sVCFset&sHeaderset)



