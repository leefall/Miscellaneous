#!/usr/bin/env python
import glob
import os
import sys



def sForANNOVAR(sFile):
	
	
	
	
	fp=open(sFile)
	fout=open("/mnt/QNAP/leefall2/Korea1K/Decompose/Decompose_"+sFile,"w")
	
	for sLine in fp.readlines():
		if sLine[0]=="#":
			pass
		else:
			sLine=sLine.strip()
			t=sLine.split("\t")
			(sChr, nPosition, sRef, sAlt, sInfo)=(t[0],t[1],t[3],t[4],t[7])
			
			lInfo=sInfo.split(";")
			for sCategory in lInfo:
				if "=" in sCategory:
					(sHeader,sData)=sCategory.split("=")
					if sHeader=="AC":
						sAC=sData
					elif sHeader=="AF":
						sAF=sData
					elif sHeader=="AN":
						sAN=sData
					else:
						pass
				else:
					pass
			lAlt=sAlt.split(",")
			lAC=sAC.split(",")
			lAF=sAF.split(",")
			#lAN=sAN.split(",")
			
			for sAltNum in range(0,len(lAlt)):
				try:
					fout.write("{0}\n".format("\t".join([sChr,nPosition,".",sRef,lAlt[sAltNum],".",".",\
					str(lAC[sAltNum])+"|"+str(lAF[sAltNum])+"|"+str(sAN)])))
				except:
					print(sLine)
					sys.exit()
	



#output form
#Chr     Start   End     Ref     Alt     AF|AC|NofSample
#10      123810032       123810032       C       T       1.0|2|1
#10      133967449       133967449       C       T       1.0|2|1
#11      124489539       124489539       G       A       1.0|2|1





if __name__=="__main__":
	#sInputFile="liftover_chr21.vcf"
	
	os.chdir("/mnt/QNAP/reference/Korea1K/hg19")
	lFilelist=glob.glob("*.vcf")
	for sInputFile in lFilelist:
		sForANNOVAR(sInputFile)
	
	
	
	
	