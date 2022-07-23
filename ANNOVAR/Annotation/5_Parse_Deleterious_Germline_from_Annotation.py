#!/usr/bin/env python
import gzip
import numpy as np
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

sTargetIDFile="/mnt/alpha/leefall2/TCGA_HNSC/mc3/allmc3Barcode.txt"

class SNV:

	def __init__(self):
		
		
		self.sChromosome=""
		self.nPosition=""
		self.sReferenceAllele=""
		self.sAlternativeAllele=""
		self.ESP6500ALL=""
		self.ESP6500EA=""
		self.sGene=""
		self.sMolecularEffects=""
		self.sMolecularChange="NA"
		self.sVariantType=""
		self.sGenotype=""
		self.sRsID=""
		self.n1000G="NA"
		self.n1000GEAS="NA"
		self.n1000GAFR="NA"
		self.n1000GAMR="NA"
		self.n1000GEUR="NA"
		self.n1000GSAS="NA"
		self.nExACALL="NA"
		self.nExACEAS="NA"
		self.nExACAFR="NA"
		self.nExACAMR="NA"
		self.nExACFIN="NA"
		self.nExACNFE="NA"
		self.nExACOTH="NA"
		self.nExACSAS="NA"
		self.ESP6500ALL="NA"
		self.ESP6500AA="NA"
		self.ESP6500EA="NA"
		self.sCosmicNoncoding=""
		self.sCosmicCoding=""
		self.nKOVA="NA"
		self.gNomadALL="NA"
		self.gNomadEAS="NA"
		self.nAllelecount=""
		self.nPolyPhen="NA"
		self.nCADD="NA"
		self.nSIFT="NA"
		self.sClinvar="NA"
		self.nTCGAAF=0
		self.nTCGAAlleleCount=0
		self.nTCGANofSamples=0
		self.nICGCID="NA"
		self.nICGCOccurrence="NA"
		self.sCIVIC="NA"
		self.sOncoKB="NA"
		self.nKorean1KAC="NA"
		self.nKorean1KAF="NA"
		self.nKorean1KAN="NA"
		self.sSomaticStatus="Somatic"
		self.nKoreanDepAF="NA"
		self.nKoreanDepHet="NA"
		self.nKoreanDepHomo="NA"
		self.sINDELcheck="No"
		self.sTCGAGene="NA"
		self.sGTEX="NA"
		





	def Parse_Zygosity(self,sZygosity,sList):
		#(sChr,nPosition,sRef,sAlt,Zygosity)=(t[0],t[1],t[3],t[4],t[5])
		self.sChromosome=sList[0]
		self.nPosition=sList[1]
		self.sReferenceAllele=sList[3]
		self.sAlternativeAllele=sList[4]
		sZygosity=sZygosity.replace("decomposeblocksub","")
		self.sGenotype=sZygosity
		lZygosity=sZygosity.split(",")
		
		nHomocount=0
		nHeterocount=0
		
		for sSample in lZygosity:
			#print(sSample)
			sSample=sSample.replace("decomposeblocksub","")
			try:
				sAllele=sSample.split("_")[1]
			except IndexError:
				print(sSample)
				sys.exit()
			if sAllele=="het":
				nHeterocount+=1
			elif sAllele=="hom":
				nHomocount+=1
			else:
				print("ERROR")
				print(sZygosity)
				sys.exit()
		
		#Homo|Hetero
		self.nAllelecount=str(nHomocount)+"|"+str(nHeterocount)
		
		


	def Parse_Scorelist(self,lList):
		self.nSIFT=lList[5]
		self.nPolyPhen=lList[12]
		self.nCADD=lList[42]

	def Parse_CosmicCoding(self,lList):
		sCosmic=lList[1]
		sCosmic=sCosmic.split(";")[0]
		sCosmic=sCosmic.split("=")[-1]
		if self.sCosmicCoding=="":
			self.sCosmicCoding=sCosmic
		else:
			self.sCosmicCoding=self.sCosmicCoding+"|"+sCosmic
			
	def Parse_CosmicNoncoding(self,lList):
		sCosmic=lList[1]
		sCosmic=sCosmic.split(";")[0]
		sCosmic=sCosmic.split("=")[-1]
		if self.sCosmicNoncoding=="":
			self.sCosmicNoncoding=sCosmic
		else:
			self.sCosmicNoncoding=self.sCosmicNoncoding+"|"+sCosmic

	def Parse_dbSNP(self,lList):
		sRSID=lList[1]
		self.sRsID=sRSID


	def Parse_SIFT(self,lList):
		nSIFT=lList[-1]
		self.nSIFT=nSIFT
		

	def Parse_Frequency(self,sLine):
		(sChr, nStartPosition, nEndPosition, sRef, sAlt, nPopFreqMax, n1000G_ALL,n1000G_AFR,n1000G_AMR,n1000G_EAS,n1000G_EUR,n1000G_SAS,nExAC_ALL,nExAC_AFR,nExAC_AMR,nExAC_EAS,nExAC_FIN,nExAC_NFE,nExAC_OTH,nExAC_SAS,nESP6500siv2_ALL,nESP6500siv2_AA,nESP6500siv2_EA,nCG46)=\
		sLine.split("\t")

		self.n1000G=n1000G_ALL
		self.n1000GEAS=n1000G_EAS
		self.n1000GAFR=n1000G_AFR
		self.n1000GAMR=n1000G_AMR
		self.n1000GEUR=n1000G_EUR
		self.n1000GSAS=n1000G_SAS
		self.nExACALL=nExAC_ALL
		self.nExACEAS=nExAC_EAS
		self.nExACAFR=nExAC_AFR
		self.nExACAMR=nExAC_AMR
		self.nExACFIN=nExAC_FIN
		self.nExACNFE=nExAC_NFE
		self.nExACOTH=nExAC_OTH
		self.nExACSAS=nExAC_SAS
		self.ESP6500ALL=nESP6500siv2_ALL
		self.ESP6500EA=nESP6500siv2_EA
		self.ESP6500AA=nESP6500siv2_EA

	def Parse_ExAC(self,sList):
		(nExACALL, nExACEAS)=(sList[5],sList[8])

	#nNonExACALL
		self.nNonExACALL=nExACALL
		self.nNonExACEAS=nExACEAS
		#self.nExACSAS=nExAC_SAS


	def Parse_gNomad(self,sList):
		(ngNomadALL, ngNomadEAS)=(sList[5],sList[13])


		self.gNomadALL=ngNomadALL
		self.gNomadEAS=ngNomadEAS
		#self.nExACSAS=nExAC_SAS



	def ExonicParse(self,t):
		#t=[lNMIDs,nExons,sCodons,Allele,sGene]

		lNMIDs=t[0]
		nExons=t[1]
		sCodons=t[2]
		Allele=t[3]
		sGene=t[4]
		lCodonNumber=[]


		if "|" in self.sHGVSProteinChange:
			lRefseqPID=self.sHGVSProteinChange.split("|")
		else:
			lRefseqPID=[self.sHGVSProteinChange]


		for n in xrange(len(sCodons)):
			try:
				sCodons[n]=lRefseqPID[n]+":"+sCodons[n]
			except IndexError:
				sCodons[n]=lRefseqPID[0]+":"+sCodons[n]

		if not lNMIDs==[]:
			for sNMID in lNMIDs:
				sCodon=sNMID.split(":")[-1]
				sCodon=sCodon.replace("c.","")
				if "ins" in sCodon:
					sCodon=sCodon.split("ins")[0]
				elif "del" in sCodon:
					sCodon=sCodon.split("del")[0]
				elif "dup" in sCodon:
					sCodon=sCodon.split("dup")[0]
				else:
					sCodon=sCodon[1:-1]


				lCodonNumber.append(sCodon)


		self.nCodon="|".join(lCodonNumber)
		self.nExon="|".join(nExons)
		self.sHGVSCodingChange="|".join(lNMIDs)
		self.sHGVSProteinChange="|".join(sCodons)

	def Parse_Gene(self,sList,dExonicDict,sKey):
		self.sGene=sList[1]
		
		
		
		if "(" in self.sGene:
			self.sGene=self.sGene.split("(")[0]
		
		
		
		
		
		dExonicVariantDict={'nonframeshift insertion':'In_Frame_Ins', 'frameshift deletion':'Frame_Shift_Del', 'nonsynonymous SNV':'Missense_Mutation', 'stoploss':'Nonsense_Mutation', 'stopgain':'Nonsense_Mutation', 'nonframeshift deletion':'In_Frame_Del', 'frameshift insertion':'Frame_Shift_Ins','frameshift substitution':'Frame_Shift_Substitution','nonframeshift substitution':'In_Frame_Shift_Substitution','frameshift block substitution':'Frame_Shift_Substitution','nonframeshift block substitution':'In_Frame_Shift_Substitution','exonic;splicing':'Splicing','splicing':'Splicing'}
		
		
		self.sMolecularEffects=sList[0]
		
		
		
		if "-" in sList[5]:
			self.sVariantType="Insertion"
		elif "-" in sList[6]:
			self.sVariantType="Deletion"
		else:
			if ((len(sList[5])==1)&(len(sList[6])==1)):
				self.sVariantType="SNV"
			elif ((len(sList[5])==2)&(len(sList[6])==2)):
				self.sVariantType="DNP"
			elif ((len(sList[5])==3)&(len(sList[6])==3)):
				self.sVariantType="TNP"
			else:
				print("##############ERROR Parse Gene############")
				print(sList)
				sys.exit()
		
		
		
		#if self.sMolecularEffects=="exonic":
		if ((self.sMolecularEffects=="exonic") or (self.sMolecularEffects=="exonic;splicing")):
			#self.sMolecularEffects=self.sMolecularEffects+"("+dExonicDict[sKey]+")"
			try:
				sMolecularEffect=dExonicDict[sKey]
				self.sMolecularChange=dExonicDict[sKey].split("|")[1]
			except KeyError:
				print("Variant unmatched ERROR!!!!!!!!!!!!!!!!!!!!!!!!")
				print(dExonicDict)
				sys.exit()
			
			#sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
			
			
			if "|" in sMolecularEffect:
				sMolecularEffect=sMolecularEffect.split("|")[0]
			
			(sChr,nPosition,sRef,sAlt)=sKey.split("_")
			
			
			
			
				
			
			
			if sMolecularEffect=="synonymous SNV":
				self.sMolecularEffects="synonymous_SNV"
			elif sMolecularEffect=="unknown":
				#if nPosition=="40904642":
					#print(self.sVariantType)
				if self.sVariantType=="SNV":
					self.sMolecularEffects="unknown"
					(nSIFT,nPolyphen,nCADD)=(self.nSIFT,self.nPolyPhen,self.nCADD)
					#print(nSIFT,nPolyphen,nCADD)
					try:
						if float(nSIFT)<0.3:
							self.sMolecularEffects="Missense_Mutation"
							
					
					except ValueError:
						pass
					
					try:
						if float(nPolyphen)>0.85:
							self.sMolecularEffects="Missense_Mutation"
							
					
					except ValueError:
						pass
					
					try:
						if float(nCADD)>15:
							self.sMolecularEffects="Missense_Mutation"
							
					except ValueError:
						pass
					
					
					
					
				elif self.sVariantType=="Deletion":
					if len(sRef)%3==0:
						self.sMolecularEffects="In_Frame_Del"
						
					else:
						self.sMolecularEffects="Frame_Shift_Del"
						
					
					
				elif self.sVariantType=="Insertion":
					
					if len(sAlt)%3==0:
						self.sMolecularEffects="In_Frame_Ins"
						
					else:
						self.sMolecularEffects="Frame_Shift_Del"
						
			else:
				if sMolecularEffect in dExonicVariantDict.keys():
					self.sMolecularEffects=dExonicVariantDict[sMolecularEffect]
					
				else:
					print("Exception Error!!!!!!!!!!!!!!!!!!!!!!")
					print(self.sMolecularEffects)
					print(sMolecularEffect)
					sys.exit()
			
			
			
			
			
			
		
	def Parse_GTEX(self,Result):
		

		self.sGTEX=Result

		
		

	def Parse_TCGA(self,sTCGA):
		
		(nTCGAAF,nTCGAAlleleCount,nTCGANofSamples)=sTCGA.split("|")
		
		self.nTCGAAF=float(nTCGAAF)
		self.nTCGAAlleleCount=int(nTCGAAlleleCount)
		self.nTCGANofSamples=int(nTCGANofSamples)

	def parsesCytoband(self,sCytoband):
		self.sCytogeneticLocation=sCytoband


	def Parse_KOVA(self,sKOVA):
		self.nKOVA=sKOVA
		
	def Parse_Korea1K(self,lList):
		
		lKorean1K=lList[1].split("|")
		
		
		self.nKorean1KAC=float(lKorean1K[0])
		self.nKorean1KAF=float(lKorean1K[1])
		self.nKorean1KAN=float(lKorean1K[2])
		
	def Parse_Dep(self,lList):
		
		lDep=lList[1].split("|")
		
		
		self.nKoreanDepAF=float(lDep[0])
		self.nKoreanDepHet=float(lDep[1])
		self.nKoreanDepHomo=float(lDep[2])
		
		


######ANNOVAR Parse######

def Zygosity_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	#fp=open("./ANNOVAR_Input/For_ANNOVAR_"+sFile.replace("vcf","maf"))
	fp=open(sFile)
	#fp.readline()
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt,Zygosity)=(t[0],t[1],t[3],t[4],t[5])
		if "0" in sRef:
			pass
		elif "0" in sAlt:
			pass
		else:
			sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
			dVariantDict[sKey]=SNV()
			Zygosity=Zygosity.replace("-N","")
			dVariantDict[sKey].Parse_Zygosity(Zygosity,t)
	return dVariantDict






def Score_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_Score_"+sFile+".hg19_multianno.txt")
	fp.readline()
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[0],t[1],t[3],t[4])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_Scorelist(t)
	#print(dVariantDict["21_26979768_C_A"].nPolyPhen)
	#print(dVariantDict["21_26979768_C_A"].nCADD)
	
	return dVariantDict
	
	
	
def Execute_Frequency_ANNOVAR(sFile):
	#print "Frequency annotation"
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/annovar/humandb/ -buildver hg19 -out ./ANNOVAR_Temporaly/Result_Frequency_"+sFile+" -remove -protocol popfreq_all_20150413  -operation f -nastring NA")



def Frequency_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_Frequency_"+sFile+".hg19_multianno.txt")
	t=fp.readline().split("\t")
	#print t
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		#print t
		(sChr, nStartPosition, nEndPosition, sRef, sAlt, nPopFreqMax, n1000G_ALL,n1000G_AFR,n1000G_AMR,n1000G_EAS,n1000G_EUR,n1000G_SAS,nExAC_ALL,nExAC_AFR,nExAC_AMR,nExAC_EAS,nExAC_FIN,nExAC_NFE,nExAC_OTH,nExAC_SAS,nESP6500siv2_ALL,nESP6500siv2_AA,nESP6500siv2_EA,nCG46)=\
		sLine.split("\t")
		
		sKey="_".join([sChr.replace("chr",""),nStartPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_Frequency(sLine)
			
	
			
			
	return dVariantDict




def Execute_KOVA_ANNOVAR(sFile):
        #print "Cosmic Coding annotation"
		os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg19 -out ./ANNOVAR_Temporaly/Result_KOVA_"+sFile+" -dbtype KOVA_AF ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb")

def KOVA_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_KOVA_"+sFile+".hg19_KOVA_AF_dropped")
	#t=fp.readline().split("\t")
	#print t
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sExACs, sChr,nPosition,sRef,sAlt)=(t[1],t[2],t[3],t[5],t[6])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_KOVA(sExACs)
			
		else:
			pass
		
		
	return dVariantDict



def Execute_Geneanno_ANNOVAR(sFile):
	#print "Gene Annotation"
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/annovar/humandb/ -build hg19 -out ./ANNOVAR_Temporaly/Result_GenneAnno_"+sFile)

def Gene_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	
	dExonicDict=dict()
	fExon=open("./ANNOVAR_Temporaly/Result_GenneAnno_"+sFile+".exonic_variant_function")
	for sLine in fExon.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[3],t[4],t[6],t[7])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			(sEffect,sHGVS)=(t[1],t[2])
			dExonicDict[sKey]=sEffect+"|"+sHGVS
	
	
	fp=open("./ANNOVAR_Temporaly/Result_GenneAnno_"+sFile+".variant_function")
	#t=fp.readline().split("\t")
	#print t
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[2],t[3],t[5],t[6])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_Gene(t,dExonicDict,sKey)
			
		else:
			pass
		
		
	return dVariantDict



def Execute_Cosmic_Coding_ANNOVAR(sFile):
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg19 -out ./ANNOVAR_Temporaly/Result_CosmicCoding_"+sFile+" -dbtype cosmic88_coding ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar")


def Execute_Cosmic_Coding_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_CosmicCoding_"+sFile+".hg19_cosmic88_coding_dropped")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[2],t[3],t[5],t[6])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_CosmicCoding(t)
		else:
			pass
	fp.close()
	
	
	
		
	return dVariantDict



def Execute_Cosmic_Noncoding_ANNOVAR(sFile):
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg19 -out ./ANNOVAR_Temporaly/Result_CosmicNonCoding_"+sFile+" -dbtype cosmic88_noncoding ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar")



def Execute_Cosmic_Noncoding_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_CosmicNonCoding_"+sFile+".hg19_cosmic88_noncoding_dropped")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[2],t[3],t[5],t[6])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_CosmicNoncoding(t)
		else:
			pass
	fp.close()
	
	
	
		
	return dVariantDict






def avSIFT(sFile):
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb -buildver hg19 -out ./ANNOVAR_Temporaly/Result_SIFT_"+sFile+" -remove -protocol avsift -operation f -nastring NA")


def SIFT_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_SIFT_"+sFile+".hg19_multianno.txt")
	#t=fp.readline().split("\t")
	#print t
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		#(sExACs, sChr,nPosition,sRef,sAlt)=(t[1],t[2],t[3],t[5],t[6])
		(sChr,nPosition,sRef,sAlt)=(t[0],t[1],t[3],t[4])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_SIFT(t)
			
		else:
			pass
		
		
	return dVariantDict



def dbSNP(sFile):
	#print "dbSNP annotation"
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg19 -out ./ANNOVAR_Temporaly/Result_dbSNP_"+sFile+" -dbtype avsnp150 ./ANNOVAR_Input/For_ANNOVAR_"+sFile+" /storage/home/leefall2/clara/ANNOVAR/annovar/humandb/")


def RS_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_dbSNP_"+sFile+".hg19_avsnp150_dropped")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[2],t[3],t[5],t[6])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_dbSNP(t)
			
		else:
			pass
	fp.close()
	return dVariantDict




def Execute_nonTCGAExAC_ANNOVAR(sFile):
	#print "Gene Annotation"
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg19 -out ./ANNOVAR_Temporaly/Result_nonTCGAExAC_"+sFile+" -dbtype exac03nontcga ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/annovar/humandb/ -otherinfo")


def NonTCGAExAC_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_nonTCGAExAC_"+sFile+".hg19_multianno.txt")
	#t=fp.readline().split("\t")
	#print t
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		#(sExACs, sChr,nPosition,sRef,sAlt)=(t[1],t[2],t[3],t[5],t[6])
		(sChr,nPosition,sRef,sAlt)=(t[0],t[1],t[3],t[4])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_ExAC(t)
			
		else:
			pass
		
		
	return dVariantDict



def gNomad_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_gNomad_"+sFile+".hg19_multianno.txt")
	#t=fp.readline().split("\t")
	#print t
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		#(sExACs, sChr,nPosition,sRef,sAlt)=(t[1],t[2],t[3],t[5],t[6])
		(sChr,nPosition,sRef,sAlt)=(t[0],t[1],t[3],t[4])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_gNomad(t)
			
		else:
			pass
		
		
	return dVariantDict

def CIVICParse():
	dCIVICGene=dict()

	fp=open("/storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/20200604_CIVIC_nightly-VariantSummaries.tsv")
	fp.readline()
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		sGene=t[2]
		dCIVICGene[sGene]=1
	fp.close()
	
	return dCIVICGene

def OncoKBParse():
	dOncoKBGene=dict()
	fp=open("/storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/oncokb_biomarker_drug_associations.tsv")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		sGene=t[1]
		dOncoKBGene[sGene]=1
	return dOncoKBGene

def TCGAGeneParse():
	dTCGAUVMGene=dict()
	fp=open("/mnt/towel/Ophthalmology/TCGA/VCF/Filter/mc3_somatic/Curated_Genelst.txt")

#	fp=open("/mnt/towel/Ophthalmology/TCGA/VCF/Filter/mc3_somatic/Controlled_Broad_PASS_UVM_Gene.txt")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		sGene=t[0]
		dTCGAUVMGene[sGene]=1
	return dTCGAUVMGene


def Execute_Korea1K_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_Korean1K_"+sFile+".hg19_KoreanGenome_dropped")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[2],t[3],t[5],t[6])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_Korea1K(t)
		else:
			pass
	fp.close()
	
	
	
		
	return dVariantDict

def Execute_Depression_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_Depression_"+sFile+".hg19_Depression_dropped")
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt)=(t[2],t[3],t[5],t[6])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_Dep(t)
		else:
			pass
	fp.close()
	
	
	
		
	return dVariantDict

def TCGA_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_TCGA_"+sFile+".hg19_TCGA_Pool_dropped")
	#t=fp.readline().split("\t")
	#print t
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sExACs, sChr,nPosition,sRef,sAlt)=(t[1],t[2],t[3],t[5],t[6])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			dVariantDict[sKey].Parse_TCGA(sExACs)
			
		else:
			pass
		
		
	return dVariantDict


def GTEX_Parse(sFile,dVariantDict):
	#Polyphen2_HDIV: 7 , CADD : 23
	fp=open("./ANNOVAR_Temporaly/Result_GTEX_"+sFile+".hg19_multianno.txt")
	#t=fp.readline().split("\t")
	#print t
	fp.readline()
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt,sResult)=(t[0],t[1],t[3],t[4],t[-1])
		if sResult=="NA":
			pass
		else:
			sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
			#sKey=sChr+"_"+nPos+"_"+sRef+"_"+sAlt
			#dGTEX[sKey]=sResult
			dVariantDict[sKey].Parse_GTEX(sResult)
		
		#(sExACs, sChr,nPosition,sRef,sAlt)=(t[1],t[2],t[3],t[5],t[6])
		#sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		#if sKey in dVariantDict.keys():
#			dVariantDict[sKey].Parse_TCGA(sExACs)
			
		#else:
		#	pass
		
		
	return dVariantDict



def InclusionCriteria(dVariantDict):
	#Variant is retained if it is present in TCGA or ICGC or gene is present in OncoKB and CIVIC
	
	dPassDict=dict()
	
	dCIVICGene=CIVICParse()
	dOncoKBGene=OncoKBParse()
	dTCGAUVMGene=TCGAGeneParse()
	for sKey in dVariantDict.keys():
		
		
		cVariant=dVariantDict[sKey]
		#sSomaticStatus
		if cVariant.sGene in dCIVICGene:
			cVariant.sCIVIC="CIVIC"
			dPassDict[sKey]=cVariant
			dVariantDict[sKey].sSomaticStatus="Somatic"
		else:
			#dVariantDict[sKey].sSomaticStatus="Somatic"
			pass
		
		
		if cVariant.sGene in dOncoKBGene:
			dPassDict[sKey]=cVariant
			cVariant.sOncoKB="OncoKB"
			dVariantDict[sKey].sSomaticStatus="Somatic"
		else:
			pass
		
		if cVariant.sGene in dTCGAUVMGene:
			dPassDict[sKey]=cVariant
			cVariant.sTCGAGene="TCGA"
			dVariantDict[sKey].sSomaticStatus="Somatic"
		else:
			pass
		
		
		
		
		if not cVariant.nTCGANofSamples==0:
			dPassDict[sKey]=cVariant
			dVariantDict[sKey].sSomaticStatus="Somatic"
			cVariant.sINDELcheck="Yes"
		else:
			pass
		
		if not cVariant.nICGCID=="NA":
			dPassDict[sKey]=cVariant
			dVariantDict[sKey].sSomaticStatus="Somatic"
			cVariant.sINDELcheck="Yes"
		else:
			pass
		
		if not cVariant.sCosmicNoncoding=="":
			dPassDict[sKey]=cVariant
			dVariantDict[sKey].sSomaticStatus="Somatic"
			cVariant.sINDELcheck="Yes"
		else:
			pass
		
		
		if not cVariant.sCosmicCoding=="":
			dPassDict[sKey]=cVariant
			dVariantDict[sKey].sSomaticStatus="Somatic"
			cVariant.sINDELcheck="Yes"
		else:
			pass
		
		
		#if sKey=="20_9376217_T_G":
#			print("PLCB4 Inclusion")
			#print(dVariantDict[sKey].sSomaticStatus)
		
		
	return dPassDict









def TestFilter(sTest,cVariant):
	
	if sTest=="1.":
		sTest=1
	
	if sTest=="0.":
		sTest=0
	
	
	try:
		nTest=float(sTest)
		if not nTest==0:
			cVariant.sINDELcheck="Yes"
		if nTest>0.01:
			return False
		else:
			return True
	except:
		return True
		

def IntergenicFilter(dVariantDict):
	#Variant is retained if it is present in TCGA or ICGC or gene is present in OncoKB and CIVIC
	
	dPassDict=dict()
	

	for sKey in dVariantDict.keys():
		
		
		cVariant=dVariantDict[sKey]
		if cVariant.sMolecularEffects=="intergenic":
			pass
		else:
			dPassDict[sKey]=cVariant

		
		
		
	return dPassDict



		
def ExclusionCriteria1(dVariantDict):
	#Variant Excluded if it is PAD > 0.2%
	
	dPassDict=dict()
	
	
	for sKey in dVariantDict.keys():
		
		
		cVariant=dVariantDict[sKey]
		lTestResult=[]
		sPADlist=[cVariant.n1000G,cVariant.n1000GEAS,cVariant.n1000GAFR,cVariant.n1000GAMR,cVariant.n1000GEUR,cVariant.n1000GSAS,cVariant.nExACALL,cVariant.nExACEAS,cVariant.nExACAFR,cVariant.nExACAMR,cVariant.nExACFIN,cVariant.nExACNFE,cVariant.nExACOTH,cVariant.nExACSAS,cVariant.ESP6500ALL,cVariant.ESP6500AA,cVariant.ESP6500EA,cVariant.nKOVA,cVariant.nKorean1KAF,cVariant.gNomadALL,cVariant.nKoreanDepAF]
		
		for sFrequency in sPADlist:
			sTestResult=TestFilter(sFrequency,cVariant)
			lTestResult.append(sTestResult)
		
		
		
		#print(cVariant.ESP6500ALL)
	
#		if cVariant.nPosition=="103988805":
#			print(lTestResult)
#			print(cVariant.sINDELcheck)
#	
#	
	
		
		if False in lTestResult:
			#print("no")
			dVariantDict[sKey].sSomaticStatus="Germline"
			#cVariant.sINDELcheck="Yes"
			#pass
		else:
			#dVariantDict[sKey].sSomaticStatus="Somatic"
			dPassDict[sKey]=cVariant
	
	return dPassDict



def ExclusionCriteria2(dVariantDict):
	#Variant Excluded if it is PAD > 0.2%
	
	dPassDict=dict()
	
	
	for sKey in dVariantDict.keys():
		
		
		cVariant=dVariantDict[sKey]
		lTestResult=[]
		#sPADlist=[cVariant.nKoreanDepAF]
		
		#for sFrequency in sPADlist:
			#sTestResult=TestFilter(sFrequency,cVariant)
		
		#	lTestResult.append(sTestResult)
		
		
		
		#print(cVariant.ESP6500ALL)
	
		
		#if False in lTestResult:
		if cVariant.nKoreanDepAF!="NA":
			#print("no")
			dVariantDict[sKey].sSomaticStatus="Germline"
			#pass
		else:
			
			dPassDict[sKey]=cVariant
	
	return dPassDict





def ExclusionCriteria3(dVariantDict):
	#Variant is Excluded if it is benign in Clinvar
	
	dPassDict=dict()
	

	for sKey in dVariantDict.keys():
		
		
		cVariant=dVariantDict[sKey]
		
		if "benign" in cVariant.sClinvar:
			dVariantDict[sKey].sSomaticStatus="Germline"
		elif "Benign" in cVariant.sClinvar:
			dVariantDict[sKey].sSomaticStatus="Germline"
		else:
			dPassDict[sKey]=cVariant
		if cVariant.sClinvar!="NA":
			cVariant.sINDELcheck="Yes"
	
	
	
	return dPassDict



def ExclusionCriteria4(dVariantDict):
	#Variant is Excluded if it is benign in Clinvar
	
	dPassDict=dict()
	
	dExonicVariantDict={'Frame_Shift_Del':1,'In_Frame_Ins':1,'Missense_Mutation':1,'Nonsense_Mutation':1,'Nonsense_Mutation':1,'In_Frame_Del':1,'Frame_Shift_Ins':1,'Frame_Shift_Substitution':1,'In_Frame_Shift_Substitution':1}
	for sKey in dVariantDict.keys():
		
		
		cVariant=dVariantDict[sKey]
		sSIFTKey=1
		sPolyphenKey=1
		sCADDKey=1
		
		
		nPolyPhen=cVariant.nPolyPhen
		CADD=cVariant.nCADD
		nSIFT=cVariant.nSIFT
		sRealMolecularEffects=cVariant.sMolecularEffects
		
		if sRealMolecularEffects in dExonicVariantDict:
			
			try:
				if float(nSIFT)<=0.05:
						pass
				else:
						sSIFTKey=0
			except:
					pass
			
			try:
				if float(nPolyPhen)>=0.85:
					pass
				else:
					sPolyphenKey=0
			except:
				pass
			
			try:
				if float(CADD)>=15:
					pass
				else:
					sCADDKey=0
			except:
				pass
				
			
#			if cVariant.nPosition=="7661549":
#				print(nSIFT)
#				print(nPolyPhen)
#				print(CADD)
#				print(sPassKey)
			
			
			
			if ((sSIFTKey==1) or (sPolyphenKey==1) or (sCADDKey==1)):
				dPassDict[sKey]=cVariant
			else:
				dVariantDict[sKey].sSomaticStatus="Germline"
		
		else:
			dVariantDict[sKey].sSomaticStatus="Germline"
	
	
	
	
	
	
	return dPassDict


def ExclusionCriteria5(dVariantDict):
	
	dPassDict=dict()
	for sKey in dVariantDict.keys():
		cVariant=dVariantDict[sKey]
		if cVariant.sGTEX=="NA":
			pass
		else:
			dPassDict[sKey]=cVariant
		
	
	return dPassDict
	




def ExclusionExonCriteria(dVariantDict):
	#Variant is Excluded if it is benign in Clinvar
	
	dPassDict=dict()
	
	dExonicVariantDict={"splicing":1,'Frame_Shift_Del':1,'In_Frame_Ins':1,'Missense_Mutation':1,'Nonsense_Mutation':1,'Nonsense_Mutation':1,'In_Frame_Del':1,'Frame_Shift_Ins':1,'Frame_Shift_Substitution':1,'In_Frame_Shift_Substitution':1}
		
		
	for sKey in dVariantDict.keys():
		
		
		cVariant=dVariantDict[sKey]
		sMolecularEffects=cVariant.sMolecularEffects.split("(")[0]
		#if "benign" in cVariant.sClinvar:
#			dVariantDict[sKey].sSomaticStatus="Germline"
		#elif "Benign" in cVariant.sClinvar:
			#dVariantDict[sKey].sSomaticStatus="Germline"
		if sMolecularEffects in dExonicVariantDict:
			dPassDict[sKey]=cVariant
		else:
			pass
		#if cVariant.sClinvar!="NA":
			#cVariant.sINDELcheck="Yes"
	
	
	
	return dPassDict







def UnifiedWriter(sFile,dVariantDict):
	
	
	lFile=sFile.split("_")
	t=lFile[-1]
	
	
	#print(len(dFrozen))
	#print(dFrozen)
	
	
	try:
		os.mkdir("./Deleterious_Assembled")
	except:
		pass
		
		
	
	lFile=sFile.split("_")
	t=lFile[-1]
	fFinal=open(sTargetIDFile)
	#fFinal=open("/mnt/towel/BRCA/FinalFreeze/ID/Non_Asian_TNBC_ID.txt")
	
	dFinal=dict()
	
	for sLine in fFinal.readlines():
		sLine=sLine.strip()
		dFinal[sLine]=0
	

	dExonicDict={'nonframeshift insertion':'In_Frame_Ins', 'frameshift deletion':'Frame_Shift_Del', 'nonsynonymous SNV':'Missense_Mutation', 'stoploss':'Nonsense_Mutation', 'stopgain':'Nonsense_Mutation', 'nonframeshift deletion':'In_Frame_Del', 'frameshift insertion':'Frame_Shift_Ins','frameshift substitution':'Frame_Shift_Substitution'}
	fout=open("./Deleterious_Assembled/Assembled_"+t,"w")

	for sKey in dVariantDict.keys():
		#print sKey
		
		sForMAFTOOL=0
		
		(sChr,nPosition,sRef,sAlt)=sKey.split("_")
		
		
		
		cVariant=dVariantDict[sKey]
		
		
		sVariantTypeEffect=cVariant.sMolecularEffects.split("(")[0]
		
		
		
		
		if 1:
			
			
			(nTemporAlleleCount,nTemporGenotype)=(cVariant.nAllelecount,cVariant.sGenotype)
			
			lTemporGenotype=nTemporGenotype.split(",")
			lTemporFFPEGenotype=[]
			lTemporFrozenGenotype=[]
			
			
			for sTempor in lTemporGenotype:
				sID=sTempor.split("_")[0]
				if sID in dFinal:
					lTemporFFPEGenotype.append(sTempor)
				else:
					pass
			
			
			if not cVariant.sVariantType=="SNV":
				if len(lTemporFFPEGenotype)==1:
					
					if cVariant.sCosmicCoding=="":
						lTemporFFPEGenotype=[]
				else:
					pass
			
			lUnionGenotype=lTemporFFPEGenotype+lTemporFrozenGenotype
			
			
			if len(lUnionGenotype)!=0:
				
	
				#nAllelecount=str(nHomo)+"|"+str(nHetero)
				
				#sGenotypelist=",".join(lUnionGenotype)
				
				#cVariant.sGene,sChr,nPosition,nPosition,cVariant.sVariant_Classification,"SNP",sRef,sAlt,sAlt, sBarcode
				
				
				#for sUnion in lUnionGenotype:
					
				
				nHomo=0
				nHetero=0
				for sUnion in lUnionGenotype:
					sGenotype=sUnion.split("_")[1]
					if sGenotype=="het":
						nHetero+=1
					else:
						nHomo+=1
				
				nAllelecount=str(nHomo)+"|"+str(nHetero)
				
				sGenotypelist=",".join(lUnionGenotype)
			
			
				
				sGenotype=sUnion.split("_")[1]
				sID=sUnion.split("_")[0]
				if sGenotype=="het":
					sAltone=sRef
				else:
					sAltone=sAlt
				
				sAAChange=cVariant.sMolecularChange
				sAAChange=sAAChange.split(",")[0]
				sAAChange=sAAChange.split(":")[-1]		
				
				fout.write("{0}\n".format(\
				#"\t".join([cVariant.sGene,sChr,nPosition,nPosition,sVariantTypeEffect,cVariant.sVariantType,sRef,sAltone,sAlt,sID,sAAChange])))
				"\t".join(map(str,[sChr,nPosition,sRef,sAlt,cVariant.sGene,cVariant.sMolecularEffects,cVariant.sVariantType,cVariant.sSomaticStatus,cVariant.sMolecularChange,cVariant.sRsID,cVariant.sCosmicCoding,cVariant.sCosmicNoncoding,cVariant.n1000G, cVariant.n1000GEAS, cVariant.n1000GAFR, cVariant.n1000GAMR, cVariant.n1000GEUR, cVariant.n1000GSAS, cVariant.nExACALL, cVariant.nExACEAS, cVariant.nExACAFR, cVariant.nExACAMR, cVariant.nExACFIN, cVariant.nExACNFE, cVariant.nExACOTH, cVariant.nExACSAS, \
				cVariant.ESP6500ALL, cVariant.ESP6500AA, cVariant.ESP6500EA,cVariant.gNomadALL,cVariant.gNomadEAS,cVariant.nKOVA,cVariant.nKorean1KAF,cVariant.nKoreanDepAF,cVariant.nSIFT,cVariant.nPolyPhen,cVariant.nCADD,cVariant.nTCGANofSamples,cVariant.nICGCID,cVariant.sOncoKB,cVariant.sCIVIC,nAllelecount,sGenotypelist,cVariant.sGTEX]))))
				
	
	
	fout.close()
		



######ANNOVAR END#######



def CountVariantDic(dVariantDict):
	
	fFinal=open(sTargetIDFile)
	
	dFinal=dict()
	
	for sLine in fFinal.readlines():
		sLine=sLine.strip()
		dFinal[sLine]=0
	
	
	
	
	
	nNumbeofRealVariant=0
	
	
	for sKey in dVariantDict.keys():
		#print sKey
		(sChr,nPosition,sRef,sAlt)=sKey.split("_")
		
		sTemporID=set()
		
		
		cVariant=dVariantDict[sKey]
		
		
		(nTemporAlleleCount,nTemporGenotype)=(cVariant.nAllelecount,cVariant.sGenotype)
		
		#print(nTemporGenotype)
		
		
		
		lTemporGenotype=nTemporGenotype.split(",")
		lTemporTargetGenotype=lTemporGenotype
		
		#if sKey=="21_33245806_A_G":
#			print(nTemporGenotype)
		#lTemporTargetGenotype=SomaticLoosenFilter(lTemporTargetGenotype,cVariant)
		if len(lTemporTargetGenotype)!=0:
			#dConversionID
			#if cVariant.sSomaticStatus=="Somatic":
				
				
			for sTempor in lTemporGenotype:
				sID=sTempor.split("_")[0]
				if sID in dFinal:
					sTemporID.add(sTempor)
				else:
					pass
				
				
				
			nNumbeofRealVariant+=len(sTemporID)
			
			
		
	
	return nNumbeofRealVariant
	
	


###########Unifying Variant informations
def Unifier(sFile):
	
	try:
		os.mkdir("Filtering_set")
	except:
		pass
		
	dVariantDict=dict()
	#sFile Example) Somatic_PMSNH0770.vcf
	#n=1
	
	sInputFile="Gustave_For_ANNOVAR_chr"+sFile+".txt"
	
	
	Zygosity_Parse(sInputFile,dVariantDict)
	#print(dVariantDict)
	Score_Parse(sInputFile,dVariantDict)
	Frequency_Parse(sInputFile,dVariantDict)
	KOVA_Parse(sInputFile,dVariantDict)
	#SIFT_Parse(sInputFile,dVariantDict)
	Gene_Parse(sInputFile,dVariantDict)
	#print(dVariantDict["21_9825834_C_T"].nKOVA)
	#NonTCGAExAC_Parse(sInputFile,dVariantDict)
	#print(dVariantDict["21_11098708_CAGCCGCCA_-"].nExACEAS)
	Execute_Cosmic_Coding_Parse(sInputFile,dVariantDict)
#	
	Execute_Cosmic_Noncoding_Parse(sInputFile,dVariantDict)
#
#	
#	
#	
	RS_Parse(sInputFile,dVariantDict)
	gNomad_Parse(sInputFile,dVariantDict)
	Execute_Korea1K_Parse(sInputFile,dVariantDict)
	Execute_Depression_Parse(sInputFile,dVariantDict)
	TCGA_Parse(sInputFile,dVariantDict)
	GTEX_Parse(sInputFile,dVariantDict)
	#print(len(dVariantDict))
	
	dNonIntergenicDict=IntergenicFilter(dVariantDict)
	
	#Classification 
	
	#InclusionCriteria(dNonIntergenicDict)
	#ExclusionCriteria1(dNonIntergenicDict)
	#ExclusionCriteria2(dNonIntergenicDict)
	
#	dInclusionDict=InclusionCriteria(dNonIntergenicDict)

	fout=open("Filtering_set/Filtering_"+sFile+".txt","w")
	dExclusionDict3=ExclusionExonCriteria(dNonIntergenicDict)
	
	dExclusionDict4=ExclusionCriteria4(dExclusionDict3)
	dExclusionDict5=ExclusionCriteria5(dExclusionDict4)
	fout.write("Non-Filtering\n")
	fout.write("{0}\n".format(CountVariantDic(dVariantDict)))
	
	fout.write("Non-Filtering-loci\n")
	fout.write("{0}\n".format(len(dVariantDict)))
	
	
	fout.write("Non-Intergenic\n")
	fout.write("{0}\n".format(CountVariantDic(dNonIntergenicDict)))
	fout.write("Non-Intergenic-loci\n")
	fout.write("{0}\n".format(len(dNonIntergenicDict)))
	
	
	
	fout.write("Exonic_variants\n")
	fout.write("{0}\n".format(CountVariantDic(dExclusionDict3)))
	fout.write("Exomic_loci\n")
	fout.write("{0}\n".format(len(dExclusionDict3)))
	
	
	#fout.write("Exclusion-1\n")
	#fout.write("{0}\n".format(CountVariantDic(dExclusionDict1)))
	#fout.write("Exclusion-2\n")
	#fout.write("{0}\n".format(CountVariantDic(dExclusionDict2)))
	#fout.write("Exclusion-3\n")
	#fout.write("{0}\n".format(CountVariantDic(dExclusionDict3)))
	fout.write("Deleterious-Variant\n")
	fout.write("{0}\n".format(CountVariantDic(dExclusionDict4)))
	fout.write("Deleterious-loci\n")
	fout.write("{0}\n".format(len(dExclusionDict4)))
	
	fout.write("GTEX-Variant\n")
	fout.write("{0}\n".format(CountVariantDic(dExclusionDict5)))
	fout.write("GTEX-loci\n")
	fout.write("{0}\n".format(len(dExclusionDict5)))
	
	
	
#	#fout.write("Exclusion-4\n")
	#fout.write("{0}\n".format(len(dExclusionDict4)))
	#UnifiedWriterAll(sInputFile,dNonIntergenicDict)
	UnifiedWriter(sInputFile,dExclusionDict5)
	#MAFtoolsWriter(sInputFile,dExclusionDict4)
	#MAFtoolsAllWriter(sInputFile,dNonIntergenicDict)









	#return dDeleteriousDict



def VCF_to_Variant(sFile):
	sTarget="/mnt/QNAP/leefall2/S07604514_Regions.bed"
	#Filter_VCF(sFile,sTarget)
	#Assemble_Data_set(sFile)
	#ANNOVAR_Annotation(sFile)
	Unifier(sFile)


def ParseVariantperChromosome(sChromosome):
	lFilelist=glob.glob("/mnt/towel/Ophthalmology/VCF/NonIntersect/Unified/*.txt")
	dVariantDict=dict()
	sChromosome="chr"+sChromosome
	for sFile in lFilelist:
		fp=open(sFile,"r")
		sID=sFile.split("/")[-1]
		sID=sID.split(".tx")[0]
		for sLine in fp.readlines():
			sLine=sLine.strip()
			(sChr, nStartPosition,nEndPosition,sRef,sAlt,sGenotype,sSomaticStatus)=sLine.split("\t")
			if sChr==sChromosome:
				sIDset=sID+"_"+sGenotype+"_"+sSomaticStatus
				sKey="_".join([sChr, nStartPosition,nEndPosition,sRef,sAlt])
				if sKey in dVariantDict.keys():
					dVariantDict[sKey].add(sIDset)
				else:
					dVariantDict[sKey]=set()
					dVariantDict[sKey].add(sIDset)
			else:
				pass
	fout=open("/mnt/towel/Ophthalmology/VCF/NonIntersect/Unified/Chromosome/For_ANNOVAR_"+sChromosome+".txt","w")
	for sKey in dVariantDict.keys():
		lIDlist=list(dVariantDict[sKey])
		fout.write("{0}\t{1}\n".format("\t".join(sKey.split("_")),",".join(lIDlist)))
	fout.close()


def consumer_task(q, cosmic_dict):
	#dDeleteriousDict=dict()
	while not q.empty():
		value=q.get(True, 0.05)
		VCF_to_Variant(value)
		#ParseVariantperChromosome(value)
		#ANNOVAR_Annotation(value)




		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))


def producer_task(q, cosmic_dict):

	#sFilelist=glob.glob("FinalMutect2TumorOnly_08-49256*.vcf.gz")
	lChrlist=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]
	#lChrlist=["21"]

	for i in lChrlist:
		value=i
		cosmic_dict[value]=None
		logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
		q.put(value)
#		else:
			#pass





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()

	number_of_cpus=30


	#os.chdir("/second_storage/leefall2/CancerPanel/vcf/Ion/Target_"+sDir)
	os.chdir("/home/leefall2/TCGA_HNSC_Work/Germline/For_ANNOVAR/Chromosome")



	manager=Manager()
	fibo_dict=manager.dict()
	producer=Process(target=producer_task, args=(data_queue, fibo_dict))
	producer.start()
	producer.join()
	consumer_list=[]
	for i in range(number_of_cpus):
		consumer=Process(target=consumer_task, args=(data_queue,fibo_dict))
		consumer.start()
		consumer_list.append(consumer)

	[consumer.join() for consumer in consumer_list]

	logger.info(fibo_dict)
	
	
	os.system("cat /home/leefall2/TCGA_HNSC_Work/Germline/For_ANNOVAR/Chromosome/Deleterious_Assembled/Assembled_* > /home/leefall2/TCGA_HNSC_Work/Assembled_Germline_Variant.txt")
	
	
