#!/usr/bin/env python
import gzip
import numpy as np
import os, subprocess, glob, time, tabix
import sys, time, random, re ,requests, logging, glob
import concurrent.futures
from multiprocessing import Process, Queue, Pool, cpu_count, current_process, Manager
sWorDir="/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA"
sMC3File=sWorDir+"/mc3/LoF_MC3.txt"
sReflectFile=sWorDir+"/Reflcetion_Assembled_Germline_Variant.txt"

logger=logging.getLogger()
logger.setLevel(logging.DEBUG)
formatter=logging.Formatter("%(asctime)s - %(message)s")

ch=logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.setFormatter(formatter)
logger.addHandler(ch)



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
		#self.sGenotype=""
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
		self.sGermlineID=set()
		self.sGermlineGenotype=""
		self.sRLOHID=set()
		self.sRLOHGenotype=""
		self.sLOHID=set()
		self.sLOHGenotype=""
		self.sSomaticID=set()
		self.sSomaticGenotype=""
		self.sGTEX=""
		self.sIDCount=''
		
		
	def ParseVarscan2(self, sLine):
		(sChr,nPosition,sRef,sAlt,sGene,sMolecularEffects,sVariantType,sSomaticStatus,sMolecularChange,sRsID,sCosmicCoding,sCosmicNoncoding,\
		n1000G, n1000GEAS, n1000GAFR, n1000GAMR, n1000GEUR, n1000GSAS, nExACALL, nExACEAS, nExACAFR, nExACAMR, nExACFIN, nExACNFE, nExACOTH, \
		nExACSAS, ESP6500ALL, ESP6500AA, ESP6500EA,gNomadALL,gNomadEAS,nKOVA,nKorean1KAF,nKoreanDepAF,nSIFT,nPolyPhen,nCADD,nTCGANofSamples,\
		nICGCID,sOncoKB,sCIVIC,nAllelecount,sLOHlist,sRLOHlist)=sLine.split("\t")
		sRLOHlist=sRLOHlist.strip()
		self.sChromosome=sChr
		self.nPosition=nPosition
		self.sReferenceAllele=sRef
		self.sAlternativeAllele=sAlt
		self.sMolecularEffects=sMolecularEffects
		self.sMolecularChange=sMolecularChange
		self.sVariantType=sVariantType
		self.sRsID=sRsID
		self.n1000G=n1000G
		self.n1000GEAS=n1000GEAS
		self.n1000GAFR=n1000GAFR
		self.n1000GAMR=n1000GAMR
		self.n1000GEUR=n1000GEUR
		self.n1000GSAS=n1000GSAS
		self.nExACALL=nExACALL
		self.nExACEAS=nExACEAS
		self.nExACAFR=nExACAFR
		self.nExACAMR=nExACAMR
		self.nExACFIN=nExACFIN
		self.nExACNFE=nExACNFE
		self.nExACOTH=nExACOTH
		self.nExACSAS=nExACSAS
		self.ESP6500ALL=ESP6500ALL
		self.ESP6500AA=ESP6500AA
		self.ESP6500EA=ESP6500EA
		self.sCosmicNoncoding=sCosmicNoncoding
		self.sCosmicCoding=sCosmicCoding
		self.nKOVA=nKOVA
		self.gNomadALL=gNomadALL
		self.gNomadEAS=gNomadEAS
		#self.nAllelecount=""
		self.nPolyPhen=nPolyPhen
		self.nCADD=nCADD
		self.nSIFT=nSIFT
		#self.sClinvar="NA"
		#self.nTCGAAF=0
		#self.nTCGAAlleleCount=0
		self.nTCGANofSamples=nTCGANofSamples
		self.nICGCID=nICGCID
		#self.nICGCOccurrence="NA"
		self.sCIVIC=sCIVIC
		self.sOncoKB=sOncoKB
		#self.nKorean1KAC="NA"
		self.nKorean1KAF=nKorean1KAF
		#self.nKorean1KAN="NA"
		#self.sSomaticStatus="Somatic"
		self.nKoreanDepAF=nKoreanDepAF
		#self.nKoreanDepHet="NA"
		#self.nKoreanDepHomo="NA"
		
		
		
		self.sRLOHGenotype=sLOHlist
		#self.sLOHID=set()
		self.sLOHGenotype=sRLOHlist
		
		#sLOHlist,sRLOHlist
		lLOHlist=sLOHlist.split(",")
		lRLOHlist=sRLOHlist.split(",")
		for LOH in lLOHlist:
			self.sLOHID.add(LOH)
		
		for RLOH in lRLOHlist:
			self.sRLOHID.add(RLOH)
	
	
	
	
	def AddGermline(self, sLine):
		(sChr,nPosition,sRef,sAlt,sGene,sMolecularEffects,sVariantType,sSomaticStatus,sMolecularChange,sRsID,sCosmicCoding,sCosmicNoncoding,n1000G, n1000GEAS, n1000GAFR, n1000GAMR, n1000GEUR, n1000GSAS, nExACALL, nExACEAS, nExACAFR, nExACAMR, nExACFIN, nExACNFE, nExACOTH, nExACSAS, \
		ESP6500ALL, ESP6500AA, ESP6500EA,gNomadALL,gNomadEAS,nKOVA,nKorean1KAF,nKoreanDepAF,nSIFT,nPolyPhen,nCADD,\
		nTCGANofSamples,nICGCID,sOncoKB,sCIVIC,sNewGermlineGenotype,sNewRLOHGenotype,sNewLOHGenotype,sNewSomaticGenotype,sNewIDCount,sGTEX)=sLine.split("\t")
		
		
		self.sRLOHGenotype=sNewRLOHGenotype
		#self.sLOHID=set()
		self.sLOHGenotype=sNewLOHGenotype
		
		#sLOHlist,sRLOHlist
		lLOHlist=sNewRLOHGenotype.split(",")
		lRLOHlist=sNewLOHGenotype.split(",")
		for LOH in lLOHlist:
			self.sLOHID.add(LOH)
		
		for RLOH in lRLOHlist:
			self.sRLOHID.add(RLOH)
		
		
		
		
		self.sGermlineGenotype=sNewGermlineGenotype
		#sTempor
		lGermlinelist=sNewGermlineGenotype.split(",")
		for sGermIDGeno in lGermlinelist:
			(sGermID,sGermGenotype)=sGermIDGeno.split("_")
			
			
			if sGermID in self.sLOHID:
				pass
			elif sGermID in self.sRLOHID:
				pass
			elif sGermID in self.sSomaticID:
				pass
			else:
				self.sGermlineID.add(sGermIDGeno)
		
		
		#self.sGermlineGenotype=sGenotypelist
		
	
	
	def ParseGermline(self, sLine):
		(sChr,nPosition,sRef,sAlt,sGene,sMolecularEffects,sVariantType,sSomaticStatus,sMolecularChange,sRsID,sCosmicCoding,sCosmicNoncoding,n1000G, n1000GEAS, n1000GAFR, n1000GAMR, n1000GEUR, n1000GSAS, nExACALL, nExACEAS, nExACAFR, nExACAMR, nExACFIN, nExACNFE, nExACOTH, nExACSAS, \
		ESP6500ALL, ESP6500AA, ESP6500EA,gNomadALL,gNomadEAS,nKOVA,nKorean1KAF,nKoreanDepAF,nSIFT,nPolyPhen,nCADD,\
		nTCGANofSamples,nICGCID,sOncoKB,sCIVIC,sNewGermlineGenotype,sNewRLOHGenotype,sNewLOHGenotype,sNewSomaticGenotype,sNewIDCount,sGTEX)=sLine.split("\t")
		self.sIDCount=sNewIDCount
		sGenotypelist=sNewGermlineGenotype
		sRLOHlist=sNewRLOHGenotype
		sLOHlist=sNewLOHGenotype
		
		self.sChromosome=sChr
		self.nPosition=nPosition
		self.sReferenceAllele=sRef
		self.sAlternativeAllele=sAlt
		self.sMolecularEffects=sMolecularEffects
		self.sMolecularChange=sMolecularChange
		self.sVariantType=sVariantType
		self.sRsID=sRsID
		self.n1000G=n1000G
		self.n1000GEAS=n1000GEAS
		self.n1000GAFR=n1000GAFR
		self.n1000GAMR=n1000GAMR
		self.n1000GEUR=n1000GEUR
		self.n1000GSAS=n1000GSAS
		self.nExACALL=nExACALL
		self.nExACEAS=nExACEAS
		self.nExACAFR=nExACAFR
		self.nExACAMR=nExACAMR
		self.nExACFIN=nExACFIN
		self.nExACNFE=nExACNFE
		self.nExACOTH=nExACOTH
		self.nExACSAS=nExACSAS
		self.ESP6500ALL=ESP6500ALL
		self.ESP6500AA=ESP6500AA
		self.ESP6500EA=ESP6500EA
		self.sCosmicNoncoding=sCosmicNoncoding
		self.sCosmicCoding=sCosmicCoding
		self.nKOVA=nKOVA
		self.gNomadALL=gNomadALL
		self.gNomadEAS=gNomadEAS
		#self.nAllelecount=""
		self.nPolyPhen=nPolyPhen
		self.nCADD=nCADD
		self.nSIFT=nSIFT
		#self.sClinvar="NA"
		#self.nTCGAAF=0
		#self.nTCGAAlleleCount=0
		self.nTCGANofSamples=nTCGANofSamples
		self.nICGCID=nICGCID
		#self.nICGCOccurrence="NA"
		self.sCIVIC=sCIVIC
		self.sOncoKB=sOncoKB
		#self.nKorean1KAC="NA"
		self.nKorean1KAF=nKorean1KAF
		#self.nKorean1KAN="NA"
		#self.sSomaticStatus="Somatic"
		self.nKoreanDepAF=nKoreanDepAF
		self.sGTEX=sGTEX
		
		self.sGermlineGenotype=sGenotypelist
		
		lGermlinelist=sGenotypelist.split(",")
		for sGermID in lGermlinelist:
			#if sGermID in self.sLOH:
#				pass
			#elif sGermID in self.sRLOH:
				#pass
			#else:
			self.sGermlineID.add(sGermID)
		
		
		
		self.sRLOHGenotype=sLOHlist
		#self.sLOHID=set()
		self.sLOHGenotype=sRLOHlist
		
		#sLOHlist,sRLOHlist
		lLOHlist=sLOHlist.split(",")
		lRLOHlist=sRLOHlist.split(",")
		for LOH in lLOHlist:
			self.sLOHID.add(LOH)
		
		for RLOH in lRLOHlist:
			self.sRLOHID.add(RLOH)
	
	def AddSomatic(self, sLine):
		(sChr,nPosition,sRef,sAlt,sGene,sMolecularEffects,sVariantType,sSomaticStatus,sMolecularChange,sRsID,sCosmicCoding,sCosmicNoncoding,n1000G, n1000GEAS, n1000GAFR, n1000GAMR, n1000GEUR, n1000GSAS, nExACALL, nExACEAS, nExACAFR, nExACAMR, nExACFIN, nExACNFE, nExACOTH, nExACSAS, \
		ESP6500ALL, ESP6500AA, ESP6500EA,gNomadALL,gNomadEAS,nKOVA,nKorean1KAF,nKoreanDepAF,nSIFT,nPolyPhen,nCADD,nTCGANofSamples,nICGCID,sOncoKB,sCIVIC,sGenotypelist)=\
		sLine.split("\t")
		
		
		self.sSomaticGenotype=sGenotypelist
		lSomaticlist=sGenotypelist.split(",")
		for sSomaticIDGenome in lSomaticlist:
			(sSomaticID,sSomaticGenotype)=sSomaticIDGenome.split("_")
			
			
			if sSomaticID in self.sLOHID:
				pass
			elif sSomaticID in self.sRLOHID:
				pass
			#elif sGermID in self.sRLOH:
#				pass
			else:
				self.sSomaticID.add(sSomaticIDGenome)
	
	def ParseSomatic(self, sLine):
		#(sChr,nPosition,sRef,sAlt,sGene,sMolecularEffects,sVariantType,sSomaticStatus,sMolecularChange,sRsID,sCosmicCoding,sCosmicNoncoding,\
		#n1000G, n1000GEAS, n1000GAFR, n1000GAMR, n1000GEUR, n1000GSAS, nExACALL, nExACEAS, nExACAFR, nExACAMR, nExACFIN, nExACNFE, nExACOTH, \
		#nExACSAS, ESP6500ALL, ESP6500AA, ESP6500EA,gNomadALL,gNomadEAS,nKOVA,nKorean1KAF,nKoreanDepAF,nSIFT,nPolyPhen,nCADD,nTCGANofSamples,\
		#nICGCID,sOncoKB,sCIVIC,nAllelecount,sGenotypelist)=sLine.split("\t")
		(sChr,nPosition,sRef,sAlt,sGene,sMolecularEffects,sVariantType,sSomaticStatus,sMolecularChange,sRsID,sCosmicCoding,sCosmicNoncoding,n1000G, n1000GEAS, n1000GAFR, n1000GAMR, n1000GEUR, n1000GSAS, nExACALL, nExACEAS, nExACAFR, nExACAMR, nExACFIN, nExACNFE, nExACOTH, nExACSAS, \
		ESP6500ALL, ESP6500AA, ESP6500EA,gNomadALL,gNomadEAS,nKOVA,nKorean1KAF,nKoreanDepAF,nSIFT,nPolyPhen,nCADD,nTCGANofSamples,nICGCID,sOncoKB,sCIVIC,sGenotypelist)=\
		sLine.split("\t")
		
		
		self.sChromosome=sChr
		self.nPosition=nPosition
		self.sReferenceAllele=sRef
		self.sAlternativeAllele=sAlt
		self.sMolecularEffects=sMolecularEffects
		self.sMolecularChange=sMolecularChange
		self.sVariantType=sVariantType
		self.sRsID=sRsID
		self.n1000G=n1000G
		self.n1000GEAS=n1000GEAS
		self.n1000GAFR=n1000GAFR
		self.n1000GAMR=n1000GAMR
		self.n1000GEUR=n1000GEUR
		self.n1000GSAS=n1000GSAS
		self.nExACALL=nExACALL
		self.nExACEAS=nExACEAS
		self.nExACAFR=nExACAFR
		self.nExACAMR=nExACAMR
		self.nExACFIN=nExACFIN
		self.nExACNFE=nExACNFE
		self.nExACOTH=nExACOTH
		self.nExACSAS=nExACSAS
		self.ESP6500ALL=ESP6500ALL
		self.ESP6500AA=ESP6500AA
		self.ESP6500EA=ESP6500EA
		self.sCosmicNoncoding=sCosmicNoncoding
		self.sCosmicCoding=sCosmicCoding
		self.nKOVA=nKOVA
		self.gNomadALL=gNomadALL
		self.gNomadEAS=gNomadEAS
		#self.nAllelecount=""
		self.nPolyPhen=nPolyPhen
		self.nCADD=nCADD
		self.nSIFT=nSIFT
		#self.sClinvar="NA"
		#self.nTCGAAF=0
		#self.nTCGAAlleleCount=0
		self.nTCGANofSamples=nTCGANofSamples
		self.nICGCID=nICGCID
		#self.nICGCOccurrence="NA"
		self.sCIVIC=sCIVIC
		self.sOncoKB=sOncoKB
		#self.nKorean1KAC="NA"
		self.nKorean1KAF=nKorean1KAF
		#self.nKorean1KAN="NA"
		#self.sSomaticStatus="Somatic"
		self.nKoreanDepAF=nKoreanDepAF
		
		self.sSomaticGenotype=sGenotypelist
		lSomaticlist=sGenotypelist.split(",")
		#self.sGermlineGenotype=sGenotypelist
		
		#lGermlinelist=sGenotypelist.split(",")
		for sSomaticID in lSomaticlist:
			#if sGermID in self.sLOH:
#				pass
			#elif sGermID in self.sRLOH:
				#pass
			#else:
			self.sSomaticID.add(sSomaticID)
		
	
#	
#sWorDir="/mnt/alpha/leefall2/TCGA_Non_TNBC_BRCA"
#sMC3File=sWorDir+"/mc3/LoF_MC3.txt"
#
#
#dSomaticDict=dict()
#
#fSomatic=open(sMC3File)
#
#for sLine in fSomatic.readlines():
#	sLine=sLine.strip()
#	t=sLine.split("\t")
#	
#	(sChr,nPosition,sRef,sAlt)=(t[0],t[1],t[2],t[3])
#	sKey="_".join([sChr,nPosition,sRef,sAlt])
#	sSomaticID=t[-1]
#	lInformation=t[4:-1]
#	
#	dSomaticDict[sKey]=dict()
#	dSomaticDict[sKey]["ID"]=sSomaticID
#	dSomaticDict[sKey]["Info"]




def ParseGermlineFile(sReflectFile,dLOHVariantDict):
	
	
	sFile=sReflectFile
	fp=open(sFile)
	
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		(sChr,nPosition,sRef,sAlt,sGene,sMolecularEffects,sVariantType,sSomaticStatus,sMolecularChange,sRsID,sCosmicCoding,sCosmicNoncoding,n1000G, n1000GEAS, n1000GAFR, n1000GAMR, n1000GEUR, n1000GSAS, nExACALL, nExACEAS, nExACAFR, nExACAMR, nExACFIN, nExACNFE, nExACOTH, nExACSAS, \
		ESP6500ALL, ESP6500AA, ESP6500EA,gNomadALL,gNomadEAS,nKOVA,nKorean1KAF,nKoreanDepAF,nSIFT,nPolyPhen,nCADD,\
		nTCGANofSamples,nICGCID,sOncoKB,sCIVIC,sNewGermlineGenotype,sNewRLOHGenotype,sNewLOHGenotype,sNewSomaticGenotype,sNewIDCount,sGTEX)=sLine.split("\t")
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		
		if sKey in dLOHVariantDict.keys():
			dLOHVariantDict[sKey].AddGermline(sLine)
		else:
			dLOHVariantDict[sKey]=SNV()
			dLOHVariantDict[sKey].ParseGermline(sLine)
		
		
	
	return dLOHVariantDict

def ParseSomaticFile(sMC3File,dLOHVariantDict):
	
	sFile=sMC3File
	#sFile="/mnt/towel/TCGA/Somatic/ANNOVAR_mc3/Chromosome_LoF/chr"+sChr+".txt"
	fp=open(sFile)
	
	
	for sLine in fp.readlines():
		sLine=sLine.strip()
		#(sChr,nPosition,sRef,sAlt,sGene,sMolecularEffects,sVariantType,sSomaticStatus,sMolecularChange,sRsID,sCosmicCoding,sCosmicNoncoding,\
		#n1000G, n1000GEAS, n1000GAFR, n1000GAMR, n1000GEUR, n1000GSAS, nExACALL, nExACEAS, nExACAFR, nExACAMR, nExACFIN, nExACNFE, nExACOTH, \
		#nExACSAS, ESP6500ALL, ESP6500AA, ESP6500EA,gNomadALL,gNomadEAS,nKOVA,nKorean1KAF,nKoreanDepAF,nSIFT,nPolyPhen,nCADD,nTCGANofSamples,\
		#nICGCID,sOncoKB,sCIVIC,nAllelecount,sGenotypelist)=sLine.split("\t")
		
		(sChr,nPosition,sRef,sAlt)=\
		(sLine.split("\t")[0],sLine.split("\t")[1],sLine.split("\t")[2],sLine.split("\t")[3])
		sKey="_".join([sChr.replace("chr",""),nPosition,sRef,sAlt])
		
		if sKey in dLOHVariantDict.keys():
			dLOHVariantDict[sKey].AddSomatic(sLine)
		else:
			dLOHVariantDict[sKey]=SNV()
			dLOHVariantDict[sKey].ParseSomatic(sLine)
		
		
	
	return dLOHVariantDict




def UnifiedWriter(sWorDir,dVariantDict):
	
	
	#
		
	
	fout=open(sWorDir+"/FinalAssembled_GermlineLOHSomatic.txt","w")

	for sKey in dVariantDict.keys():
		#print sKey
		
		#sForMAFTOOL=0
		
		(sChr,nPosition,sRef,sAlt)=sKey.split("_")
		
		
		
		cVariant=dVariantDict[sKey]
		#self.sGermlineID
		#self.sRLOHID
		#self.sLOHID
		#self.sSomaticID
		
		sIDCount=cVariant.sIDCount
		lIDCount=sIDCount.split("|")
		sGermlineGenotype=",".join(list(cVariant.sGermlineID))
		sRLOHGenotype=",".join(list(cVariant.sRLOHID))
		sLOHGenotype=",".join(list(cVariant.sLOHID))
		sSomaticGenotype=",".join(list(cVariant.sSomaticID))
		#sVariantTypeEffect=cVariant.sMolecularEffects.split("(")[0]
		if sIDCount=='':
			sRealIDCount="|".join(lIDCount[0:3])+"|"+str(len(cVariant.sSomaticID))
		else:
			sRealIDCount="0|0|0"+"|"+str(len(cVariant.sSomaticID))
			
		
		#sIDCount=str(len(cVariant.sGermlineID))+"|"+str(len(cVariant.sRLOHID))+"|"+str(len(cVariant.sLOHID))+"|"+str(len(cVariant.sSomaticID))
		
		
		
		fout.write("{0}\n".format(\
				#"\t".join([cVariant.sGene,sChr,nPosition,nPosition,sVariantTypeEffect,cVariant.sVariantType,sRef,sAltone,sAlt,sID,sAAChange])))
				"\t".join(map(str,[sChr,nPosition,sRef,sAlt,cVariant.sGene,cVariant.sMolecularEffects,cVariant.sVariantType,cVariant.sSomaticStatus,cVariant.sMolecularChange,cVariant.sRsID,cVariant.sCosmicCoding,cVariant.sCosmicNoncoding,cVariant.n1000G, cVariant.n1000GEAS, cVariant.n1000GAFR, cVariant.n1000GAMR, cVariant.n1000GEUR, cVariant.n1000GSAS, cVariant.nExACALL, cVariant.nExACEAS, cVariant.nExACAFR, cVariant.nExACAMR, cVariant.nExACFIN, cVariant.nExACNFE, cVariant.nExACOTH, cVariant.nExACSAS, \
				cVariant.ESP6500ALL, cVariant.ESP6500AA, cVariant.ESP6500EA,cVariant.gNomadALL,cVariant.gNomadEAS,cVariant.nKOVA,cVariant.nKorean1KAF,cVariant.nKoreanDepAF,cVariant.nSIFT,cVariant.nPolyPhen,cVariant.nCADD,cVariant.nTCGANofSamples,cVariant.nICGCID,cVariant.sOncoKB,cVariant.sCIVIC,sGermlineGenotype,sRLOHGenotype,sLOHGenotype,sSomaticGenotype,sRealIDCount,cVariant.sGTEX]))))
				
	





dVariantDict=dict()
	
	
#dLOHVariantDict=ParseLOHFile(sChr,dVariantDict)
	#dGermLOHVariantDict=ParseGermlineFile(sFile,dLOHVariantDict)
dSomaticLOHVariantDict=ParseSomaticFile(sMC3File,dVariantDict)
dGermSomaticLOHVariantDict=ParseGermlineFile(sReflectFile,dSomaticLOHVariantDict)
UnifiedWriter(sWorDir,dGermSomaticLOHVariantDict)