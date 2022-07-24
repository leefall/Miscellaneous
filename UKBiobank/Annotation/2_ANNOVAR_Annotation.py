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




class SNV:

	def __init__(self):
		self.sHGVSGenomicChange=""
		self.sHGVSCodingChange=""
		self.sHGVSProteinChange=""
		self.sHGNCGeneSymbol=""
		self.nEntrezGeneID=""
		self.sEnsemblGeneID=""
		self.sChromosome=""
		self.nChromosome=""
		self.nPosition=""
		self.sReferenceAllele=""
		self.sAlternativeAllele=""
		self.sCytogeneticLocation=""
		self.sStrandOrientation=""
		self.nCodon=""
		self.nExon=""
		self.sMolecularEffects=""
		self.sVariantType=""
		self.sGenotype=""
		self.sdbSNP=""
		self.sClinvar=""
		self.sCosmic=""
		self.sFunctionalDomain=""


	def Parse_ScoreResult(self,sLine):
		sLine=sLine.strip()
		t=sLine.split("\t")

		(self.nSIFT, self.nCADD, self.nPolyPhen)=(t[5],t[22],t[7])




	def Parse_Zygosity(self,sZygosity):
		self.sZygosity=sZygosity


	def Parse_Scorelist(self,lList):
		self.nPolyPhen=lList[7]
		self.nCADD=lList[23]

	def Parse_Cosmic(self,lList):
		sCosmic=lList[1]
		sCosmic=sCosmic.split(";")[0]
		sCosmic=sCosmic.split("=")[-1]
		if self.sCosmic=="":
			self.sCosmic=sCosmic
		else:
			self.sCosmic=self.sCosmic+"|"+sCosmic

	def Parse_dbSNP(self,lList):
		sRSID=lList[1]
		self.sRsID=sRSID




	def Parse_Frequency(self,sLine):
		(sChr, nStartPosition, nEndPosition, sRef, sAlt, nPopFreqMax, n1000G_ALL,n1000G_AFR,n1000G_AMR,n1000G_EAS,n1000G_EUR,n1000G_SAS,nExAC_ALL,nExAC_AFR,nExAC_AMR,nExAC_EAS,nExAC_FIN,nExAC_NFE,nExAC_OTH,nExAC_SAS,nESP6500siv2_ALL,nESP6500siv2_AA,nESP6500siv2_EA,nCG46)=\
		sLine.split("\t")

		self.n1000GALL=n1000G_ALL
		self.n1000GEAS=n1000G_EAS
		#self.nExACALL=nExAC_ALL
		#self.nExACEAS=nExAC_EAS
		#self.nExACSAS=nExAC_SAS
		self.ESP6500ALL=nESP6500siv2_ALL
		self.ESP6500EA=nESP6500siv2_EA

	def Parse_ExAC(self,sExACs):
		(nExACALL, nExACEAS)=(sExACs.split(",")[0],sExACs.split(",")[3])


		self.nExACALL=nExACALL
		self.nExACEAS=nExACEAS
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




	def parsesCytoband(self,sCytoband):
		self.sCytogeneticLocation=sCytoband




	def Parse_Deleterious(self,lList):
		self.sGene=lList[5]
		self.sMolecularEffect=lList[6]
		self.nSIFT=lList[9]
		self.nCADD=lList[11]
		self.nPolyPhen=lList[10]
		self.n1000GALL=lList[12]
		self.n1000GEAS=lList[13]
		self.nExACALL=lList[14]
		self.nExACEAS=lList[15]
		self.nExACSAS="."
		self.ESP6500ALL="."
		self.ESP6500EA="."
		self.sZygosity=lList[4]
		self.sCosmic=lList[8]
		self.sRsID=lList[7]


######ANNOVAR######

def Execute_Score_ANNOVAR(sFile):
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_Score_"+sFile+" -remove -protocol dbnsfp35c -operation f -nastring NA")

def Execute_ALL_Frequency_ANNOVAR(sFile):
	#print "Frequency annotation"
	#os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_ALLFrequency_"+sFile+" -remove -protocol ALL.sites.2015_08 -operation f -nastring NA")
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg38 "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -out ./ANNOVAR_Temporaly/Result_1KGFrequency_"+sFile+" ")
	
	
def Execute_AFR_Frequency_ANNOVAR(sFile):
	#print "Frequency annotation"
	#os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_AFRFrequency_"+sFile+" -remove -protocol AFR.sites.2015_08 -operation f -nastring NA")
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_afr -buildver hg38 "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -out ./ANNOVAR_Temporaly/Result_1KGAFRFrequency_"+sFile+" ")
	
def Execute_AMR_Frequency_ANNOVAR(sFile):
	#print "Frequency annotation"
	#os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_AMRFrequency_"+sFile+" -remove -protocol AMR.sites.2015_08 -operation f -nastring NA")
	#os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_amr -buildver hg38 "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -out ./ANNOVAR_Temporaly/Result_1KGAMRFrequency_"+sFile+" ")
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_amr -buildver hg38 "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -out ./ANNOVAR_Temporaly/Result_1KGAMRFrequency_"+sFile+" ")
	
def Execute_EAS_Frequency_ANNOVAR(sFile):
	#print "Frequency annotation"
	#os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_EASFrequency_"+sFile+" -remove -protocol EAS.sites.2015_08  -operation f -nastring NA")
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_eas -buildver hg38 "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -out ./ANNOVAR_Temporaly/Result_1KGEASFrequency_"+sFile+" ")
	
def Execute_EUR_Frequency_ANNOVAR(sFile):
	#print "Frequency annotation"
	#os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_EURFrequency_"+sFile+" -remove -protocol EUR.sites.2015_08  -operation f -nastring NA")
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_eur -buildver hg38 "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -out ./ANNOVAR_Temporaly/Result_1KGEURFrequency_"+sFile+" ")
	
	
def Execute_SAS_Frequency_ANNOVAR(sFile):
	#print "Frequency annotation"
	#os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_SASFrequency_"+sFile+" -remove -protocol SAS.sites.2015_08  -operation f -nastring NA")
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -dbtype 1000g2015aug_sas -buildver hg38 "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -out ./ANNOVAR_Temporaly/Result_1KGSASFrequency_"+sFile+" ")
	

def Execute_KOVA_ANNOVAR(sFile):
        #print "Cosmic Coding annotation"
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg38 -out ./ANNOVAR_Temporaly/Result_KOVA_"+sFile+" -dbtype KOVA_AF ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb")


def Execute_Geneanno_ANNOVAR(sFile):
	#print "Gene Annotation"
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -build hg38 -out ./ANNOVAR_Temporaly/Result_GenneAnno_"+sFile)

def Execute_Cosmic_Coding_ANNOVAR(sFile):
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg19 -out ./ANNOVAR_Temporaly/Result_CosmicCoding_"+sFile+" -dbtype cosmic88_coding ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar")

def Execute_Cosmic_Noncoding_ANNOVAR(sFile):
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg19 -out ./ANNOVAR_Temporaly/Result_CosmicNonCoding_"+sFile+" -dbtype cosmic88_noncoding ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar")

def avSIFT(sFile):
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb -buildver hg19 -out ./ANNOVAR_Temporaly/Result_SIFT_"+sFile+" -remove -protocol avsift -operation f -nastring NA")


def dbSNP(sFile):
	#print "dbSNP annotation"
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/annotate_variation.pl -filter -buildver hg38 -out ./ANNOVAR_Temporaly/Result_dbSNP_"+sFile+" -dbtype avsnp150 ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/")

def Execute_nonTCGAExAC_ANNOVAR(sFile):
	#print "Gene Annotation"
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_nonTCGAExAC_"+sFile+" -remove -protocol exac03nontcga -operation f -nastring NA")
	#os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl -buildver hg19 -out ./ANNOVAR_Temporaly/Result_nonTCGAExAC_"+sFile+" -dbtype exac03nontcga ./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/annovar/humandb/ -remove -protocol -operation f -nastring NA")


def gNomad(sFile):
	os.system("/storage/home/leefall2/clara/ANNOVAR/annovar/table_annovar.pl "+"./"+sFile+" /storage/home/leefall2/clara/ANNOVAR/2018_ANNOVAR/annovar/humandb/ -buildver hg38 -out ./ANNOVAR_Temporaly/Result_gNomad_"+sFile+" -remove -protocol gnomad30_genome -operation f -nastring NA")






######ANNOVAR END#######



def ANNOVAR_Annotation(sFile):

	try:
		os.mkdir("ANNOVAR_Temporaly")
	except:
		pass


	sInputfile="For_ANNOVAR_"+sFile+".txt"
#	
#	Execute_Score_ANNOVAR(sInputfile)
#	
#	Execute_ALL_Frequency_ANNOVAR(sInputfile)
#	Execute_AFR_Frequency_ANNOVAR(sInputfile)
#	Execute_AMR_Frequency_ANNOVAR(sInputfile)
#	Execute_EAS_Frequency_ANNOVAR(sInputfile)
#	Execute_EUR_Frequency_ANNOVAR(sInputfile)
#	Execute_SAS_Frequency_ANNOVAR(sInputfile)
##
##	
#	Execute_Geneanno_ANNOVAR(sInputfile)
###
##

##

#	dbSNP(sInputfile)
	Execute_nonTCGAExAC_ANNOVAR(sInputfile)

	gNomad(sInputfile)


def Assemble_Data_set(sFile):
	dVariantDict=dict()
	#print "Deleterious_SNP/Deleterious_"+sFile
	#sID=sFile.split("_")[-1]
	#sID=sID.split(".")[0]
	fp=open("./Temporaly/For_ANNOVAR_Intersected_"+sFile.split(".vcf")[0]+".vcf","r")
	for sLine in fp.readlines():

		sLine=sLine.strip()
		t=sLine.split("\t")
		#print sLine
		(sChr, nPosition,nEndPosition, sRef, sAlt,sGenotype)=\
		(t[0],t[1],t[2],t[3],t[4],t[5])

		sKey="_".join([sChr,nPosition,nEndPosition,sRef,sAlt])
		if sKey in dVariantDict.keys():
			#dDeleteriousDict[sKey].sTCGAID.add(sID)
			#dVariantDicts[sKey]=[sGenotype,"Somatic"]
			print(sKey)
			pass

		else:
			dVariantDict[sKey]=[sGenotype,"Somatic"]
	fp.close()
	sHaplotypeFile="Sureselect6_vcf_"+sFile.split("_")[1]
	fp=open("/mnt/towel/Ophthalmology/VCF/NonIntersect/Haplotypecaller/Temporaly/For_ANNOVAR_Intersected_"+sHaplotypeFile.split(".vcf")[0]+".vcf")
	#print(dVariantDict)
	for sLine in fp.readlines():
		sLine=sLine.strip()
		t=sLine.split("\t")
		#print sLine
		(sChr, nPosition,nEndPosition, sRef, sAlt,sGenotype,sFilter)=\
		(t[0],t[1],t[2],t[3],t[4],t[5],t[14])
		if sFilter=="LowQual":
			pass
		else:
			sKey="_".join([sChr,nPosition,nEndPosition,sRef,sAlt])
			if sKey in dVariantDict.keys():
				sSomaticGenotype=dVariantDict[sKey][0] # het
				if ((sSomaticGenotype=="het") and (sGenotype=="hom")):
					sSomaticStatus="Germline"
					dVariantDict[sKey]=[sGenotype,sSomaticStatus]
				elif ((sSomaticGenotype=="het") and (sGenotype=="het")):
					sSomaticStatus="Germline"
					dVariantDict[sKey]=[sGenotype,sSomaticStatus]
				elif ((sSomaticGenotype=="hom") and (sGenotype=="het")):
					pass
				elif ((sSomaticGenotype=="hom") and (sGenotype=="hom")):
					sSomaticStatus="Germline"
					dVariantDict[sKey]=[sGenotype,sSomaticStatus]

			else:
				dVariantDict[sKey]=[sGenotype,"Germline"]

	fout=open("/mnt/towel/Ophthalmology/VCF/NonIntersect/Unified/"+sFile.split("_")[-1].replace(".vcf.gz","")+".txt","w")

	for sKey in dVariantDict.keys():
		[sGenotype,sSomaticStatus]=dVariantDict[sKey]
		fout.write("{0}\t{1}\t{2}\n".format(sKey.replace("_","\t"),sGenotype,sSomaticStatus))





	#return dDeleteriousDict




def Filter_VCF(sFile,sTarget):
	try:
		os.mkdir("Temporaly")
	except:
		pass

	os.system("intersectBed -a ./"+sFile+" -b "+sTarget+" > ./Temporaly/Intersected_"+sFile.split(".vcf")[0]+".vcf")
#	fp=open("./Temporaly/Intersected_"+sFile,"r")
#	fout=open("./Temporaly/Filtered_"+sFile,"w")
#	#dFilterdic={"clustered_events":1,"clustered_read_position":1,"homologous_mapping_event":1,"multi_event_alt_allele_in_normal":1,"str_contraction":1,"strand_artifact":1,"t_lod_fstar":1}
#	#dGermline={
#	for sLine in fp.xreadlines():
#		if sLine[0]=="#":
#			fout.write(sLine)
#		else:
#			sFilter=sLine.split("\t")[6]
#			sGenotype=sLine.split("\t")[-2]
#			sGenotype=sGenotype.split(":")[0]
#			if sFilter=="PASS":
#				if not ((sGenotype=="./.") or (sGenotype=="0/0") or (sGenotype=="./0") or (sGenotype=="0/.")):
#					fout.write(sLine)
#
#	fout.close()


def VCF_to_Variant(sFile):
	sTarget="/mnt/QNAP/leefall2/S07604514_Regions.bed"
	#Filter_VCF(sFile,sTarget)
	#Assemble_Data_set(sFile)
#	ANNOVAR_Annotation(sFile)


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
		#VCF_to_Variant(value)
		#ParseVariantperChromosome(value)
		ANNOVAR_Annotation(value)




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

	number_of_cpus=10


	#os.chdir("/second_storage/leefall2/CancerPanel/vcf/Ion/Target_"+sDir)
	os.chdir("/mnt/towel/leefall2/UKBiobank/ANNOVAR/Chromosome")



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
