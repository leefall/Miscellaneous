#!/usr/bin/env python
import os, subprocess, glob, time
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





	
def bwa(sFile):
#	fout=open("Code_"+sFile+".txt","w")
#	try:
#		os.mkdir("Preprocess")
#	except:
#		pass
#	try:
#		os.mkdir("Final_Bam")
#	except:
#		pass
	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''
	
	

bwa mem -t 8 -R "@RG\\tID:'''+sID+'''\\tPL:illumina\\tLB:Exome\\tSM:'''+sID+'''" /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta '''+sFile+''' '''+sPartner+''' > /mnt/towel/BRCA/Preprocess/'''+sID+'''.sam

samtools view -bhS /mnt/towel/BRCA/Preprocess/'''+sID+'''.sam > /mnt/towel/BRCA/Preprocess/'''+sID+'''.bam
samtools sort /mnt/towel/BRCA/Preprocess/'''+sID+'''.bam -o /mnt/towel/BRCA/Preprocess/sorted'''+sID+'''.bam


''',shell=True)

def MarkDuplicate(sFile):

	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	
	subprocess.call('''

java -Xmx8G -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary -jar /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/Tools/picard_tools_1.108/MarkDuplicates.jar \
AS=TRUE \
I=/mnt/towel/BRCA/Preprocess/sorted'''+sID+'''.bam \
O=/mnt/towel/BRCA/Preprocess/dedup_sorted'''+sID+'''.bam \
METRICS_FILE=/mnt/towel/BRCA/Preprocess/'''+sID+'''_duplicates \
REMOVE_DUPLICATES=true \
CREATE_INDEX=True


''',shell=True)

def BaseRecalibrate(sFile):

	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''


### Base Quality Score Recalibration
##  1) Base Recalibrator


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" BaseRecalibrator \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-I /mnt/towel/BRCA/Preprocess/dedup_sorted'''+sID+'''.bam \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
--known-sites /mnt/QNAP/reference/GATK/hg19/dbsnp_138.hg19.vcf \
-O /mnt/towel/BRCA/Preprocess/recal_'''+str(sID)+'''.table

## 2) Apply Recalibrate

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" ApplyBQSR \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-I /mnt/towel/BRCA/Preprocess/dedup_sorted'''+sID+'''.bam \
--bqsr-recal-file /mnt/towel/BRCA/Preprocess/recal_'''+str(sID)+'''.table \
-O /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sID)+'''.bam


rm /mnt/towel/BRCA/Preprocess/'''+sID+'''.sam


''',shell=True)



def IntersectBed(sFile):

	
	
	
	sID=sFile.split(".")[0]
	sID=sID.split("_")[1]
	#sTumorFile=sID+"-T"
	#sNormalFile=sID+"-N"
	

	subprocess.call('''
	


bedtools intersect -a ./'''+sFile+''' -b /mnt/towel/Ophthalmology/Code/BLAST/KMT2B.bed > /mnt/towel/Ophthalmology/Preprocess/postrecal/KMT2B/KMT2B_'''+str(sID)+'''.bam


''',shell=True)



def GATKHaplotypecaller(sFile):

	sID=sFile.split("-")[0]
	sTail=sFile.split("_")[1]
	sTN=sFile.split("_")[0]
	sTN=sTN.split("-")[-1]
	sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")
	
	sPartner=sFile.split("_")[0]+"_"+sTail
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" HaplotypeCaller \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I /mnt/towel/BRCA/Preprocess/Intersected_postrecal_'''+str(sID)+'''.bam \
-O /mnt/towel/BRCA/VCF/Sureselect6_vcf_'''+str(sID)+'''.vcf.gz 



#Remove temporal data


''',shell=True)
#java -Xmx2G -Djava.io.tmpdir=/storage/home/leefall2/mypro/Cancer_Panel_Package/Package/temporary -jar /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/Tools/GenomeAnalysisTK_2.8_1_g932cd3a/GenomeAnalysisTK.jar \
#-R /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/ucsc.hg19.fasta \
#-T UnifiedGenotyper \
#-I ./Preprocess/recal_realign_dedup_sorted_Reordered_'''+str(sID)+'''.bam \
#-o ./GATK_'''+str(sID)+'''.vcf \
#--dbsnp /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19/dbsnp_137.hg19.vcf \
#-stand_call_conf 30 \
#-stand_emit_conf 10 \
#-glm BOTH

def CollectOxogArtifact(sFile):

	sID=sFile.split("_")[0]
	sTail=sFile.split("_")[1]
	#sTN=sFile.split("_")[0]
	#sTN=sTN.split("-")[-1]
	#sID=sID+"-"+sTN
	sTail=sTail.replace("1","2")

	sPartner=sFile.split("_")[0]+"_"+sTail


	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.0.0/gatk --java-options "-Xmx8G" CollectSequencingArtifactMetrics \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sID)+'''.bam \
-O /mnt/towel/BRCA/Preprocess/Preprocess/oxog_'''+str(sID)+''' \
--FILE_EXTENSION .txt



''',shell=True)


def Mutect2(sFile):

	sID=sFile.split("-")[0]
	sID=sID.split("_")[1]
	sTumorFile=sID+"-T"
	sNormalFile=sID+"-N"
	
	#print('''/mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam''')
	#print('''/mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam''')
	subprocess.call('''

/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" Mutect2 \
--tmp-dir=/mnt/QNAP/leefall2/gatk-4.1.0.0/tmp \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam \
-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam \
--normal-sample '''+str(sNormalFile)+''' \
--germline-resource /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
--read-filter NotSecondaryAlignmentReadFilter \
--read-filter NotSupplementaryAlignmentReadFilter \
-L /mnt/towel/BRCA/SureSelect_kinome.bed \
-pon /mnt/QNAP/reference/GATK/hg19/lifted_over_gatk4_mutect2_4136_pon.vcf.gz \
--f1r2-tar-gz /mnt/towel/BRCA/Preprocess/f1r2_'''+str(sID)+'''.tar.gz \
-O /mnt/towel/BRCA/VCF/GATK4rawVCF/rawMutect2_'''+str(sID)+'''.vcf.gz


#Remove temporal data


''',shell=True)


def PileupSummary(sFile):

	sID=sFile.split("-")[0]
	sID=sID.split("_")[1]
	sTumorFile=sID+"-T"
	sNormalFile=sID+"-N"
	

	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" GetPileupSummaries \
-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam \
-V /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
-L /mnt/towel/BRCA/SureSelect_kinome.bed \
-O /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam.targeted_sequencing.table 

/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" GetPileupSummaries \
-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam \
-V /mnt/QNAP/reference/GATK/hg19/af-only-gnomad.raw.sites.hg19.vcf.gz \
-L /mnt/towel/BRCA/SureSelect_kinome.bed \
-O /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam.targeted_sequencing.table 






''',shell=True)

def CalculateContamination(sFile):

	sID=sFile.split("-")[0]
	sID=sID.split("_")[1]
	sTumorFile=sID+"-T"
	sNormalFile=sID+"-N"
	

	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" CalculateContamination \
-I /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam.targeted_sequencing.table \
-matched /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sNormalFile)+'''.bam.targeted_sequencing.table \
-tumor-segmentation /mnt/towel/BRCA/Preprocess/'''+str(sTumorFile)+'''.segments.table \
-O /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam.contamination.table



''',shell=True)

def LearnReadOrientationModel(sFile):

	sID=sFile.split("-")[0]
	sID=sID.split("_")[1]
	sTumorFile=sID+"-T"
	sNormalFile=sID+"-N"
	

	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" LearnReadOrientationModel \
-I /mnt/towel/BRCA/Preprocess/f1r2_'''+str(sID)+'''.tar.gz \
-O /mnt/towel/BRCA/Preprocess/'''+str(sID)+'''_read-orientation-model.tar.gz


''',shell=True)

def FilterMutactCall(sFile):

	sID=sFile.split("-")[0]
	sID=sID.split("_")[1]
	sTumorFile=sID+"-T"
	sNormalFile=sID+"-N"
	

	subprocess.call('''


/mnt/QNAP/leefall2/gatk-4.1.4.1/gatk --java-options "-Xmx8G" FilterMutectCalls \
-R /mnt/QNAP/reference/GATK/hg19/ucsc.hg19.fasta \
-V /mnt/towel/BRCA/VCF/GATK4rawVCF/rawMutect2_'''+str(sID)+'''.vcf.gz \
--contamination-table /mnt/towel/BRCA/Preprocess/postrecal_'''+str(sTumorFile)+'''.bam.contamination.table \
--tumor-segmentation /mnt/towel/BRCA/Preprocess/'''+str(sTumorFile)+'''.segments.table \
--orientation-bias-artifact-priors /mnt/towel/BRCA/Preprocess/'''+str(sID)+'''_read-orientation-model.tar.gz \
-O /mnt/towel/BRCA/VCF/GATK4FilterVCF/OrientationBiasFilter/OriFiltered_'''+str(sID)+'''.vcf.gz


''',shell=True)



def BamToFASTQ(sFile):

	sID=sFile.split(".")[0]
	sID=sID.split("_")[1]
	#sTumorFile=sID+"-T"
	#sNormalFile=sID+"-N"
	

	subprocess.call('''

samtools bam2fq /mnt/towel/Ophthalmology/Preprocess/postrecal/KMT2B/KMT2B_'''+str(sID)+'''.bam > FASTQ/KMT2B_'''+str(sID)+'''.fastq






''',shell=True)



def FASTQtoFASTA(sFile):

	sID=sFile.split(".")[0]
	sID=sID.split("_")[1]
	#sTumorFile=sID+"-T"
	#sNormalFile=sID+"-N"
	

	subprocess.call('''



sed -n '1~4s/^@/>/p;2~4p' /mnt/towel/Ophthalmology/Preprocess/postrecal/FASTQ/'''+str(sID)+'''.fastq > /mnt/towel/Ophthalmology/Preprocess/postrecal/FASTA/'''+str(sID)+'''.fasta









''',shell=True)


def ParseFASTA(sFile):

	sID=sFile.split(".")[0]
	sID=sID.split("_")[1]
	#sTumorFile=sID+"-T"
	#sNormalFile=sID+"-N"
	
#/mnt/towel/Ophthalmology/Preprocess/postrecal/FASTQ/BLAST/Exact_"+str(sID)+".txt
	
	fp=open("/mnt/towel/Ophthalmology/Preprocess/postrecal/FASTQ/BLAST/Exact_"+str(sID)+".txt")
	fFASTA=open("/mnt/towel/Ophthalmology/Preprocess/postrecal/FASTA/"+str(sID)+".fasta","r")
	
	dReadDict=dict()
	for sLine in fp.readlines():
		t=sLine.split("\t")
		dReadDict[t[0]]=1
	
	
	fout=open("/mnt/towel/Ophthalmology/Preprocess/postrecal/FASTA/ForTublar_"+str(sID)+".fasta","w")
	nCount=0
	sKey=0
	for sLine in fFASTA.readlines():
		nCount+=1
		if nCount%2==0:
			if sKey==1:
				fout.write(sLine)
			else:
				pass
		if nCount%2==1:
			sStrip=sLine.strip()
			sStrip=sStrip[1:]
			if sStrip in dReadDict:
				fout.write(sLine)
				sKey=1
			else:
				sKey=0
				pass
	
	
	
	
	
	
#	subprocess.call('''
#
#
#
#sed -n '1~4s/^@/>/p;2~4p' /mnt/towel/Ophthalmology/Preprocess/postrecal/FASTQ/'''+str(sID)+'''.fastq > /mnt/towel/Ophthalmology/Preprocess/postrecal/FASTA/'''+str(sID)+'''.fasta
#
#
#
#
#
#
#
#
#
#''',shell=True)







def Blastn(sFile):

	sID=sFile.split(".")[0]
	sID=sID.split("_")[1]
	#sTumorFile=sID+"-T"
	#sNormalFile=sID+"-N"
	

	subprocess.call('''

blastn -query /mnt/towel/Ophthalmology/Preprocess/postrecal/FASTA/'''+str(sID)+'''.fasta -db /mnt/towel/Ophthalmology/Code/BLAST/KMT2B_target_100.fa -out /mnt/towel/Ophthalmology/Preprocess/postrecal/FASTQ/BLAST/Exact_'''+str(sID)+'''.txt -evalue 0.00001 -outfmt 6 -task blastn -num_threads 20


''',shell=True)

def producer_task(q, cosmic_dict):
	
	
	#s089=glob.glob("PMSNH089*-T_1.fastq.gz") #10
	#s090=glob.glob("PMSNH090*-T_1.fastq.gz") #10
	#s089=glob.glob("PMSNH089*-N_1.fastq.gz") #4
	#s090=glob.glob("PMSNH090*-N_1.fastq.gz") #6
	
	#s091=glob.glob("PMSNH091*-T_1.fastq.gz")
	#s092=glob.glob("PMSNH092*_1.fastq.gz")#13
	#s085=glob.glob("postrecal_PMSNH077*-T.bam")#11
	#s086=glob.glob("postrecal_PMSNH078*-T.bam")
	#s087=glob.glob("postrecal_PMSNH079*-T.bam")
	#s088=glob.glob("postrecal_PMSNH080*-T.bam")
	#s089=glob.glob("postrecal_PMSNH081*-T.bam")
	#s090=glob.glob("postrecal_PMSNH084*-T.bam")
	#s078=glob.glob("postrecal_PMSNH078*-T.bam")#11
#	s085=glob.glob("PMSNH085*_1.fastq.gz")#11
	#s099=glob.glob("PMSNH099*-T_1.fastq.gz")
	#s096N=glob.glob("PMSNH096*-N_1.fastq.gz")
	#s097N=glob.glob("PMSNH097*-N_1.fastq.gz")
	
	
#	sPiano=["postrecal_08-54084.bam"]
	sPiano=glob.glob("postrecal_*bam")
	
	#sPiano.remove("postrecal_08-54084.bam")
	sFilelist=sPiano
	
	sFilelist=sPiano
	
	#sFilelist=["IonXpress_001_R_2017_07_05_12_17_40_user_S5XL-0059-39-20170704_KPCDx_Test_Auto_user_S5XL-0059-39-20170704_KPCDx_Test_103.bam"]
	for i in sFilelist:
		#value=random.randint(1,len(sFilelist))
		#sCosmicFile=glob.glob("Cosmic_"+i)
		#if not sCosmicFile:
		value=i
		cosmic_dict[value]=None
	
		logger.info("Producer [%s] putting value [%s] into queue.." % (current_process().name, value))
		q.put(value)
#		else:
			#pass


def consumer_task(q, cosmic_dict):
	while not q.empty():
		value=q.get(True, 0.05)
		#a,b=0,1
#		for item in range(value):
#			a,b=b,a+b
		
		
		#bwa(value)
		#MarkDuplicate(value)
		#BaseRecalibrate(value)
		#IntersectBed(value)
		#GATKHaplotypecaller(value)
		#Mutect2(value)
		#CollectOxogArtifact(value)
		#PileupSummary(value)
		#CalculateContamination(value)
		#LearnReadOrientationModel(value)
		#FindTumorSampleName(value)
		#FilterMutactCall(value)
		#Mutect2OnlyTumor(value)
		#FilterVCFcall(value)
		#FilterfromVCFcall(value)
		#FilterBias(value)
		#IntersectBed(value)
		#BamToFASTQ(value)
		#FASTQtoFASTA(value)
		
		#Blastn(value)
		ParseFASTA(value)
		#BlastnPairwise(value)
		
		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))





if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
	os.chdir("/mnt/towel/Ophthalmology/Preprocess/postrecal")
#	number_of_cpus=cpu_count()-2
	number_of_cpus=8
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
	
	
	
	
	
#	os.system("cp /home/Pathology/Ion_Torrent/Cancer_Panel/code/ANNOVAR_annotation_cancer_Panel.py .")
#	os.system("python ANNOVAR_annotation_cancer_Panel.py")
	
	
	
	
	
	
	
	
	print("Start Time")
	print(StartTime)
	print("End Time")
	print(time.ctime())





