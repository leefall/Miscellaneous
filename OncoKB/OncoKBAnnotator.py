#!/usr/bin/env python
import numpy as np
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



def VCF2MAF(sFile):
	
	os.system('''
	perl /storage/home/leefall2/tools/vcf2maf/vcf2maf-1.6.19/vcf2maf.pl --ref-fasta /storage/home/leefall2/mypro/Cancer_Panel_Package/Package/data/hg19.fasta --input-vcf /mnt/towel/BRCA/VCF/Parse/ForOncoKB/'''+sFile+''' --output-maf /mnt/towel/BRCA/VCF/Parse/ForOncoKB/MAF/'''+sFile.replace(".vcf",".maf")
	)
	


def MAFtoOncokbInput(sFile):
	fmaf=open("/mnt/towel/BRCA/VCF/Parse/ForOncoKB/MAF/"+sFile.replace(".vcf",".maf"))
	fout=open("/mnt/towel/BRCA/VCF/Parse/ForOncoKB/OncoKBInput/"+sFile.replace(".vcf",".maf"),"w")
	fmaf.readline()
	fmaf.readline()
	sID=sFile.split(".vcf")[0]
	fout.write("NCBI_Build	Hugo_Symbol	Variant_Classification	Tumor_Sample_Barcode	HGVSp_Short	HGVSp	Chromosome	Start_Position	End_Position	Reference_Allele	Tumor_Seq_Allele1	Tumor_Seq_Allele2\n")

	for sLine in fmaf.readlines():
		t=sLine.split("\t")
		(sBuild,sHugo,sClass,sHGVSpShort,sHGVSp,sChr,sStart,sEnd,sRef,sAlt1,sAlt2)=\
		(t[3],t[0],t[8],t[36],t[35],t[4],t[5],t[6],t[10],t[11],t[12])
		fout.write("{0}\n".format("\t".join([sBuild,sHugo,sClass,sID,sHGVSpShort,sHGVSp,sChr,sStart,sEnd,sRef,sAlt1,sAlt2])))
		
		
	fout.close()


def OncoKBannotator(sFile):
	
	
	os.system("python /storage/home/leefall2/tools/oncokb-annotator/MafAnnotator.py -i /mnt/towel/BRCA/VCF/Parse/ForOncoKB/OncoKBInput/"+sFile+" -o /storage/home/leefall2/mypro/Cancer_Panel_Package/vcf/ForOncoKB/OncoKBannotated/oncokb_"+sFile+" -c /mnt/towel/BRCA/Code/OncoKB/Clinical_Information.txt -b a30ced5d-10d2-4358-a586-97a14dab5341")



def consumer_task(q, cosmic_dict):
	dDeleteriousDict=dict()
	while not q.empty():
		value=q.get(True, 0.05)


		VCF2MAF(value)
		MAFtoOncokbInput(value)
		OncoKBannotator(value)
		cosmic_dict[value]="complete"
		logger.info("consumer [%s] getting value [%s] from queue..." % (current_process().name, value))




if __name__=="__main__":
	StartTime=(time.ctime())
	data_queue=Queue()
#	number_of_cpus=cpu_count()-2
	os.chdir("/mnt/towel/BRCA/VCF/Parse/ForOncoKB/OncoKBInput")
	number_of_cpus=10
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

