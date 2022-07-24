#!/usr/bin/env python

from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
matrices = matGen.SigProfilerMatrixGeneratorFunc("SNU", "GRCh37", "/mnt/towel/leefall2/Analysis/Signature/SNU", exome=False, bed_file="/mnt/towel/leefall2/Noheader_Sureselect6_Regions.bed", chrom_based=False, plot=True, tsb_stat=False, seqInfo=False)
sig.sigProfilerExtractor("text", "DeleteriousLOHResults", "/mnt/towel/leefall2/Analysis/Signature/SNU/output/SBS/SNU.SBS96.region", reference_genome="GRCh37")



