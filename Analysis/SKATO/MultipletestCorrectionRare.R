#!/usr/bin/env R
library("optparse")

option_list = list(
make_option(c("-I", "--inputFile"), type="character",help="Input file Name", metavar="character"),
make_option(c("-O", "--outdir"), type="character",help="Output Directory", metavar="character")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

Input.File<-opt$inputFile
output.dir<-opt$outdir

sFisher<-data.frame(read.table(Input.File,header=FALSE,stringsAsFactors = FALSE))

colnames(sFisher)<-c("Gene","Pvalue")

sFisher$FDR<-p.adjust(sFisher$Pvalue,method="fdr")
sFisher$bon<-p.adjust(sFisher$Pvalue,method="bonferroni")


Sig.Vanila<-sFisher[which(sFisher$Pvalue<0.05),]
Sig.FDR<-sFisher[which(sFisher$FDR<0.05),]
Sig.bon<-sFisher[which(sFisher$bon<0.05),]


Input.File.Header<-strsplit(strsplit(Input.File,"/")[[1]][length(strsplit(Input.File,"/")[[1]])],fixed=TRUE,".")[[1]][1]


write.table(sFisher,file=paste(output.dir,"/",Input.File.Header,"_All.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(Sig.Vanila,file=paste(output.dir,"/",Input.File.Header,"_JustP.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(Sig.FDR,file=paste(output.dir,"/",Input.File.Header,"_FDR.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
write.table(Sig.bon,file=paste(output.dir,"/",Input.File.Header,"_bon.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
