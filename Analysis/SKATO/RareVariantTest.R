#install.packages("SKAT")
library(SKAT)
library("optparse")

option_list = list(
  make_option(c("-I", "--inputFile"), type="character",help="Input file Name", metavar="character"),
  make_option(c("-O", "--outdir"), type="character",help="Output Directory", metavar="character"),
  make_option(c("-G", "--group"), type="character",help="group File", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

Input.File<-opt$inputFile
output.dir<-opt$outdir
Group.File<-opt$group


Input.Matrix<-data.frame(read.table(Input.File,header=FALSE,stringsAsFactors = FALSE))
Group.variable<-data.frame(read.table(Group.File,header=FALSE,stringsAsFactors = FALSE))

aaaa<-list(Input.Matrix,Group.variable[,1])
names(aaaa)<-c("GenomeMatrix","GroupPhenotype")


bbbb<-t(data.matrix(aaaa$GenomeMatrix))


obj<-SKAT_Null_Model(aaaa$GroupPhenotype ~ 1, out_type="D")

Result.p.value<-SKAT(bbbb,obj,method="SKATO")$p.value

Input.File.Header<-strsplit(strsplit(Input.File,"/")[[1]][length(strsplit(Input.File,"/")[[1]])],fixed=TRUE,".")[[1]][1]

write.table(data.frame(Input.File.Header,Result.p.value),file=paste(output.dir,"/Result_",Input.File.Header,".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")










