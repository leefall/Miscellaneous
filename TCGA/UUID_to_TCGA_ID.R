#!/usr/bin/env R
#BiocManager::install("GenomicDataCommons")
library(GenomicDataCommons)
manifest <- read.table("./gdc_manifest_20211028_130109.txt",header=TRUE,stringsAsFactors=FALSE)
head(manifest)
files.uuid<-manifest$id
TCGAtranslateID = function(file_ids, legacy = TRUE) {
  info = files(legacy = legacy) %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  # The mess of code below is to extract TCGA barcodes
  # id_list will contain a list (one item for each file_id)
  # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  # so we can later expand to a data.frame of the right size
  barcodes_per_file = sapply(id_list,length)
  # And build the data.frame
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    submitter_id = unlist(id_list)))
}

res = TCGAtranslateID(files.uuid)
#print(res)
write.table(res,file="./UCEC_TCGAID_table.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")




