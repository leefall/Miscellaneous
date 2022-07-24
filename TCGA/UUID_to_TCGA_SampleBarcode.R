#!/usr/bin/env R
#BiocManager::install("GenomicDataCommons")
#BiocManager::install("TCGAutils")
library(TCGAutils)
manifest <- read.table("NoC239_ForReview_mc3_NonNA_gdc_manifest_20211116_094540.txt",header=TRUE,stringsAsFactors=FALSE)
#exampleUUID
QueryResult<-UUIDtoBarcode(manifest$id, from_type = "file_id",legacy=TRUE)
Result<-data.frame(manifest,aliquotID=QueryResult$associated_entities.entity_submitter_id )

#exampleUUID<-c("792c45e6-c89d-4abb-a353-7869ef6d752f","9881b111-7b51-4f22-a06a-ce442edc6f37")

write.table(Result,file="./NoC239_TCGAaliquotID_table.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")




