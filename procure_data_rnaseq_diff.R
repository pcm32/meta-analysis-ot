#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(RCurl))

option_list = list(
  make_option(
    c("-i", "--inputdir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Input directory"
  ),
  make_option(
    c("-o", "--outdir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Output directory"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

rna_seq_diff<-readRDS(paste0(opt$inputdir,"/rna_seq_diff.rds"))

for( accession in rna_seq_diff ) {
  cat("Getting ",accession,"\n")
  expNormFileName<-paste0(accession,"-raw-counts.tsv")
  fread(input=paste0("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/",accession,"/",expNormFileName))->exprWAnnot
  # skip Gene name on second column
  fwrite(exprWAnnot[,-2], file=paste0(opt$outdir,"/",expNormFileName,".undecorated"), sep = "\t")
  exprWAnnot[,c('Gene ID','Gene Name')]->annot
  saveRDS(annot, file=paste0(opt$outdir,"/",accession,".annot.rds"))
}