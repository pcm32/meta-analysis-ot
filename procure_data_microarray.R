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

microarrays<-readRDS(paste0(opt$inputdir,"/microarray.rds"))

for( accession in microarrays ) {
  strsplit(getURL(paste0("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/",accession,"/"),ftp.use.epsv=TRUE, dirlistonly = TRUE),split = '\n')->listing
  expNormFileName<-listing[[1]][grep(x=listing[[1]],pattern=".*normalized-expressions.*")]
  
  fread(input=paste0("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/",accession,"/",expNormFileName))->exprWAnnot
  fwrite(exprWAnnot[,3:ncol(exprWAnnot)], file=paste0(opt$outdir,"/",expNormFileName,".undecorated"), sep = "\t")
  exprWAnnot$mean<-rowMeans(exprWAnnot[,4:ncol(exprWAnnot)])
  # sort by gene id and mean, ascending
  setkeyv(exprWAnnot, c("Gene ID","mean"))
  # retrieve the last element for each gene (highest mean), leave only probes
  exprWAnnot[,.SD[.N],by=`Gene ID`][,DesignElementAccession]->highestMeanProbePerGene
  exprWAnnot[,c('Gene ID','Gene Name','DesignElementAccession')]->annot
  annot[,DesignElementAccession:=factor(DesignElementAccession)]
  saveRDS(annot, file=paste0(opt$outdir,"/",accession,".annot.rds"))
  saveRDS(highestMeanProbePerGene, file=paste0(opt$outdir,"/",accession,".probes.rds"))
}





