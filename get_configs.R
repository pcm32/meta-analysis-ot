#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(ExpressionAtlasInternal))

option_list = list(
  make_option(
    c("-d", "--datasets_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Experiments."
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

fread(input=opt$datasets_file)->datasets

rna_seq_diff<-c()
microarray<-c()

for( accession in datasets$accession ) {
  xmlConfigFilename<-paste0(opt$outdir,"/",accession,"-configuration.xml")
  download.file(url=paste0("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/",accession,"/",accession,"-configuration.xml"), destfile = xmlConfigFilename)
  experimentConfig <- parseAtlasConfig( xmlConfigFilename )
  experimentType <- experimentConfig$experimentType
  cat(experimentType)
  if(grepl( pattern = "rnaseq", experimentType )) {
    rna_seq_diff<-append(rna_seq_diff, accession)
  } else if(grepl( pattern = "array", experimentType )) {
    microarray<-append(microarray, accession)
  }
}

saveRDS(rna_seq_diff, file=paste0(opt$outdir,"/rna_seq_diff.rds") )
write.table(rna_seq_diff, file = paste0(opt$outdir, "/rna_seq_diff.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
saveRDS(microarray, file=paste0(opt$outdir,"/microarray.rds"))
write.table(microarray, file = paste0(opt$outdir, "/microarray.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

