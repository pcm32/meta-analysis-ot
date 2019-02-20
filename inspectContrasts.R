#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(ExpressionAtlasInternal))

option_list = list(
  make_option(
    c("-a", "--accession"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Accession to inspect contrasts for."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

accession <- opt$accession

#xmlConfigFilename<-paste0(opt$outdir,"/",accession,"-configuration.xml")
xmlConfigFilename<-tempfile(fileext = "xml")
download.file(url=paste0("ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/",accession,"/",accession,"-configuration.xml"), destfile = xmlConfigFilename)
experimentConfig <- parseAtlasConfig( xmlConfigFilename )
experimentType <- experimentConfig$experimentType
print(paste0(accession, " is of type ", experimentType) )

analytics <- experimentConfig$allAnalytics[[1]]
expContrasts <- atlas_contrasts( analytics )
expAssayGroups <- assay_groups( analytics )

cat(paste('Accession',"C.ID", "Size", "C. Name", sep="\t"),"\n")

for( expContrast in expContrasts ) {
  # Get the sample sizes for the contrast at hand.
    refAssayGroupID <- reference_assay_group_id( expContrast )
    testAssayGroupID <- test_assay_group_id( expContrast )
    
    refAssayGroup <- expAssayGroups[[ refAssayGroupID ]]
    testAssayGroup <- expAssayGroups[[ testAssayGroupID ]]
    
    referenceSamplesSize<-length(refAssayGroup@biological_replicates)
    testSamplesSize<-length(testAssayGroup@biological_replicates)
    
    cat(paste(accession, expContrast@contrast_id, paste0(referenceSamplesSize,":",testSamplesSize), expContrast@contrast_name, sep="\t"), "\n")
}
      
