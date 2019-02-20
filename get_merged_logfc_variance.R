#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(limma))
suppressMessages( library( ExpressionAtlasInternal ) )

option_list = list(
  make_option(
    c("-d", "--datasets_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Experiments."
  ),
  make_option(
    c("-i", "--inputdir"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Input directory with R objects"
  ),
  make_option(
    c("-o", "--outdir"),
    help = "Output directory"
  )
)


opt <- parse_args(OptionParser(option_list = option_list))

fread(input=opt$datasets_file)->datasets

microarrays<-readRDS(paste0(opt$inputdir,"/microarray.rds"))

results<-list()

for( accession in microarrays ) {
  
  experimentConfig <- parseAtlasConfig( paste0( opt$inputdir,"/",accession,"-configuration.xml" ))
  analytics <- experimentConfig$allAnalytics[[1]]
  expContrasts <- atlas_contrasts( analytics )
  expAssayGroups <- assay_groups( analytics )

  readRDS(paste0(opt$inputdir,"/",accession,".annot.rds"))->annot
  readRDS(paste0(opt$inputdir,"/",accession,".probes.rds"))->highestMeanProbePerGene

  for( contrast in unlist(strsplit(x=datasets$contrast[datasets$accession==accession], split = ",")) ) {
    
    for( expContrast in expContrasts ) {
      # Get the sample sizes for the contrast at hand.
      if( expContrast@contrast_id == contrast ) {
        refAssayGroupID <- reference_assay_group_id( expContrast )
        testAssayGroupID <- test_assay_group_id( expContrast )
        
        refAssayGroup <- expAssayGroups[[ refAssayGroupID ]]
        testAssayGroup <- expAssayGroups[[ testAssayGroupID ]]
        
        referenceSamplesSize<-length(refAssayGroup@biological_replicates)
        testSamplesSize<-length(testAssayGroup@biological_replicates)
      }
    }
    
    readRDS(paste0(opt$inputdir,"/",accession,"_",contrast,"_fit.rds"))->fit
    
    setkey(annot,DesignElementAccession)
    as.data.table(topTable(fit, confint=TRUE, number=Inf))->expTTable
    expTTable[,variance:=((CI.R-CI.L)/3.92)^2,]
    setkey(expTTable,designElements)
    expTTable$Accession<-accession
    expTTable$Contrast<-contrast
    expTTable$refSampleSize<-referenceSamplesSize
    expTTable$testSampleSize<-testSamplesSize
    setnames(expTTable, "adj.P.Val", "padj")
    
    columns<-c('Accession', 'Contrast', 'Gene ID', 'Gene Name', 'logFC', 'variance', 'padj', 'refSampleSize', 'testSampleSize')
    results[[paste0(accession,"_",contrast)]]<-expTTable[annot, on=c(designElements="DesignElementAccession")][, columns, with=FALSE]
    summary(results[[paste0(accession,"_",contrast)]])
  }
}
print("Done microarrays\n")

readRDS(paste0(opt$inputdir,"/rna_seq_diff.rds"))->rna_seq_diff

for( accession in rna_seq_diff) {
  
  experimentConfig <- parseAtlasConfig( paste0( opt$inputdir,"/",accession,"-configuration.xml" ))
  analytics <- experimentConfig$allAnalytics[[1]]
  expContrasts <- atlas_contrasts( analytics )
  expAssayGroups <- assay_groups( analytics )
  
  readRDS(paste0(opt$inputdir,"/",accession,".annot.rds"))->annot
  setkey(annot, `Gene ID`)
  
  for( contrast in unlist(strsplit(x=datasets$contrast[datasets$accession==accession], split = ",")) ) {
    
    for( expContrast in expContrasts ) {
    # Get the sample sizes for the contrast at hand.
      if( expContrast@contrast_id == contrast ) {
        refAssayGroupID <- reference_assay_group_id( expContrast )
        testAssayGroupID <- test_assay_group_id( expContrast )
        
        refAssayGroup <- expAssayGroups[[ refAssayGroupID ]]
        testAssayGroup <- expAssayGroups[[ testAssayGroupID ]]
        
        referenceSamplesSize<-length(refAssayGroup@biological_replicates)
        testSamplesSize<-length(testAssayGroup@biological_replicates)
      }
    }
    
    readRDS(paste0(opt$inputdir,"/",accession,"_",contrast,"_deseq.rds"))->deseq2Res
    data.table(Accession=accession, Contrast=contrast, `Gene ID`=deseq2Res@rownames, 
               logFC=deseq2Res@listData$log2FoldChange, 
               variance=((deseq2Res@listData$lfcSE)^2*(referenceSamplesSize+testSamplesSize)),
               refSampleSize=referenceSamplesSize,
               testSampleSize=testSamplesSize,
               padj=deseq2Res@listData$padj
               )->expTTable
    setkey(expTTable,`Gene ID`)
    columns<-c('Accession', 'Contrast', 'Gene ID', 'Gene Name', 'logFC', 'variance', 'padj', 'refSampleSize', 'testSampleSize')
    results[[paste0(accession,"_",contrast)]]<-expTTable[annot, on=c(`Gene ID`="Gene ID")][!is.na(padj), columns, with=FALSE]
    summary(results[[paste0(accession,"_",contrast)]])
  }
  
}

summary(results)

saveRDS(results, file = paste0(opt$outdir,"/","expression_tables.rds"))

results<-rbindlist(results)

fwrite(file=paste0(opt$outdir,"/","merged_logfc_variance.tsv"), 
       # merge the tables and produce desired output.
       results,sep='\t')
