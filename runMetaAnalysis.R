#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(MetaVolcanoR))
# these will be remove once metavolcano resolves dependency issues on install
# and uses direct method calls from packages.
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(metafor))
suppressPackageStartupMessages(require(plotly))
suppressPackageStartupMessages(require(parallel))

option_list = list(
  make_option(
    c("-r", "--rds_expression_tables"),
    action = "store",
    default = NA,
    type = 'character',
    help = "Path to RDS object with list of data.tables for each experiment-contrast"
  ),
  make_option(
    c("-o", "--outdir"),
    help = "Output directory"
  ),
  make_option(
    c("-n", "--n_cores"),
    help = "Number of cores to use for metafor",
    default = 1,
    type = 'integer'
  ),
  make_option(
   c("-p", "--adj_pvalue_threshold"),
   help = "Adjusted p-value threshold, float between 0 and 1",
   default = 0.05,
   type = 'numeric'
  ),
  make_option(
    c("-l", "--logFC_threshold"),
    help = "LogFC threshold above 0",
    default = 1.0,
    type = 'numeric'
  ),
  make_option(
    c("-j", "--job_name"),
    help = "Job name, for metafor file naming",
    default = "meta-analysis_metavolcano",
    type = 'character'
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

readRDS(opt$rds_expression_tables)->expressionTables

meta_degs_metafor <- do.metafor(geo2r_res=expressionTables, pcriteria='padj', 
                                genenamecol='Gene Name', geneidcol='Gene ID', 
                                foldchangecol='logFC', 
                                pvalue=opt$adj_pvalue_threshold, logfc=opt$logFC_threshold, ncores=opt$n_cores,
                                collaps=FALSE, jobname=opt$job_name, draw = TRUE,
                                llcol='CI.L', rlcol='CI.R', vcol='variance',
                                outputfolder=opt$outdir, cvar=TRUE)

write.table(meta_degs_metafor, paste0(outputfolder=opt$outdir, "metafor_degs_", jobname=opt$job_name, "_summarizing.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
