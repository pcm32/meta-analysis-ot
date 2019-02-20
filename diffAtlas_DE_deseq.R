#!/usr/bin/env Rscript

# diffAtlas_DE_deseq2.R
# RNA-seq differential expression statistics computation for Expression Atlas.
#

suppressMessages( library( DESeq2 ) )

suppressMessages( library( ExpressionAtlasInternal ) )

# diffAtlas_DE_deseq2
# - Differential expression analysis (2-group comparison) using DESeq2.
# Arguments:
# 	expAcc <- ArrayExpress accession of experiment.
# 	atlasProcessingDirectory <- path to Atlas processing directory.
# 	countsMatrixFile <- path to counts matrix from iRAP.
diffAtlas_DE_deseq2 <- function( expAcc, atlasProcessingDirectory, countsMatrixFile ) {

	e <- try( {

		# Make the config filename.
		xmlConfigFilename <- paste( expAcc, "-configuration.xml", sep = "" )
		xmlConfigFilename <- file.path( atlasProcessingDirectory, xmlConfigFilename )

		# Parse the config.
		cat( paste( "Reading XML from", xmlConfigFilename, "...\n" ) )
		experimentConfig <- parseAtlasConfig( xmlConfigFilename )
		cat( "Successfully read XML config.\n" )

		# Get the experiment type.
		experimentType <- experimentConfig$experimentType

		cat( paste( "Experiment type is", experimentType, "\n" ) )

		# Check this is an RNA-seq experiment.
		if( !grepl( "rnaseq", experimentType ) ) {
			stop( paste(
						"Experiment type",
						experimentType,
						"does not look like an RNA-seq experiment. Cannot continue."
			) )
		}

		# Get the list of analytics objects from the config.
		allAnalytics <- experimentConfig$allAnalytics

		# There should only be one analytics object for an RNA-seq experiment,
		# so make sure this is the case.
		if( length( allAnalytics ) > 1 ) {
			stop( "More than one analytics element found in XML. Cannot continue." )
		}
		
		# Get the analytics object.
		analytics <- allAnalytics[[ 1 ]]

		if( platform( analytics ) != "rnaseq" ) {
			stop( paste( 
						"Don't know what to do with analytics of type",
						platform( analytics )
			) )
		}

		cat( paste( "Reading raw counts from", countsMatrixFile, "...\n" ) )

		# Read in the raw counts.
		countsMatrix <- read.delim(
								   countsMatrixFile,
								   header = TRUE,
								   stringsAsFactors = FALSE,
								   row.names = 1
								   )
		
		cat( "Successfully read raw counts.\n" )
		
		# Get the contrasts, assay groups, and batch effects.
		expContrasts <- atlas_contrasts( analytics )
		expAssayGroups <- assay_groups( analytics )

		cat( paste( "Found", length( expContrasts ), "contrasts and", length( expAssayGroups ), "assay groups.\n" ) )

		# Go through the contrasts...
		invisible( sapply( expContrasts, function( expContrast ) {

						# Get the contrast info.
			contrastID <- contrast_id( expContrast )
			contrastName <- contrast_name( expContrast )
			refAssayGroupID <- reference_assay_group_id( expContrast )
			testAssayGroupID <- test_assay_group_id( expContrast )
			contrastBatchEffects <- batch_effects( expContrast )

			if( Sys.getenv("ONLY_CONTRAST_ID") != "" && !(contrastID %in% unlist(strsplit(Sys.getenv("ONLY_CONTRAST_ID"), split="," )) ) ) {
				cat( paste0("Skipping contrast ",contrastID," ",contrastName,"\n") )
				return(0)
		        }

			cat( paste( 
					   "Processing contrast",
					   contrastID,
					   "(",
					   contrastName,
					   ")...\n"
			) )

			# Get the two assay groups for this contrast as a list.
			contrastAssayGroups <- select_assay_groups_for_contrast( 
																 expAssayGroups, 
																 refAssayGroupID, 
																 testAssayGroupID 
																 )
			
			cat( "Creating biological replicate annotations...\n" )

			# Create the data frame containing annotations for each
			# biological replicate, grouping technical replicates.
			bioRepAnnotations <- make_biorep_annotations( contrastAssayGroups, contrastBatchEffects )

			cat( "Annotations created successfully.\n" )

			cat( "Subsetting normalized data for this contrast...\n" )

			# Next, subset the normalized data frame so that it
			# contains only the columns of data for this contrast.
			countsForContrast <- countsMatrix[ , bioRepAnnotations$AssayName ]
			
			cat( "Subsetting successful.\n" )

			cat( "Checking for technical replicates...\n" )

			# Get the technical replicate group IDs -- these are the
			# ones in the BioRepName column that are duplicated.
			techRepGroupIDs <- unique( bioRepAnnotations$BioRepName[ which( duplicated( bioRepAnnotations$BioRepName ) ) ] )

			# If there are any technical replicates...
			if( length( techRepGroupIDs ) > 0 ) {
				
				cat( "Technical replicates found. Calculating sum of each set of technical replicates...\n" )

				# Replace the columns for technical replicates with their sums in the normalized data.
				countsForContrast <- add_tech_rep_sums( countsForContrast, bioRepAnnotations, techRepGroupIDs )

				cat( "Technical replicate averaging successful.\n" )
			
			} else {

				cat( "No technical replicates found.\n" )
			}
			
			cat( "Setting up biological replicate annotations for DESeqDataSet object...\n" )

			# Remove the AssayName column from the annotations data frame now we no longer need it.
			bioRepAnnotations$AssayName <- NULL

			# Remove duplicated technical replicate rows, if any.
			duplicatedRowIndices <- which( duplicated( bioRepAnnotations$BioRepName ) )
			if( length( duplicatedRowIndices ) > 0 ) {
				bioRepAnnotations <- bioRepAnnotations[ -duplicatedRowIndices, ]
			}

			# Add the biological replicate names as the row names.
			rownames( bioRepAnnotations ) <- bioRepAnnotations$BioRepName
			
			# Remove the bio rep names column as we don't need it.
			bioRepAnnotations$BioRepName <- NULL

			# Make the column names R-safe.
			colnames( bioRepAnnotations ) <- make.names( colnames( bioRepAnnotations ) )

			cat( "Annotations modified successfully.\n" )

			cat( "Re-ordering this contrast's counts matrix...\n" )

			# Put the columns of the counts matrix into the same order as the rows in bioRepAnnotations.
			countsForContrast <- countsForContrast[ , rownames( bioRepAnnotations ) ]

			cat( "Counts matrix re-ordered successfully.\n" )

			# Convert countsMatrix (data frame) into a matrix.
			countsForContrast <- as.matrix( countsForContrast )

			cat( "Checking for batch effects...\n" )

			# Make the formula to create the design matrix, adding
			# batch effects if necessary.
			# If we have batch effects...
			if( length( contrastBatchEffects ) > 0 ) {

				cat( paste( length( contrastBatchEffects ), "batch effects found.\n" ) )
				
				cat( "Making R-safe batch effect names...\n" )

				# Get the batch effect names -- these should be the same as
				# the column headings in the ExpressionSet's info.
				batchEffectNames <- make.names( sapply( contrastBatchEffects, function( batchEffect ) { effect_name( batchEffect ) } ) )
				
				cat( "Successfully created R-safe batch effect names.\n" )

				cat( "Creating formula for design matrix...\n" )

				# Create the formula string.
				formulaString <- paste( "~", paste( batchEffectNames, collapse = " + " ), "+ groups" )

				cat( paste( "Formula is:", formulaString, "\n" ) )
			
			} else {
				
				cat( "No batch effects found.\n" )

				# If we don't have any batch effects, create the formula string without them.
				formulaString <- "~ groups"

				cat( paste( "Formula is:", formulaString, "\n" ) )
			}

			# Now we have all the info we need, we can create the DESeqDataSet object.
			cat( "Creating DESeqDataSet object...\n" )
			
			deseqDataSet <- DESeqDataSetFromMatrix(
												   countData = countsForContrast,
												   colData = bioRepAnnotations,
												   design = as.formula( formulaString )
												   )

			cat( "DESeqDataSet created successfully.\n" )

			cat( "Running differential expression analysis...\n" )

			deseqDataSet <- DESeq( deseqDataSet )

			cat( "Differential expression analysis successful.\n" )

			cat( "Performing independent filtering and creating results table...\n" )

			res <- results( deseqDataSet, contrast = c( "groups", "test", "ref" ) )

			if( Sys.getenv("DESEQ_RDS_OUTPATH") != "" ) {
			  saveRDS(res, file=file.path( Sys.getenv( "DESEQ_RDS_OUTPATH" ), paste0(expAcc,"_",contrastID,"_deseq.rds") ))	
			}	

			# Make sure that all the columns of res are numeric. If not maybe
			# something went wrong; non-numeric values cause problems later.
			if( !all( sapply( res, is.numeric ) ) ) {
				stop( "Non-numeric values found in DESeq results. Cannot continue." )
			}

			cat( "Independent filtering and results table creation successful.\n" )

			
			# Here need to select wanted columns for outut files.
			# Stats results:
			contrastResults <- data.frame(
										  id = rownames( res ),
										  log2FoldChange = res$log2FoldChange,
										  padj = res$padj
										  )
			# MvA plot data:
			plotData <- data.frame(
								   geneID = rownames( res ),
								   avgExpr = res$baseMean,
								   logFC = res$log2FoldChange,
								   adjPval = res$padj
								   )

			# Create filename to write to.
			# Stats results:
			resFile <- paste( expAcc, contrastID, "analytics", "tsv", sep = "." )
			resFile <- file.path( Sys.getenv( "HOME" ), "tmp", resFile )
			# Data for MvA plot:
			plotDataFile <- paste( expAcc, contrastID, "plotdata", "tsv", sep = "." )
			plotDataFile <- file.path( Sys.getenv( "HOME" ), "tmp", plotDataFile )

			cat( paste( "Writing differential expression results to", resFile, "...\n" ) )

			write.table( contrastResults, file=resFile, row.names=FALSE, quote=FALSE, sep="\t" )

			cat( paste( "Results written successfully.\n" ) )

			cat( paste( "Writing data for MvA plot to", plotDataFile, "...\n" ) )

			write.table( plotData, file=plotDataFile, row.names=FALSE, quote=FALSE, sep="\t" )

			cat( "Plot data written successfully.\n" )

			cat( paste( "Successfully completed differential expression analysis for contrast", contrastID, "\n" ) )

		} ) )

	} ) # try
	
	# Die if we got an error.
	if( class( e ) == "try-error" ) {
		stop( e )
	}
}


# select_assay_groups_for_contrast
# 	- Given a list of experiment assay groups, a test ID and a reference ID,
# 	return a list containing just the requested test and reference assay
# 	groups.
select_assay_groups_for_contrast <- function( expAssayGroups, refAssayGroupID, testAssayGroupID ) {

	# Get the assay groups.
	refAssayGroup <- expAssayGroups[[ refAssayGroupID ]]
	testAssayGroup <- expAssayGroups[[ testAssayGroupID ]]

	# Return them in a list.
	return( 
		   list( 
				 refAssayGroup = refAssayGroup,
				 testAssayGroup = testAssayGroup
				 ) 
	)
}


# make_biorep_annotations
# 	- Given the contrast's assay groups, the list of batch effects (if any),
# 	create a data frame representing the biological replicates and their
# 	corresponding assay names and which group (test or ref) they belong to.
make_biorep_annotations <- function( contrastAssayGroups, contrastBatchEffects ) {

	refAssaysToBioReps <- make_assays_to_bioreps_df( biological_replicates( contrastAssayGroups$refAssayGroup ) )
	testAssaysToBioReps <- make_assays_to_bioreps_df( biological_replicates( contrastAssayGroups$testAssayGroup ) )
	
	# Mappings of assay names to biological replicate names. Biological
	# replicate names are most commonly the same as the assay names except in
	# the case of technical replicates, in which case the biological replicate
	# name is the technical replicate group ID.
	allAssaysToBioReps <- rbind( refAssaysToBioReps, testAssaysToBioReps )

	# Remove the row names as they're not useful.
	rownames( allAssaysToBioReps ) <- NULL

	# Get the contrast assay names to subset the batch effect lists with.
	contrastAssayNames <- allAssaysToBioReps$AssayName
	
	# Make the column for the group (i.e. test or reference) that the assays
	# belong to.
	groupsCol <- c( rep( "ref", nrow( refAssaysToBioReps ) ), rep( "test", nrow( testAssaysToBioReps ) ) )

	# Add the groups column to the data frame.
	allAssaysToBioReps$groups <- as.factor( groupsCol )

	if( length( contrastBatchEffects ) > 0 ) {
	
		# Get a list of data frames mapping assay names to batch effects.
		allBatchEffectDfs <- lapply( contrastBatchEffects, function( batchEffect ) {
			batch_effect_to_df( batchEffect, contrastAssayNames )
		})
		
		# Combine all batch effect data frames into one.
        batchEffectsDf <- do.call( "cbind", allBatchEffectDfs )

		# Sort the rows of batchEffectsDf in the same order as the assay names in
		# allAssaysToBioReps.
		batchEffectsDf <- batchEffectsDf[ allAssaysToBioReps$AssayName , , drop = FALSE ]

		# Combine the batch effects data frame and the assays-to-bioreps data frame
		# to make the final data frame.
		return( cbind( allAssaysToBioReps, batchEffectsDf ) )
	
	} else {
		
		return( allAssaysToBioReps )
	}
}


# batch_effect_to_df
# 	- Given a batch effect object and a vector of assay names, create a data
# 	frame mapping the batch names to the corresponding assay names.
batch_effect_to_df <- function( batchEffect, contrastAssayNames ) {

	assayBatches <- batches( batchEffect )
	effectName <- effect_name( batchEffect )

	# Get a list of data frames with the assays from each batch.
	batchDataFrames <- lapply( names( assayBatches ), function( batchName ) {
		
		assayNames <- assayBatches[[ batchName ]]

		data.frame( AssayNames = assayNames, batch = batchName, stringsAsFactors = FALSE )
	} )

	# Merge the data frames in to one data frame.
	allBatches <- do.call( "rbind", batchDataFrames )

	# Make the columns into factors.
	allBatches <- data.frame( apply( allBatches, 2, function( x ) { as.factor( x ) } ) )

	# Change the column heading from "batch" to the actual effect name.
	colnames( allBatches ) <- c( "AssayName", effectName )

	# Subset only the contrast's assay names.
	allBatches <- allBatches[ which( allBatches$AssayName %in% contrastAssayNames ), ]

	# Re-order the rows based on the assay names.
	allBatches <- allBatches[ order( allBatches$AssayName ), ]

	rownames( allBatches ) <- allBatches$AssayName
	allBatches$AssayName <- NULL

	return( allBatches )
}


# make_assays_to_bioreps_df
# 	- Given a list of biological replicate objects create a data frame mapping
# 	assay names to their biological replicate names. These will either be the
# 	same as the assay name (no technical replicates), or the technical
# 	replicate group ID (technical replicates).
make_assays_to_bioreps_df <- function( bioReps ) {

	assaysToBioRepsList <- lapply( bioReps, function( bioRep ) {

		assayNames <- biorep_assay_names( bioRep )
		techRepId <- technical_replicate_id( bioRep )

		if( length( techRepId ) > 0 ) {
			
			data.frame( AssayName = assayNames, BioRepName = techRepId, stringsAsFactors = FALSE )

		} else {
			
			data.frame( AssayName = assayNames, BioRepName = assayNames, stringsAsFactors = FALSE )
		}
	})

	assaysToBioRepsDf <- do.call( "rbind", assaysToBioRepsList )

	return( assaysToBioRepsDf )
}


# add_tech_rep_sums
# 	- Given a data frame of expression values, the data frame of biological
# 	replicate annotations, technical replicates group IDs, calculate the sum
# 	for each set of technical replicate columns in the data frame, and replace
# 	the original columns with the new column(s) containing the sums.
add_tech_rep_sums <- function( expressionData, bioRepAnnotations, techRepGroupIDs ) {

	# Create a list with one element per column of sums.
	techRepSumsList <- lapply( techRepGroupIDs, function( techRepGroupID ) {

		# Get the row indices for the technical replicates.
		techRepRows <- which( bioRepAnnotations$BioRepName == techRepGroupID )

		# Get the assay names on these rows.
		techRepAssayNames <- bioRepAnnotations$AssayName[ techRepRows ]

		# Get the columns of normalized data for these assays.
		techRepColumns <- expressionData[ , techRepAssayNames ]

		# Get the sum of each row.
		techRepSums <- apply( techRepColumns, 1, function( x ) { sum( x ) } )

		techRepSumsDF <- structure( data.frame( techRepSums ), names = techRepGroupID )
	} )

	# Turn the list into a data frame.
	techRepSumsDF <- do.call( "cbind", techRepSumsList )

	# Next need to remove the columns for the technical replicate assays from
	# the normalized data frame.
	# First get all the assay names we want to remove.
	allTechRepRows <- which( bioRepAnnotations$BioRepName %in% techRepGroupIDs )

	allTechRepAssays <- bioRepAnnotations$AssayName[ allTechRepRows ]

	# Remove the columns for these assay names from the normalized data frame.
	expressionData <- expressionData[ , !( colnames( expressionData ) %in% allTechRepAssays ), drop = FALSE ]

	# Add the columns of sums instead.
	expressionData <- cbind( expressionData, techRepSumsDF )

	return( expressionData )
}


###################
# RUN MAIN FUNCTION
###################

# Run with arguments if there are any, otherwise don't do anything. Having this
# lets us source this file in an R session and load all the functions so we can
# run them separately if desired.
args <- commandArgs( TRUE )
if( length( args ) > 0) {
	do.call( diffAtlas_DE_deseq2, as.list( args ) )
}
