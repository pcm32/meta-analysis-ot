#!/bin/env bash

working_directory=${2}
datasets=${1}

get_configs.R -d $datasets -o $working_directory

procure_data_microarray.R -i $working_directory -o $working_directory
procure_data_rnaseq_diff.R -i $working_directory -o $working_directory

while IFS='' read -r accession || [[ -n "$accession" ]]; do
	export ONLY_CONTRAST_ID=$( grep $accession $datasets | awk -F'\t' '{ print $2 }' )
	export DESEQ_RDS_OUTPATH=$working_directory
	diffAtlas_DE_deseq.R $accession $working_directory $working_directory/$accession\-raw-counts.tsv.undecorated 
done < $working_directory/rna_seq_diff.txt

while IFS='' read -r accession || [[ -n "$accession" ]]; do
	export ONLY_CONTRAST_ID=$( grep $accession $datasets | awk -F'\t' '{ print $2 }' )
	export FIT_RDS_OUTPATH=$working_directory
	diffAtlas_DE_limma.R $accession $working_directory 
done < $working_directory/microarray.txt

get_merged_logfc_variance.R -d $datasets -i $working_directory -o $working_directory
