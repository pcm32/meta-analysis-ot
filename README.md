# meta-analysis-ot

This repo holds scripts for running meta-analysis with EBI Expression Atlas Microarray and RNA-Seq differential studies. 
Given a list of Expression Atlas accessions and desired contrasts for them, these scripts will fetch data and metadata 
from Expression Atlas's ftp site at EBI, do some preprocessing, call differential expression and the merge results for usage
with meta-analysis packages.

# Define desired contrasts

For each experiment of choice, explore the contrasts by executing

```
inspectContrasts.R -a <accession>
```

This will provide the list of contrast for that experiment, for instance:


| Accession | C.ID | Size | C. Name |
|-----------|------|------|---------|
| E-GEOD-16879 | g13_g1 | 6:7 | 'after first infliximab treatment; no response to infliximab treatment; Crohn's disease' vs 'control' in 'colon' |
| E-GEOD-16879 | g14_g2 | 6:10 | 'after first infliximab treatment; no response to infliximab treatment; Crohn's disease' vs 'control' in 'ileum' |
| E-GEOD-16879 | g13_g3 | 6:16 | 'after first infliximab treatment; no response to infliximab treatment; ulcerative colitis' vs 'control' in 'colon' |

For each experiment/accession copy the accession and contrasts ids for the desired contrasts, to a tab separated file that 
looks like (make sure that columns are named exactly as shown, more than one contrast per experiment can be given, separated by commas):

| accession | contrast |
|-----------|----------|
| E-GEOD-16879 | g13_g1,g14_g2 |
| E-METAB-2175 | g4_g1 |

and save the file (referred to as datasets file later on).

# Run

```
runProcess.sh <path-to-datasets-file> <path-to-output-directory>
```

A number of intermediate files are generated, merged results are left in `merged_logfc_variance.tsv` in the output directory.
This file will include, per (experiment, contrast, gene): logFC, variance, adj p-value and sample sizes (reference and test)

A conda package and container will be available shortly. In the future this functionality will probably move to the 
Expression Atlas bioconductor package.
