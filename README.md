# Delta_Smelt_Eigenetics_Transcriptomics

## CTMs

*Scripts*: Fig1_CTM_Analyses.R

*Input Files*: 2021_CTM_data_Trial1-120.txt, Family_DIs.txt, all_pedigree_meta_DI.txt

*Description*: This script contains code for creating Fig 1 results as well as running the ANOVA for estimating differences in CTMs between rearing temperature and domestication index. This input file, 2021_CTM_data_Trial1-120.txt, contains experimental meta data for each fish and it's corresponding CTM. This input file, Family_DIs.txt, contains information on each family's numerical domestication index. This input file, all_pedigree_meta_DI.txt, contains pedigree information for each fish, including it's parents IDs and family ID.

## RNA-Seq

*Scripts*: Voom_delta_smelt.R

*Input files*: RNA_metadata2.txt, all_quant.sf

*Description*: This script contains the code for filtering and normalizing the read count matrix. It then performs an analysis using the limma voom package and creates heatmaps of results. The input file, RNA_metadata2.txt, contains the metadata for each sequenced sample. This input file, all_quant.sf, is the read count matrix for each sequenced sample (counts are Transcripts Per Million, TPM).

For generation of the read count matrix (all_quant.sf), see: https://github.com/JoannaGriffiths/RNASeq-snakemake-pipeline

## WGBS

*Scripts*: Delta-smelt-Bisulfite-pipeline.Rmd (View the [online version](https://htmlpreview.github.io/?https://github.com/JoannaGriffiths/Delta_Smelt_Eigenetics_Transcriptomics/blob/main/WGBS/Delta-smelt-Bisulfite-pipeline.html)), DMR_Interaction_test.R, DMR_DEG_correlation.R

*Description*: The Delta-smelt-Bisulfite-pipeline R markdown script contains code for quality control of fastq files, aligning reads to the reference genome, and calling methylated sites. This scripts also contains the code for running differential methylation analysis. The script, DMR_DEG_correlation.R identifies genes that are both differentially methylated and differentially expressed and whether percent methylated is significantly correlated to expression levels. The script, DMR_Interaction_test.R, determines whether there is a statistically significant temperature by domestication index interaction effect for DMRs. 

## Functional Enrichment

*Scripts*: Generate_input_files_for_GO_MWU.R, Run_GO_MWU.R, DMRichR_func_enrich.R

*Description*: The script, Generate_input_files_for_GO_MWU.R, coverts DEGs and DMRs output into the correct format for input of functional enrichment analysis. The script, Run_GO_MWU.R, then runs the functional enrichment analysis using a slightly modified code from Matz GitHib). This script, DMRichR_func_enrich.R, identifies genes that are found within (or overlapping with) DMRs. These genes are then used for input for functional enrichment analyses. 
