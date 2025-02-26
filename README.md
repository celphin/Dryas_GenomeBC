# Dryas_GenomeBC
Differential DNA methylation and RNAseq analyses
 
 ## 1_RNAseq - mRNA
* RSEM for mapping to transcriptome 
* EdgeR for transcript count normalization, differential expression analyses and plotting heat maps

 ## 2_Methylseq
* nf-core Methylseq for WGBS mapping and calling methylated C's 
* Run on wild, phenology, and seedling data
 
 ## 3_DMRs
* Metilene and Methylkit for detecting DMRs
* Run for detecting various Wild, Seedling, and Phenology DMRS
* Scripts for preparing files for Metilene and running metilene

## 4_Bedgraphs_Intersections
* Bedtools for comparisons
* Crosses over DMRs, to find intersecting or non-intersecting DMRS

## 5_Blast
* Script for extending DMRS, and mapping to fasta file
* Notes on Blasting DMRS to reference and to NCBI

## 6_Interproscan
* Scripts for running interproscan on a set of files
* Script on running interproscan on divided subseq files
* Script for formatting file to run on interproscan
* Notes on running interproscan on both Methylation and RNA data, and on comparing/combining the two

## 7_Gene_ontology
* ErmineJ for GO term enrichment

## 8_Plotting_Comparisons
* Several scripts for plotting
* General notes on heatmap plotting (RNA and methylation)
* Notes on inheritence plotting

## 9_SNP_calling
* BS-SNP for SNP calling from bisulphite treated DNA
* Basic popgen analyses

## 10_Phenotypes
* Comparing phenotypes and traits of the warmed and control plants

 ## Old_Files_Semi_Sorted - Temporary
* Notes on original runs of certain anlyses + analyses not yet used/completed

## Unused scripts
* Phenogram script for interproscan + annotation file to format to phenogram
* gff3 to cds conversion

