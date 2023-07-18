#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-1:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120000m

module load bedtools/2.30.0
module load bedops/2.4.41

gff3_file=HAN412_Eugene_curated_v1_1.gff3
full_fasta=Ha412HOv2.0-20181130.fasta
output_name=Ha412HOv2.0_CDS.fasta

#Converts to bed 
srun gff2bed < $gff3_file > annotation.bed
#Optional (required for files in Sunflower database)
grep "CDS" annotation.bed > cleaned_annotation.bed
#Get fasta based on bed file
srun bedtools getfasta -fi $full_fasta -bed cleaned_annotation.bed -fo $output_name


