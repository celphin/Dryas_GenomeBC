#########################################
#In progress notes (July 2023)
#########################################
module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/
R 

#work through Dryas_Merge_Notes.R


cut -f1,14 interproscan_dryas_full.tsv| grep "GO" >  dryas_goterm_file.tsv 

Rscript Combine_Dryas_Blast.R
