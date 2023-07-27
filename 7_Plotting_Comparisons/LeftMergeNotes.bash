module load StdEnv/2020
module load r/4.2.1
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/4.2.1/

cp /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Blast_Output/Blast_ref_output/cleaned_blast_ref/cleaned_blast_ref_total_subtract_W_C_Mat_Sen.out .


Rscript TotalMerge.R
