###############################################################################
#Find Parent DMRS with Metilene: 
    #Warming vs Control (Site Specific)
#Based on Notes2a_WildGBS_MetileneDMRS_July2022.txt
################################################################################
#Requires:
#methylseq output, unzipped
#metile_prep.sh:
    #adjusts all the files names based on categories,
    #prepares them in the correct format for metilene
    #Adjust categories, dir names, renaming convention
#metilene_run.sh
    #runs metilene for specified parameters
    #Adjust categories/dir names
#metilene_filter_qval.sh
    #Filters metilene output
    #Adjust categories/dir names
#NOTE:
    #Save all in metilene files for later 
##################################################################################

# copy old ref genome bedGraph data to scratch

# aim to compare the high lat sites with low lat sites and find consistient differences

# which of these overlap with the W/C DMRs?







###############################################
#Nunavut:
mkdir Nunavut_Metilene
cd Nunavut_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .

tmux new-session -s Nunavut_DMRS
tmux attach-session -t Nunavut_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account def-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0

#Nunavut adjustments
#h1,h2: 
    #h1="W", h2="C"
#Input directories:
    #input_dir=Nunavut_${h1}_${h2}_input_files 
    #in_metilene="Nunavut_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Copy Loop:
    #cp ${methylseq_output_dir}/*FERT*.deduplicated.bedGraph ./${input_dir}
    #cp ${methylseq_output_dir}/*DRY*.deduplicated.bedGraph ./${input_dir}
    #cp ${methylseq_output_dir}/*MEAD*.deduplicated.bedGraph ./${input_dir}
    #cp ${methylseq_output_dir}/*WILL*.deduplicated.bedGraph ./${input_dir}
    #cp ${methylseq_output_dir}/*CASS*.deduplicated.bedGraph ./${input_dir}
#Comment out rename for loops 


sh metilene_prep.sh 

#------------------------------------------------------------------

#Nunavut specific adjustments,
#h1,h2: 
    #h1="W", h2="C"
#in_metilene="Nunavut_metilene_"$h1"_"$h2".input"
#input_dir=Nunavut_${h1}_${h2}_input_files 
#outputname=Nunavut_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
#parameters: maxdist, mincpgs, mindiff

#metilene:
# params: maxdist, mincpgs, mindiff
module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7
#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value
#h1="W", h2="C"
#input_dir="/home/msandler/scratch/Wild_Metilene/Nunavut_Metilene/Nunavut_${h1}_${h2}_input_files"
#outputname=Nunavut_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
module load nixpkgs/16.09 
module load gcc/7.3.0
module load r/3.6.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/3.6/
sh metilene_filter_qval.sh 70 5 4 0.9 0.001
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001





####################################################################################
#Svalbard:
mkdir Svalbard_Metilene
cd Svalbard_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .

tmux new-session -s Svalbard_DMRS
tmux attach-session -t Svalbard_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account def-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0


#Svalbard adjustments
#h1,h2: 
    #h1="W", h2="C"
#Input directories:
    #input_dir=Svalbard_${h1}_${h2}_input_files 
    #in_metilene="Svalbard_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Copy Loop:
    #cp ${methylseq_output_dir}/*Sval*.deduplicated.bedGraph ./${input_dir}
#Comment out rename for loops 


sh metilene_prep.sh 

#---------------------------------------------

#Svalbard specific adjustments,
#h1,h2: 
    #h1="W", h2="C"
#in_metilene="Svalbard_metilene_"$h1"_"$h2".input"
#input_dir=Svalbard_${h1}_${h2}_input_files 
#outputname=Svalbard_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
#parameters: maxdist, mincpgs, mindiff

#metilene:
# params: maxdist, mincpgs, mindiff
module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7


#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value
#h1="W", h2="C"
#input_dir="/home/msandler/scratch/Wild_Metilene/Svalbard_Metilene/Svalbard_${h1}_${h2}_input_files"
#outputname=Svalbard_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
module load nixpkgs/16.09 
module load gcc/7.3.0
module load r/3.6.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/3.6/
sh metilene_filter_qval.sh 70 5 4 0.9 0.001
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001

###################################################################################
#Sweden:
mkdir Sweden_Metilene
cd Sweden_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .

tmux new-session -s Sweden_DMRS
tmux attach-session -t Sweden_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account def-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0


#Sweden adjustments
#h1,h2: 
    #h1="W", h2="C"
#Input directories:
    #input_dir=Sweden_${h1}_${h2}_input_files 
    #in_metilene="Sweden_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Copy Loop:
    #cp ${methylseq_output_dir}/*LAT*.deduplicated.bedGraph ./${input_dir}
#Comment out rename for loops 


sh metilene_prep.sh 

#---------------------------------------------

#Sweden  specific adjustments,
#h1,h2: 
    #h1="W", h2="C"
#in_metilene="Sweden_metilene_"$h1"_"$h2".input"
#input_dir=Sweden_${h1}_${h2}_input_files 
#outputname=Sweden_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
#parameters: maxdist, mincpgs, mindiff

#metilene:
# params: maxdist, mincpgs, mindiff
module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7

#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value
#h1="W", h2="C"
#input_dir="/home/msandler/scratch/Wild_Metilene/Sweden_Metilene/Sweden_${h1}_${h2}_input_files"
#outputname=Sweden_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
module load nixpkgs/16.09 
module load gcc/7.3.0
module load r/3.6.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/3.6/
sh metilene_filter_qval.sh 70 5 4 0.9 0.001
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001






#################################################################################
#Alaska 
mkdir Alaska_Metilene
cd Alaska_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .

tmux new-session -s Sweden_DMRS
tmux attach-session -t Sweden_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account def-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0


#Alaska adjustments
#h1,h2: 
    #h1="W", h2="C"
#Input directories:
    #input_dir=Alaska_${h1}_${h2}_input_files 
    #in_metilene="Alaska_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Copy Loop:
    #cp ${methylseq_output_dir}/*ALAS*.deduplicated.bedGraph ./${input_dir}
#Comment out rename for loops 


sh metilene_prep.sh 

#----------------------------------------
#Still not run
#Alaska specific adjustments,
#h1,h2: 
    #h1="W", h2="C"
#in_metilene="Alaska_metilene_"$h1"_"$h2".input"
#input_dir=Alaska_${h1}_${h2}_input_files 
#outputname=Alaska_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
#parameters: maxdist, mincpgs, mindiff

#metilene:
# params: maxdist, mincpgs, mindiff
module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7
#-------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value
#h1="W", h2="C"
#input_dir="/home/msandler/scratch/Wild_Metilene/Alaska_Metilene/Alaska_${h1}_${h2}_input_files"
#outputname=Alaska_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
module load nixpkgs/16.09 
module load gcc/7.3.0
module load r/3.6.0
module load gdal
module load udunits
module load python
export R_LIBS_USER=/home/msandler/R/x86_64-pc-linux-gnu-library/3.6/
sh metilene_filter_qval.sh 70 5 4 0.9 0.001
sh metilene_filter_qval.sh 150 5 4 0.9 1e-5
sh metilene_filter_qval.sh 70 5 0.7 0.7 0.001

##############################################################
cd ..
#Copy betools files to directory
cp Nunavut_Metilene/Nunavut_W_C_70_5_4_0.9_qval.0.001.bedgraph Nunavut_W_C.bedGraph
cp Svalbard_Metilene/Svalbard_W_C_70_5_4_0.9_qval.0.001.bedgraph Svalbard_W_C.bedGraph
cp Sweden_Metilene/Sweden_W_C_70_5_4_0.9_qval.0.001.bedgraph Sweden_W_C.bedGraph
cp Alaska_Metilene/Alaska_W_C_70_5_4_0.9_qval.0.001.bedgraph Alaska_W_C.bedGraph

#Bedtools intersect all 4
module load bedtools/2.30.0
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALL_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Nunavut
bedtools intersect -u -a Svalbard_W_C.bedGraph -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > SVAL_SWED_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Svalbard
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Sweden_W_C.bedGraph | bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALEX_SWED_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Sweden
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Alaska_W_C.bedGraph > ALEX_SVAL_ALAS_Sites_intersect_DMRs.bedgraph

#Bedtools intersect all but Alaska
bedtools intersect -u -a Nunavut_W_C.bedGraph -b Svalbard_W_C.bedGraph| bedtools intersect -u -a stdin -b Sweden_W_C.bedGraph > ALEX_SVAL_SWED_Sites_intersect_DMRs.bedgraph

