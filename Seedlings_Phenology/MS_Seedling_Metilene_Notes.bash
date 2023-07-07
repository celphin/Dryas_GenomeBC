########################################
#Find seedling DMRS with Metilene: 
    #Warming vs Control
    #High vs low arctic conditions
#Based on Notes2_Seedling_MetileneDMRS_Sept2022.txt
########################################
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
########################################
#Warming vs control DMRS

cd scratch
mkdir Seedling_Metilene
cd Seedling_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .
#In cedar5
tmux new-session -s Seedling_Warming_DMRS
tmux attach-session -t Seedling_Warming_DMRS

#----------------------------------------
#Metilene input prep:

#Seedling, Warming control specific adjustments:
#h1,h2: 
    #h1="W", h2="C"
#Input directories:
    #input_dir=SE_${h1}_${h2}_input_files 
    #in_metilene="SE_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Seedling_metilene_input_bedGraphs"
#Adjust rename for loop to be:
    #for bg in SE_*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in SE_*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

module load StdEnv/2020
module load bedtools/2.30.0
sh metilene_prep.sh 

#----------------------------------------
#Run metilene:
salloc -c32 --time 2:50:00 --mem 120000m --account rpp-rieseber

#Seedling,Warming control specific adjustments,
#h1,h2: 
    #h1="W", h2="C"
#in_metilene="SE_metilene_"$h1"_"$h2".input"
#input_dir=SE_${h1}_${h2}_input_files 
#outputname=SE_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
#parameters: maxdist, mincpgs, mindiff

#metilene:
# params: maxdist, mincpgs, mindiff
module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7

#----------------------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value
#h1="W", h2="C"
#input_dir="/home/msandler/scratch/Seedling_Metilene/SE_${h1}_${h2}_input_files"
#outputname=SE_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
#Seedling, High vs low DMRS 
#Metilene input prep:
cd ~/scratch/Seedling_Metilene
#In cedar5
tmux new-session -s Seedling_HighLow_DMRS
tmux attach-session -t Seedling_HighLow_DMRS

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0


#Seedling, Low/High Arctic - specific adjustments:
#h1,h2: 
    #h1="L", h2="H"
#Input directories (same as for Seedling Warming):
    #input_dir=SE_${h1}_${h2}_input_files 
    #in_metilene="SE_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Seedling_metilene_input_bedGraphs"
#Adjust rename for loop to be:
    #for bg in SE_*_*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in SE_*_*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

sh metilene_prep.sh 
cd ..
#----------------------------------------
#Run metilene:
salloc -c32 --time 2:50:00 --mem 120000m --account rpp-rieseber

#Seedling, Low/High arctic control specific adjustments
#h1,h2: 
    #h1="L", h2="H"
#in_metilene="SE_metilene_"$h1"_"$h2".input"
#input_dir=SE_${h1}_${h2}_input_files 
#outputname=SE_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
#parameters: maxdist, mincpgs, mindiff

#metilene:
# params: maxdist, mincpgs, mindiff
module load StdEnv/2020
module load bedtools/2.30.0 
sh metilene_run.sh 70 5 4
sh metilene_run.sh 150 5 4 
sh metilene_run.sh 70 5 0.7

#-------------------------------------
#Filter based on qval
#params: maxdist, mincpgs, mindiff, minmeandif, q-value
#h1="L", h2="H", in file
#input_dir="/home/msandler/scratch/Seedling_Metilene/SE_${h1}_${h2}_input_files"
#outputname=SE_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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