########################################
#Find Parent DMRS with Metilene: 
    #Warming vs Control (Non site specific)
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
mkdir Parent_Metilene
cd Parent_Metilene
#In cedar5
tmux new-session -s Parent_Warming_DMRS
tmux attach-session -t Parent_Warming_DMRS
salloc -c32 --time 2:50:00 --mem 120000m --account rpp-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0


#Seedling, Warming control specific adjustments:
#h1,h2: 
    #h1="W", h2="C"
#Input directories:
    #input_dir=P_${h1}_${h2}_input_files 
    #in_metilene="P_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Comment out rename for loops (already with proper prefix)
    #for bg in _*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in _*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

sh metilene_prep.sh 
cd ..

#----------------------------------------
#Run metilene:
salloc -c32 --time 2:50:00 --mem 120000m --account rpp-rieseber

#Seedling,Warming control specific adjustments,
#h1,h2: 
    #h1="W", h2="C"
#in_metilene="P_metilene_"$h1"_"$h2".input"
#input_dir=P_${h1}_${h2}_input_files 
#outputname=P"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
#input_dir="/home/msandler/scratch/Seedling_Metilene/P_${h1}_${h2}_input_files"
#outputname=P_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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

####################################################################
#Site specific DMRS:

