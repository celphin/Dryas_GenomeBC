########################################
#Find Phenology DMRS with Metilene: 
    #Mature flower vs Senesence (Non site specific)
    #Notes: haven't run this yet
#Based on Notes4_Phenology_MetileneDMRS_Sept2022.txt
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

cd scratch
mkdir Phenology_Metilene
cd Phenology_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .

#In cedar5
tmux new-session -s Phenology_DMRS
tmux attach-session -t Phenology_DMRS
salloc -c32 --time 2:50:00 --mem 120000m --account rpp-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0


#Phenology specific adjustments:
#h1,h2: 
    #h1="Mat", h2="Sen"
#Input directories:
    #input_dir=${h1}_${h2}_input_files 
    #in_metilene="metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/Dryas/CE_Phenology_metilene_input_bedGraphs"
#change copy process because Senesence bedgraphs don't have deduplicated:
    #Under: cp -r ${methylseq_output_dir}/*.deduplicated.bedGraph ./${input_dir}
    #Add:   cp -r ${methylseq_output_dir}/Sen_*.bedGraph ./${input_dir}
#Comment out rename for loops shown below (already with proper prefix)
    ##for bg in _*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    ##for bg in _*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

sh metilene_prep.sh 
cd ..

#----------------------------------------
#Run metilene:
salloc -c32 --time 2:50:00 --mem 120000m --account rpp-rieseber

#Phenology specific adjustments,
#h1,h2: 
    #h1="Mat", h2="Sen"
#in_metilene="metilene_"$h1"_"$h2".input"
#input_dir=${h1}_${h2}_input_files 
#outputname="$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
#h1="Mat", h2="Sen"
#input_dir="/home/msandler/scratch/Phenology_Metilene/${h1}_${h2}_input_files"
#outputname="$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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

