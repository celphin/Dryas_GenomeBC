###############################################################################
#Find Parent DMRS with Metilene: 
    #Warming vs Control (Non site specific)
#Based on Notes2_Seedling_MetileneDMRS_Sept2022.txt
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
###################################################################################
#High plants Warming Vs Control
cd ~/scratch/
mkdir Wild_Metilene
cd Wild_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .
#In cedar5
tmux new-session -s Wild_Warming_DMRS
tmux attach-session -t Wild_Warming_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account rpp-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0


#Wild, Warming control specific adjustments:
#h1,h2: 
    #h1="W", h2="C"
#Input directories:
    #input_dir=Wild_${h1}_${h2}_input_files 
    #in_metilene="Wild_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Comment out rename for loops (already with proper prefix)
    #for bg in _*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in _*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

sh metilene_prep.sh 
cd ..

#----------------------------------------
#Run metilene:
salloc -c32 --time 2:50:00 --mem 120000m --account rpp-rieseber

#Wild,Warming control specific adjustments,
#h1,h2: 
    #h1="W", h2="C"
#in_metilene="Wild_metilene_"$h1"_"$h2".input"
#input_dir="Wild_${h1}_${h2}_input_files"
#outputname=Wild_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
#input_dir="/home/msandler/scratch/Wild_Metilene/Wild_${h1}_${h2}_input_files"
#outputname=Wild_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
#Low Plant Warming Vs Control












####################################################################################
#D. Octopetala Warming Vs Control





####################################################################################
#D. integrifolia Warming Vs. Control




















######################################################################################

# All Wild plants (not matching seedlings)
#Warming vs control DMRS


#

cd scratch
mkdir Wild_Metilene
cd Wild_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .
#In cedar5
tmux new-session -s Wild_Warming_DMRS
tmux attach-session -t Wild_Warming_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account rpp-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0


#Wild, Warming control specific adjustments:
#h1,h2: 
    #h1="W", h2="C"
#Input directories:
    #input_dir=Wild_${h1}_${h2}_input_files 
    #in_metilene="Wild_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Comment out rename for loops (already with proper prefix)
    #for bg in _*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in _*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

sh metilene_prep.sh 
cd ..

#----------------------------------------
#Run metilene:
salloc -c32 --time 2:50:00 --mem 120000m --account rpp-rieseber

#Wild,Warming control specific adjustments,
#h1,h2: 
    #h1="W", h2="C"
#in_metilene="Wild_metilene_"$h1"_"$h2".input"
#input_dir="Wild_${h1}_${h2}_input_files"
#outputname=Wild_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
#input_dir="/home/msandler/scratch/Wild_Metilene/Wild_${h1}_${h2}_input_files"
#outputname=Wild_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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



############################################################################################
#Wild:  Low Arctic(Alaska + Sweden) vs High Arctic (Svalbard + Alex sites)
#Low vs high DMRS 

cd scratch
mkdir Wild_Lat_Metilene
cd Wild_Lat_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .
#In cedar5
tmux new-session -s Wild_Lat_DMRS
tmux attach-session -t Wild_Lat_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account rpp-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0

#Note: see if there is a one of {a,b} option

#Wild, Low/High Latitude adjustments
#h1,h2: 
    #h1="L", h2="H"
#Input directories:
    #input_dir=Wild_Lat_${h1}_${h2}_input_files 
    #in_metilene="Wild_Lat_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Include following rename for loops 
    #for bg in *LAT*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in *ALAS*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in *CASS*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *FERT*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *DRY*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *MEAD*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *WILL*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *SVAL*; do mv "$bg" "${h2}_${bg}"; done


sh metilene_prep.sh 
cd ..

#------------------------------------------------------------------

#Wild, Latitude control specific adjustments,
#h1,h2: 
    #h1="L", h2="H"
#in_metilene="Wild_Lat_metilene_"$h1"_"$h2".input"
#input_dir=Wild_Lat_${h1}_${h2}_input_files 
#outputname=Wild_Lat_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
#h1="L", h2="H"
#input_dir="/home/msandler/scratch/Wild_Lat_Metilene/Wild_Lat_${h1}_${h2}_input_files"
#outputname=Wild_Lat_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
################################################################################################
#Wild: Species comparison: D. integrifolia (at Alex) vs D. octapetela (at Svalbard, Sweden, Alaska)

cd scratch
mkdir Wild_Species_Metilene
cd Wild_Species_Metilene
cp ~/projects/def-rieseber/Dryas_shared_data/MS_scripts/metilene*.sh .
#In cedar1
tmux new-session -s Wild_Spec_DMRS
tmux attach-session -t Wild_Spec_DMRS
salloc -c32 --time 4:50:00 --mem 120000m --account def-rieseber

#----------------------------------------
#Metilene input prep:
module load StdEnv/2020
module load bedtools/2.30.0

#Note: see if there is a one of {a,b} option

#Wild, Low/High Latitude adjustments
#h1,h2: 
    #h1="DO", h2="DI"
#Input directories:
    #input_dir=Wild_Species_${h1}_${h2}_input_files 
    #in_metilene="Wild_Species_metilene_"$h1"_"$h2".input"
    #methylseq_output_dir="/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Wild_metilene_input_bedGraphs"
#Include following rename for loops 
    #for bg in *LAT*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in *ALAS*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in *SVAL*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in *CASS*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *FERT*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *DRY*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *MEAD*; do mv "$bg" "${h2}_${bg}"; done
    #for bg in *WILL*; do mv "$bg" "${h2}_${bg}"; done
 


sh metilene_prep.sh 
cd ..

#------------------------------------------------------------------

#Wild, Species specific adjustments,
#h1,h2: 
    #h1="DO", h2="DI"
#in_metilene="Wild_Species_metilene_"$h1"_"$h2".input"
#input_dir=Wild_Species_${h1}_${h2}_input_files 
#output_name=Wild_Species_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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
#h1="DO", h2="DI"
#input_dir="/home/msandler/scratch/Wild_Species_Metilene/Wild_Species_${h1}_${h2}_input_files"
#outputname=Wild_Species_"$h1"_"$h2"_"${maxdist}"_"${mincpgs}"_"${mindiff}"
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