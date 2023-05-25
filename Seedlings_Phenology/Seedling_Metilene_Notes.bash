########################################
#Find seedling DMRS with Metilenes: 
    #Warming vs Control
    #High vs low arctic conditions
#Based on Notes2_Seedling_MetileneDMRS_Sept2022.txt
########################################
##metile_prep.sh:
    #adjusts all the files names based on categories,
    #prepares them in the correct format for metilene
########################################
#Copy Methylseq output to scratch directory:

cp ~/projects/def-rieseber/Dryas_shared_data/CE_Seedling_metilene_input_bedGraphs/ ~/scratch/Seedling_Metilene


########################################
#Warming vs control DMRS

#Metilene input prep:
cd ~/scratch/Seedling_Metilene
tmux new-session -s Seedling_Warming_DMRS
tmux attach-session -t Seedling_Warming_DMRS

module load StdEnv/2020
module load bedtools/2.30.0


#Warming control specific adjustments:
#h1,h2: 
    #h1="W", h2="C"
#Input directorys:
    #input_dir=SE_${h1}_${h2}_input_files 
    #in_metilene="SE_metilene_"$h1"_"$h2".input"
#Adjust rename for loop to be:
    #for bg in SE_*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in SE_*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

sh metilene_prep.sh 

#Run metilene:

#

######################################
#High vs low DMRS
#Metilene input prep:
cd ~/scratch/Seedling_Metilene
tmux new-session -s Seedling_Warming_DMRS
tmux attach-session -t Seedling_Warming_DMRS

module load StdEnv/2020
module load bedtools/2.30.0

#metile_prep.sh:
    #adjusts all the files names based on categories,
    #prepares them in the correct format for metilene
#Warming control specific adjustments:
#h1,h2: 
    #h1="L", h2="C"
#Input directorys:
    #input_dir=SE_${h1}_${h2}_input_files 
    #in_metilene="SE_metilene_"$h1"_"$h2".input"
#Adjust rename for loop to be:
    #for bg in SE_*_${h1}*; do mv "$bg" "${h1}_${bg}"; done
    #for bg in SE_*_${h2}*; do mv "$bg" "${h2}_${bg}"; done

sh metilene_prep.sh 

#Run metilene:

