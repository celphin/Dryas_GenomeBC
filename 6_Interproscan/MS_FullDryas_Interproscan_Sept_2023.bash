#MS_FullDryas_Interproscan_Sept_2023.bash
    #Runs interproscan on full Dryas protein file
    #Needed scripts: run_interproscan_full_dryas.sh
    #Needed input files:
        #Dryas_octopetala_H1.protein.fa
###################################################
cd ~/scratch/interproscan
mkdir Dryas_Interproscan

cp  /home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.protein.fa . 
awk '{ gsub(/\*/, "X"); print }' Dryas_octopetala_H1.protein.fa > Dryas_protein.fa

#Checking size
grep ">" Dryas_octopetala_H1.protein.fa | wc -l 
#41181 - okay as single file

sbatch run_interproscan_full_dryas.sh
