#############################################################
#Running Interproscan on DMRS blasted to reference
#August 2023 (most recent version)
#Requires:
    #prep_interproscan.sh:
        #reference protien file
        #blast to reference file
        #directory for input
    #run_interproscan_dir.sh - runs interproscan on a folder of files
#############################################################
cd scratch
mkdir DMR_Interproscan
cd DMR_Interproscan
#############################################################
#cedar5
tmux new-session -s Interproscan
tmux attach-session -t Interproscan

#############################################################
#Running Interproscan for whole directory
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
cd ..

#---------------------------------------------------
#Prep interproscan files
#load modules
module load StdEnv/2020
module load gcc/9.3.0
module laod r-bundle-bioconductor/3.16 

#ref_genome=home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.protein.fa
#input_dir=home/msandler/projects/def-rieseber/Dryas_shared_data/MS_blast_output/Blast_ref_output/cleaned_blast_ref
#output_dir=interproscan_input
sh prep_interproscan.sh
##############################################################################
#Run interproscan:

sbatch run_interproscan_dir.sh

#############################################################################
