#############################################################
#Running Interproscan on DMRS blasted to reference
#Note: Installed version does not work well with compute canada system + available dependencies, use module
#Notes on installation removed (have copy locally), as they do not work on the compute canada system
#Running interproscan on:
    #all files in blast_ref directory, using scripts 
    #on a single file (intersection of Wild W C and Phenology shown)
#Requires:
    #prep_interproscan.sh : takes directory of blast outputs to reference, and outputs a folder ready to be run through interproscan
    #run_interproscan.sh: takes directory of correctly formatted fasta files, and runs interproscan with -goterms on all of them
    #Directory with blast .out files, to reference genome
#############################################################
cd scratch
mkdir interproscan
cd interproscan
#############################################################
#cedar1
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
#input_dir=home/msandler/projects/def-rieseber/Dryas_shared_data/MS_blast_output/Blast_ref_output/
#output_dir=interproscan_input
sh prep_interproscan.sh
#----------------------------------------------------
#Run interproscan:
module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0

#input_dir=home/msandler/scratch/interproscan/interproscan_input
#output_dir=home/msandler/scratch/interproscan/interproscan_output
sh run_interproscan.sh 



############################################################
#Running interproscan for single file
#Prepare files: 
#TO DO: make into script for whole directory
#TO DO: pipe file prep


#Copy files:
cp ~/projects/def-rieseber/Dryas_shared_data/MS_blast_output/Blast_ref_output/blast_ref_intersect_Wild_W_C_Mat_Sen.out .
cp ~/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.protein.fa

#Get seqtk
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
cd ..

#load modules
module load StdEnv/2020
module load gcc/9.3.0
module laod r-bundle-bioconductor/3.16 

#Include line to remove differing chromosome #'s:
awk 'index($2, substr($1, 1, index($1, ":")-1))' blast_ref_intersect_Wild_W_C_Mat_Sen.out > no_mismatch_intersect_Wild_W_C_Mat_Sen.out
#Grab 2nd column
cut -f2 no_mismatch_intersect_Wild_W_C_Mat_Sen.out > blast_out.txt
#Map list of blasted genes onto protien
seqtk/seqtk subseq Dryas_octopetala_H1.protein.fa blast_out.txt > intersect_map.fasta
#Replace any gaps with X in file 
awk '{ gsub(/\*/, "X"); print }' intersect_map.fasta > intersect_map_clean.fasta

#---------------------------------------------------------------------------------
#Run interproscan:
#To do: make into script to run on all files in a directory

#Modules to load - don't use most recent version
module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0

interproscan.sh -i intersect_map_clean.fasta -f tsv -o interproscan_intersect_Wild_W_C_Mat_Sen.out
##############################################################

