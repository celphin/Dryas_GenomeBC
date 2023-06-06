#############################################################
#Trying to run Interproscan on DMRS:
#Note: Installed version does not work well with compute canada system + available dependencies, use module
#Notes on installation removed (have copy locally), as they do not work on the compute canada system
#Using: blast_ref_intersect_SE_W_C_SE_L_H.out , for testing -> later update to all necessary files
#############################################################
cd scratch
mkdir interproscan
cd interproscan
#############################################################
#cedar1
tmux new-session -s Interproscan
tmux attach-session -t Interproscan

#############################################################
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

#Grab 2nd column
cut -f2 blast_ref_intersect_Wild_W_C_Mat_Sen.out > blast_out.txt
#Map list of blasted genes onto protien
seqtk/seqtk subseq Dryas_octopetala_H1.protein.fa blast_out.txt > intersect_map.fasta
#Replace any gaps with X in file 
awk '{ gsub(/\*/, "X"); print }' intersect_map.fasta > intersect_map_clean.fasta

##############################################################
#Run interproscan:
#To do: make into script to run on all files in a directory

#Modules to load - don't use most recent version
module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0

interproscan.sh -i intersect_map_clean.fasta -f tsv -o interproscan_intersect_Wild_W_C_Mat_Sen.out
##############################################################




