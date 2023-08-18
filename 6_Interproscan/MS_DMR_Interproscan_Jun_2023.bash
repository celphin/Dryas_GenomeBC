#############################################################
#Old Version: See MS_DMRs_Interproscan_Aug_2023.bash, for more recent version of this 
#Running Interproscan on DMRS blasted to reference
#Note: Installed version does not work well with compute canada system + available dependencies, use module
#Running interproscan on:
    #Prepares files to run interproscan on, using interproscan_prep.sh
    #all (small) files in blast_ref directory, using sbatch
    #Splitting up large files into multiple, to speed up run - interproscan gets stuck otherwise
#filtering output:
    #TODO: Removing duplicates (not done)  
    #TODO: make files by go terms for each file, make one group
#Notes on running interproscan on single file included at the bottom of notes
#Requires: 
    #prep_interproscan.sh: 
        #reference protien file
        #blast to reference file
        #directory for input
#Notes on fasplit:
    #TODO: Figure out issue with intersecting faSplit outputs
    #Some genes very clearly repeat themselves
    #faSplit was reccomended (and all according dependencies) by compute canada for splitting fasta files for blasting
#############################################################
cd scratch
mkdir Interproscan
cd Interproscan
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

ls interproscan_input > input_files.txt

#for f in `seq 1 50`; do
#cat << EOF > tet5_${i}.sh
ls interproscan_input > input_files.txt
mkdir run_scripts

while read f; do

no_path="${f##*/}"
file_name="${no_path%.fasta}"
cat << EOF > run_scripts/interproscan_run_"$file_name".sh
#!/bin/sh
#SBATCH --account=rpp-rieseber
#SBATCH --time=3-0:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=160000m

module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0
    
srun interproscan.sh -i interproscan_input/$f -f tsv -o interproscan_output/interproscan_$file_name.tsv  --goterms

EOF

sbatch run_scripts/interproscan_run_"$file_name".sh

done < "input_files.txt"

cp *.tsv /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_interproscan_output
#-----------------------------------------------------
#Files that are too big for full fasta processing:
    #blast_ref_intersect_SE_W_C_SE_L_H.fasta 
    #blast_ref_Mat_Sen.fasta     
    #blast_ref_total_subtract_SE_W_C_P_W_C.fasta
    #blast_ref_total_subtract_SE_W_C_SE_L_H.fasta

#Processeed below by spltting fasta's

#---------------------------------------------------
#blast_ref_intersect_SE_W_C_SE_L_H.fasta 

cd ~/scratch/interproscan
mkdir interproscan_intersect_SE_W_C_SE_L_H
cd interproscan_intersect_SE_W_C_SE_L_H
cp ~/scratch/interproscan/interproscan_input/blast_ref_intersect_SE_W_C_SE_L_H.fasta . 

module load nixpkgs/16.09 
module load intel/2018.3
module load kentutils/20180716
faSplit sequence blast_ref_intersect_SE_W_C_SE_L_H.fasta 10 seq

mkdir subseq
mkdir run_scripts
cp *.fa subseq

for i in `seq 0 9`; do
cat << EOF > run_scripts/run_interproscan_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=160000m

module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0
   
srun interproscan.sh -i "subseq/seq0$i.fa" -f tsv -o interproscan_intersect_SE_W_C_SE_L_H_$i.tsv  --goterms

EOF

sbatch run_scripts/run_interproscan_${i}.sh
done

#wait until all done:
cat interproscan_intersect_SE_W_C_SE_L_H_{0..9}.tsv > interproscan_blast_ref_intersect_SE_W_C_SE_L_H.tsv
cp interproscan_blast_ref_intersect_SE_W_C_SE_L_H.tsv ../interproscan_output
cd ..

#---------------------------------------------------------
#blast_ref_Mat_Sen.fasta  

cd ~/scratch/interproscan
mkdir interproscan_Mat_Sen
cd interproscan_Mat_Sen
cp ~/scratch/interproscan/interproscan_input/blast_ref_Mat_Sen.fasta .

module load nixpkgs/16.09 
module load intel/2018.3
module load kentutils/20180716
faSplit sequence blast_ref_Mat_Sen.fasta 300 seq

mkdir subseq
mkdir run_scripts
cp *.fa subseq

for i in `seq 0 277`; do
cat << EOF > run_scripts/run_interproscan_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=160000m

module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0
   
if [ $i -lt 10 ]; then
  srun interproscan.sh -i "subseq/seq00$i.fa" -f tsv -o interproscan_Mat_Sen_$i.tsv  --goterms
elif [ $i -lt 100 ]; then
  srun interproscan.sh -i "subseq/seq0$i.fa" -f tsv -o interproscan_Mat_Sen_$i.tsv  --goterms
else
  srun interproscan.sh -i "subseq/seq$i.fa" -f tsv -o interproscan_Mat_Sen_$i.tsv  --goterms
fi


EOF

sbatch run_scripts/run_interproscan_${i}.sh
done

cat interproscan_Mat_Sen_{0..277}.tsv > interproscan_blast_ref_Mat_Sen.tsv
cp interproscan_blast_ref_intersect_SE_W_C_SE_L_H.tsv ../interproscan_output
cd ..

#----------------------------------------------------------------------------------
#blast_ref_total_subtract_SE_W_C_P_W_C.fasta #in progress
#wc -l blast_ref_total_subtract_SE_W_C_P_W_C.fasta
#1870 -> split into 187?

cd ~/scratch/interproscan
mkdir interproscan_total_subtract_SE_W_C_P_W_C
cd interproscan_total_subtract_SE_W_C_P_W_C
cp ~/scratch/interproscan/interproscan_input/blast_ref_total_subtract_SE_W_C_P_W_C.fasta .

module load nixpkgs/16.09 
module load intel/2018.3
module load kentutils/20180716
faSplit sequence blast_ref_total_subtract_SE_W_C_P_W_C.fasta 100 seq
#Check to see how many made (sometimes not perfect split)

mkdir subseq
mkdir run_scripts
cp *.fa subseq

for i in `seq 0 93`; do
cat << EOF > run_scripts/run_interproscan_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=160000m

module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0
  
 
if [ $i -lt 10 ]; then
  srun interproscan.sh -i "subseq/seq00$i.fa" -f tsv -o interproscan_intersect_total_subtract_SE_W_C_P_W_C_$i.tsv.tsv  --goterms
else 
  srun interproscan.sh -i "subseq/seq0$i.fa" -f tsv -o interproscan_intersect_total_subtract_SE_W_C_P_W_C_$i.tsv  --goterms
fi

EOF

sbatch run_scripts/run_interproscan_${i}.sh
done


cat {0..93}.tsv > interproscan_blast_ref_total_subtract_SE_W_C_P_W_C.tsv
cp interproscan_blast_ref_total_subtract_SE_W_C_P_W_C.tsv ../interproscan_output
cd ..

#-------------------------------------------------------------------------------
#blast_ref_total_subtract_SE_W_C_SE_L_H.fasta #not done yet
#100 parts 
cd ~/scratch/interproscan
mkdir interproscan_total_subtract_SE_W_C_SE_L_H
cd interproscan_total_subtract_SE_W_C_SE_L_H
cp ~/scratch/interproscan/interproscan_input/blast_ref_total_subtract_SE_W_C_SE_L_H.fasta .

module load nixpkgs/16.09 
module load intel/2018.3
module load kentutils/20180716
faSplit sequence blast_ref_total_subtract_SE_W_C_SE_L_H.fasta 100 seq
#Check to see how many made

mkdir subseq
mkdir run_scripts
cp *.fa subseq

for i in `seq 0 93`; do
cat << EOF > run_scripts/run_interproscan_${i}.sh
#!/bin/bash
#SBATCH --account=rpp-rieseber
#SBATCH --time=0-7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=160000m

module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0
  
 
if [ $i -lt 10 ]; then
  srun interproscan.sh -i "subseq/seq00$i.fa" -f tsv -o interproscan_intersect_total_subtract_SE_W_C_SE_L_H_$i.tsv  --goterms
else 
  srun interproscan.sh -i "subseq/seq0$i.fa" -f tsv -o interproscan_intersect_total_subtract_SE_W_C_SE_L_H_$i.tsv  --goterms
fi

EOF

sbatch run_scripts/run_interproscan_${i}.sh
done



cp interproscan_blast_ref_total_subtract_SE_W_C_SE_L_H.tsv ../interproscan_output
cd ..
##########################################################################
#Removing duplicates:


#group interproscan results by GO terms
#grep "GO:0005515" go_terms_interproscan_intersect_Wild_W_C_Mat_Sen.out | cut -f 1,6 --output-delimiter=_ |sort|uniq
#grep "GO:" go_terms_interproscan_intersect_Wild_W_C_Mat_Sen.out |cut -f14 |sort | uniq -c | sort -bgr


#Experimenting:
#wc: 403
sort interproscan_blast_ref_intersect_Wild_W_C_Mat_Sen.tsv | awk '!a[$1]++ {print}' > no_first_field_dubs.tsv


#------------------------------------------------------------------
#Grouping by Go terms:
cd interproscan_output


ls interproscan_input > output_files.txt
while read f; do 
    #grep "GO:" go_terms_interproscan_intersect_Wild_W_C_Mat_Sen.out |cut -f14 |sort | uniq -c | sort -bgr > f
done < "output_files.txt"




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

interproscan.sh -i intersect_map_clean.fasta -f tsv -o interproscan_intersect_Wild_W_C_Mat_Sen.out -goterms
##############################################################
