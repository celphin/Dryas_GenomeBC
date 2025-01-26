##################################
#Sept 2023
#MS_Blast_Notes_Sept2023.bash:
    #running reference bash on all DMR types
###########################################
#blast_prep_extended.sh
    #Maps bedGraph file to fasta, with DMRs extended
    #Param: string, .bedGraph file to be converted 
###########################################
#In cedar5
tmux new-session -s Blast
tmux attach-session -t Blast

cd scratch
mkdir blast
cd blast
##########################################
#Add copy files line
mkdir blast_bedgraphs
cd blast_bedgraphs
cp ~/projects/def-rieseber/Dryas_shared_data/MS_Bedgraph_Intersections/BedGraphs_From_Metilene/* . 
cd ..

##########################################
#preparing files for blasting:
module load bedtools/2.30.0
cp ~/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.supercontigs.fa .


sh blast_prep_extended.sh Alaska_W_C
sh blast_prep_extended.sh Nunavut_W_C
sh blast_prep_extended.sh Svalbard_W_C
sh blast_prep_extended.sh Sweden_W_C
sh blast_prep_extended.sh SE_L_H 
sh blast_prep_extended.sh Wild_Lat_L_H
sh blast_prep_extended.sh Mat_Sen
sh blast_prep_extended.sh Parent_W_C
sh blast_prep_extended.sh SE_W_C
sh blast_prep_extended.sh Wild_W_C


mkdir fasta_dir
mv *.fasta fasta_dir/
#################################################
#Running blast to reference
#In cedar5:
tmux new-session -s Blast
tmux attach-session -t Blast

cd ~/scratch/blast
#module spider blast+ for most recent version: 2.13.0 right now 
module load StdEnv/2020
module load gcc/9.3.0
module load blast+/2.13.0 
fasta_dir="fasta_dir"

cp ~/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/*.transcript.fa .
#If any of above files not in blast directory:

blastname="Alaska_W_C"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="Nunavut_W_C"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="Svalbard_W_C"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="Sweden_W_C"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="SE_L_H"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="Wild_Lat_L_H"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="Mat_Sen"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="Parent_W_C"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="SE_W_C"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
blastname="Wild_W_C"
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6



cp *.out ~/projects/def-rieseber/Dryas_shared_data/MS_blast_output
####################################################################################################################################
#Leaves only blasts matching to the same part of the chromosome eg "Do1_01_a0004"
cd ~/projects/def-rieseber/Dryas_shared_data/MS_blast_output/

for file in *
do awk 'index($2, substr($1, 1, index($1, ":")-1))' "$file" > cleaned_"$file"
done

mkdir raw_blast_ref
mkdir cleaned_blast_ref
mv cleaned_* cleaned_blast_ref
mv blast_ref_* raw_blast_ref
####################################################################################################################################
