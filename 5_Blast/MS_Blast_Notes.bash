##################################
#Blast_Notes.bash:
    #Make intersection and difference bedgraphs
    #Prepare bedgraphs for blast
    #Blast bedgraphs
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
##########################################
#preparing files for blasting:
module load bedtools/2.30.0
cp ~/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.supercontigs.fa .

sh blast_prep_extended.sh "Mat_Sen"
sh blast_prep_extended.sh "intersect_Wild_W_C_Mat_Sen"
sh blast_prep_extended.sh "total_subtract_W_C_Mat_Sen"
sh blast_prep_extended.sh "intersect_SE_W_C_SE_L_H"
sh blast_prep_extended.sh "total_subtract_SE_W_C_SE_L_H"
sh blast_prep_extended.sh "intersect_SE_W_C_P_W_C"
sh blast_prep_extended.sh "intersect_SE_W_C_Wild_W_C"
sh blast_prep_extended.sh "total_subtract_SE_W_C_P_W_C"
sh blast_prep_extended.sh "intersect_SE_L_H_Wild_L_H"

mkdir fasta_dir
mv *.fasta fasta_dir/
#################################################
#Running blast: ncbi Rosaceae, Arabidopsis, and to reference
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


blastname="Mat_Sen"
#output format6, against rosaceaea 
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against  (running in Blast)
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#against reference: 
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6


blastname="intersect_Wild_W_C_Mat_Sen"
#ncbi, against rosaceaea
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against Arabidopsis:
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#against reference: 
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6


blastname="total_subtract_W_C_Mat_Sen"
#ncbi, against rosaceae
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against Arabidopsis:
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#against reference: 
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6

blastname="intersect_SE_W_C_SE_L_H"
#ncbi, against rosaceaea
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against Arabidopsis:
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#against reference: 
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6

blastname="total_subtract_SE_W_C_SE_L_H"
#ncbi, against rosaceaea
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against Arabidopsis: 
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#against reference: 
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6

blastname="intersect_SE_W_C_P_W_C"
#ncbi against rosaceaea 
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against Arabidopsis: (Running in Blast2)
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#agaisnt reference
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6

blastname="intersect_SE_W_C_Wild_W_C"
#agaisnt reference
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
#ncbi against rosaceaea 
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against Arabidopsis: 
##TO DO: rerun
blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"


blastname="total_subtract_SE_W_C_P_W_C"
#agaisnt reference
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
#ncbi against rosaceaea 
#blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against Arabidopsis: 
##TO DO: rerun
#blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"


blastname="intersect_SE_L_H_Wild_L_H"
#agaisnt reference
blastn -query "${fasta_dir}/${blastname}.fasta" -out "blast_ref_${blastname}.out" -subject Dryas_octopetala_H1.transcript.fa -outfmt 6
#ncbi against rosaceaea 
#blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_rosaceae_${blastname}.out" -entrez_query "Rosaceae [Family]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"
#ncbi, against Arabidopsis: 
#blastn -db nt -query "${fasta_dir}/${blastname}.fasta" -out "blast_ncbi_arabidopsis_${blastname}.out" -entrez_query "Arabidopsis [Genus]" -remote -outfmt "6 qseqid sseqid pident stitle length mismatch gapopen qstart qend sstart send evalue bitscore"


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
