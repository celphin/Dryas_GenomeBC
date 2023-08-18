##############################################
#prep_interproscan.sh: 
#Input: directory of files blasted to reference genome

#Output: directory of fasta files ready to load into Interproscan
#Required modules: 
    #module load StdEnv/2020
    #module load gcc/9.3.0
    #module laod r-bundle-bioconductor/3.16 
##############################################
ref_protein=/home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.protein.fa
input_dir=/home/msandler/projects/def-rieseber/Dryas_shared_data/MS_blast_output/Blast_ref_output/cleaned_blast_ref
output_dir=interproscan_input

mkdir $output_dir

#For every file in a directory
for f in $input_dir/*; do
    #Gets file name
    no_path="${f##*/}"
    file_name="${no_path%.out}"
    #Substrings the gene name 
    awk 'index($2, substr($1, 1, index($1, ":")-1))' $f | cut -f2 > blast_out.txt
    #Interlaps with fasta and sends to new file
    seqtk/seqtk subseq $ref_protein blast_out.txt| awk '{ gsub(/\*/, "X"); print }' > $output_dir/$file_name.fasta
done

rm blast_out.txt

#awk 'index($2, substr($1, 1, index($1, ":")-1))' blast_ref_intersect_Wild_W_C_Mat_Sen.out | cut -f2 |\
#seqtk/seqtk subseq Dryas_octopetala_H1.protein.fa | awk '{ gsub(/\*/, "X"); print }' 

