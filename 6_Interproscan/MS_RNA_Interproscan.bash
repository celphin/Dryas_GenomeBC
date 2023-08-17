#MS_RNA_Interproscan.bash
    #Getting list of real genes from RNASeq results

###################################################
cd ~/scratch/interproscan
mkdir RNAseq_mapping

cp  /home/msandler/projects/def-rieseber/Dryas_shared_data/CE_RNAseq_DERs/RNA_DER_May2023_W_C_Total.txt . 
cp  /home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.protein.fa . 

#-------------------------------------------------
#Overlap to fasta
cut -f1 RNA_DER_May2023_W_C_Total.txt > col1.txt
seqtk/seqtk subseq Dryas_octopetala_H1.protein.fa col1.txt |awk '{ gsub(/\*/, "X"); print }' > rna_intersect_map.fasta
rm col1.txt
mv rna_intersect_map.fasta RNAseq_mapping
cd RNAseq_mapping

#-------------------------------------------------
#Run Interproscan
module load nixpkgs/16.09 
module load intel/2018.3
module load kentutils/20180716
mkdir seq_dir
faSplit sequence rna_intersect_map.fasta 20 seq_dir/seq
#Check to see how many made (sometimes not perfect split)

mkdir run_scripts

for i in `seq 0 18`; do
cat << EOF > run_scripts/run_interproscan_${i}.sh
#!/bin/bash
#SBATCH --account=def-rieseber
#SBATCH --time=0-6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=160000m

module load nixpkgs/16.09
module load intel/2016.4
module load interproscan/5.23-62.0
  
 
if [ $i -lt 10 ]; then
  srun interproscan.sh -i "seq_dir/seq0$i.fa" -f tsv -o interproscan_rna_$i.tsv  --goterms
else 
  srun interproscan.sh -i "seq_dir/seq$i.fa" -f tsv -o interproscan_rna_$i.tsv  --goterms
fi

EOF

sbatch run_scripts/run_interproscan_${i}.sh
done

cat interproscan_rna_{0..9}.tsv > interproscan_rna.tsv

#######################################################################################################
#sort interproscan_blast_ref_intersect_Wild_W_C_Mat_Sen.tsv | 
#All Id's with go terms
grep 'GO' interproscan_rna.tsv | awk  '!a[$1$NF]++ {print}'|awk '!a[$1]++ {print $1}' > list_true_ids.txt
grep -f list_true_ids.txt interproscan_blast_ref_intersect_Wild_W_C_Mat_Sen.tsv > subset_intersect_Wild_W_C_Mat_Sen.tsv #Overlap Intersect W_C + Phenology
grep -f list_true_ids.txt interproscan_blast_ref_total_subtract_W_C_Mat_Sen.tsv > subset_subtract_W_C_Mat_Sen.tsv #Overlap Warming-Phenology
#######################################################################
#Getting all RNAseq interlaps with 
#List all RNAseq id's
awk  '!a[$1]++ {print $1}' interproscan_rna.tsv > rnaseq_gene_ids.txt

#For Warming-Phenology: all ID's with go terms, all Id's without go terms
grep 'GO' interproscan_blast_ref_total_subtract_W_C_Mat_Sen.tsv| awk  '!a[$1]++ {print $1}'  > go_subtract_W_C_Mat_Sen_id.tsv
grep -v -f go_subtract_W_C_Mat_Sen_id.tsv interproscan_blast_ref_total_subtract_W_C_Mat_Sen.tsv| awk  '!a[$1]++ {print}' > no_go_subtract_W_C_Mat_Sen_id.tsv

#Intersection with rnaseq gene id's
grep -f rnaseq_gene_ids.txt go_subtract_W_C_Mat_Sen_id.tsv > overlapped_go_subtract_W_C_Mat_Sen_id.tsv
grep -f rnaseq_gene_ids.txt no_go_subtract_W_C_Mat_Sen_id.tsv > overlapped_no_goterms_overlapped_warming.tsv

#Get these overlapped genes back to interproscan - > put that output in shared_dir

grep -f overlapped_go_subtract_W_C_Mat_Sen_id.tsv interproscan_blast_ref_total_subtract_W_C_Mat_Sen.tsv |grep 'GO'| awk  '!a[$1$NF]++ {print}' > interproscan_goterms_overlap_rna_subtract_W_C_Mat_Sen.tsv