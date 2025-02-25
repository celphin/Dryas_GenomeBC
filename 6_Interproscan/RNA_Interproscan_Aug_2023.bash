#MS_RNA_Interproscan.bash
    #Getting list of genes from RNASeq results
    #Needed scripts: 
        #run_interproscan_full_dryas.sh
    #Needed input files:
      #RNA_DER_Aug2023_W_C_Total.txt 
      #CE_Dryas_reference_genome/Dryas_octopetala_H1.protein.fa
###################################################
cd ~/scratch/interproscan
mkdir RNA_Interproscan

cp  /home/msandler/projects/def-rieseber/Dryas_shared_data/MS_RNAseq_DERs/RNA_DER_Aug2023_W_C_Total.txt . 
cp  /home/msandler/projects/def-rieseber/Dryas_shared_data/CE_Dryas_reference_genome/Dryas_octopetala_H1.protein.fa . 

#-------------------------------------------------
#Install seqtk
git clone https://github.com/lh3/seqtk.git;
cd seqtk; make
cd ..


#-------------------------------------------------
#Overlap to fasta
cut -f1 RNA_DER_Aug2023_W_C_Total.txt > col1.txt
seqtk/seqtk subseq Dryas_octopetala_H1.protein.fa col1.txt |awk '{ gsub(/\*/, "X"); print }' > RNA_DER_W_C_Total_protein.fa
rm col1.txt

#-------------------------------------------------
#Run Interproscan

sbatch run_interproscan_rna.sh 
#######################################################################################################
