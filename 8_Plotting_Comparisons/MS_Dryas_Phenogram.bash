#########################################################
#November 2023
#Steps:
    #Get 
    #Convert to Gene: start-end
    #Repeat with 
    #Get genome file for phenogram
#########################################################


wget https://ritchielab.org/files/RL_software/ruby_install.sh
wget https://ritchielab.org/files/RL_software/pheno_gram.rb
sh ruby_install.sh
ruby pheno_gram.rb
#######################################################
#Acquire Dryas_Genome.txt:

module load bioawk
bioawk -c fastx '{ print $name, length($seq) }' < Dryas_octopetala_H1.supercontigs.fa

#Hand Edit to be tsv file, with header:
#######################################################
grep -E "SE_W_C|Wild_W_C|Parent_W_C|POS|RNA" Dryas_Total_Origin_Phenogram_Chr.txt > Dryas_W_C_SE_Wild_P_RNA_Phenogram.txt   
grep -E "SE_L_H|Wild_Lat_L_H|POS" Dryas_Total_Origin_Phenogram_Chr.txt > Dryas_L_H_SE_Wild_Phenogram.txt   
grep -E "Sval|Alas|Swed|Nunavut|POS" Dryas_Total_Origin_Phenogram_Chr.txt > Dryas_Loc_W_C_Phenogram.txt   
grep -E "Mat_Sen|POS" Dryas_Total_Origin_Phenogram_Chr.txt > Dryas_Phenology_Phenogram.txt  
#######################################################
ruby pheno_gram.rb -i Dryas_W_C_SE_Wild_P_RNA_Phenogram.txt -g Dryas_Genome.txt -t "Dryas Warming Phenogram" -o Dryas_W_C_SE_Wild_P_RNA_Phenogram -f jpg
ruby pheno_gram.rb -i Dryas_Loc_W_C_Phenogram.txt -g Dryas_Genome.txt -t "Dryas Location Phenogram" -o Dryas_Loc_W_C_Phenogram -f jpg
ruby pheno_gram.rb -i Dryas_Phenology_Phenogram.txt -g Dryas_Genome.txt -t "Dryas Phenology Phenogram" -o Dryas_Phenology_Phenogram -f jpg
#Too many, figure out way to subset
#ruby pheno_gram.rb -i Dryas_L_H_SE_Wild_Phenogram.txt -g Dryas_Genome.txt -t "Dryas Latitude Phenogram" -o Dryas_L_H_SE_Wild_Phenogram -f jpg
#######################################################\
scp -v msandler@graham.computecanada.ca:/home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Dryas_Phenogram/*.jpg .

scp -v msandler@graham.computecanada.ca:/home/msandler/scratch/Dryas_W_C_SE_Wild_P_RNA_Phenogram.jpg .

########################################################

ruby pheno_gram.rb -i Dryas_Intersecting_W_C_Phenogram_Chr.txt -g Dryas_Genome.txt -t "Dryas Warming Intersections Phenogram" -o Dryas_Intersect_W_C_SE_Wild_P_RNA_Phenogram -f jpg

scp -v msandler@graham.computecanada.ca:/home/msandler/projects/def-rieseber/Dryas_shared_data/MS_Dryas_Phenogram/Dryas_Intersect_W_C_SE_Wild_P_RNA_Phenogram.jpg .

##################################
# April 4th Phenograms

#------------------------
# Install
# wget https://ritchielab.org/files/RL_software/ruby_install.sh
# wget https://ritchielab.org/files/RL_software/pheno_gram.rb
# sh ruby_install.sh
# ruby pheno_gram.rb
#
#-----------------------------
# Make input files
#
# Genome file format:
# ID	size	centromere
# 1       Chr_Length      start,end
#
# Input file format:
# CHR	POS	PHENOTYPE
# 1       63160   GO:0005515

#------------------------------
# GO terms

cd /home/celphin/projects/def-rieseber/Dryas_shared_data/MS_Dryas_Phenogram
cp *  ~/scratch/Dryas/MS_Dryas_Merged_Data/Phenograms/
cd  ~/scratch/Dryas/MS_Dryas_Merged_Data/Phenograms/

# test

ruby pheno_gram.rb -i Dryas_Intersecting_W_C_Phenogram_Chr.txt -g Dryas_Genome.txt -t "Dryas Warming Intersections Phenogram" -o Test -f jpg

#---------------------------------
# Make image of defense genes

# DMR Warming
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Wild_W_C"  | awk '{print $1 " " $11 " " $14}'  | sort | uniq  > Pheno_go_terms_Wild.txt
sed 's/a0000//g' Pheno_go_terms_Wild.txt | sed 's/Do1_0//g' | sed 's/_//' | sed 's/G0....//' | sed 's/ /\t/g'  > Phenogram_GOterms_Wild_w_c.txt
# add header and remove contigs
nano Phenogram_GOterms_Wild_w_c.txt

cp Phenogram_GOterms_Wild_w_c.txt Phenograms/
cd Phenograms/
ruby pheno_gram.rb -i Phenogram_GOterms_Wild_w_c.txt -g Dryas_Genome.txt -t "Dryas DMR Warming GO terms" -o DMR_Warming_GO -f jpg

# copy to local machine
scp -v celphin@cedar.computecanada.ca:/home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/Phenograms/DMR_Warming_GO.jpg .

#------------------------
# RNA
grep "GO:" Gene_RNA_Total_GO_Merged_table.tsv  | awk '{print $1 " " $7 " " $10}'  | sort | uniq  > Pheno_go_terms_Wild_RNA.txt
sed 's/a0000//g' Pheno_go_terms_Wild_RNA.txt | sed 's/Do1_0//g' | sed 's/_//' | sed 's/G0....//' | sed 's/ /\t/g'  > Phenogram_RNA_GOterms_Wild_w_c.txt
# add header and remove contigs
nano Phenogram_RNA_GOterms_Wild_w_c.txt

cp Phenogram_RNA_GOterms_Wild_w_c.txt Phenograms/
cd Phenograms/
ruby pheno_gram.rb -i Phenogram_RNA_GOterms_Wild_w_c.txt -g Dryas_Genome.txt -t "Dryas RNA Warming GO terms" -o RNA_Warming_GO -f jpg

# copy to local machine
scp -v celphin@cedar.computecanada.ca:/home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/Phenograms/RNA_Warming_GO.jpg .


#---------------------------
# DMR_RNA
# need to first extract the rows that contain RNA

grep "GO:" Gene_DMR_RNA_GO_Left_join_table.tsv | grep "Wild_W_C"  | awk '{print $1 " " $11 " " $14}'  | sort | uniq  > Pheno_go_terms_Wild.txt
sed 's/a0000//g' Pheno_go_terms_Wild.txt | sed 's/Do1_0//g' | sed 's/_//' | sed 's/G0....//' | sed 's/ /\t/g'  > Phenogram_GOterms_Wild_w_c_DMR-RNA.txt
# add header and remove contigs
nano Phenogram_GOterms_Wild_w_c_DMR-RNA.txt

cp Phenogram_GOterms_Wild_w_c_DMR-RNA.txt Phenograms/
cd Phenograms/
ruby pheno_gram.rb -i Phenogram_GOterms_Wild_w_c_DMR-RNA.txt -g Dryas_Genome.txt -t "Dryas DMR-RNA Warming GO terms" -o DMR_RNA_Warming_GO -f jpg

# copy to local machine
scp -v celphin@cedar.computecanada.ca:/home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/Phenograms/DMR_RNA_Warming_GO.jpg .

#---------------------------
# DMR Phenology
grep "GO:" Gene_DMR_Total_GO_Merged_table.tsv | grep "Mat_Sen"  | awk '{print $1 " " $11 " " $14}'  | sort | uniq  > Pheno_go_terms_Mat_Sen.txt
sed 's/a0000//g' Pheno_go_terms_Mat_Sen.txt | sed 's/Do1_0//g' | sed 's/_//' | sed 's/G0....//' | sed 's/ /\t/g'  > Phenogram_GOterms_Mat_Sen.txt
# add header and remove contigs
nano Phenogram_GOterms_Mat_Sen.txt

cp Phenogram_GOterms_Mat_Sen.txt Phenograms/
cd Phenograms/

ruby pheno_gram.rb -i Phenogram_GOterms_Mat_Sen.txt -g Dryas_Genome.txt -t "Dryas DMR Phenology GO terms" -o DMR_Mat_Sen_GO -f jpg

# copy to local machine
scp -v celphin@cedar.computecanada.ca:/home/celphin/scratch/Dryas/MS_Dryas_Merged_Data/Phenograms/DMR_Mat_Sen_GO.jpg .
